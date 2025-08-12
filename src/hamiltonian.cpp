#include "hamiltonian.hpp"
#include "constants.hpp"
#include "chance_of_jump.hpp"

namespace hamiltonian {
tuple<int, int, int> from_number_to_tuple(int n, int n_x_max, int n_z_max) {
    int z = n % n_z_max;
    int xs = n / n_z_max;
    int x = xs % n_x_max;
    int s = xs / n_x_max;
    return {s, x, z};
}

States diagonalize_hamiltonian(const States& states) {
    States new_states = states;
    auto H_ground = states.H_ground;
    auto H_excited = states.H_excited;
    auto G = states.G;
    // G is [3][n_excited_states][n_ground_states]
    // Diagonalize the ground and excited Hamiltonians, and modify G accordingly
    Eigen::SelfAdjointEigenSolver<ComplexMat> solver_ground(H_ground);
    if (solver_ground.info() != Eigen::Success) {
        throw std::runtime_error("Diagonalization of ground Hamiltonian failed.");
    }
    new_states.H_ground = solver_ground.eigenvalues().asDiagonal();
    new_states.G[0] = (new_states.G[0] * solver_ground.eigenvectors()).cwiseAbs();
    new_states.G[1] = (new_states.G[1] * solver_ground.eigenvectors()).cwiseAbs();
    new_states.G[2] = (new_states.G[2] * solver_ground.eigenvectors()).cwiseAbs();
    Eigen::SelfAdjointEigenSolver<ComplexMat> solver_excited(H_excited);
    if (solver_excited.info() != Eigen::Success) {
        throw std::runtime_error("Diagonalization of excited Hamiltonian failed.");
    }
    new_states.H_excited = solver_excited.eigenvalues().asDiagonal();
    new_states.G[0] = (solver_excited.eigenvectors().adjoint() * new_states.G[0]).cwiseAbs();
    new_states.G[1] = (solver_excited.eigenvectors().adjoint() * new_states.G[1]).cwiseAbs();
    new_states.G[2] = (solver_excited.eigenvectors().adjoint() * new_states.G[2]).cwiseAbs();
    return new_states;
}

ComplexVec define_partition_hamiltonian(const ComplexMat& H, int n_x_max, int n_z_max, const Params& params) {
    int n_states = H.rows();
    ComplexVec hamiltonian = ComplexVec::Zero(n_states * n_x_max * n_z_max);
    for (int i = 0; i < n_states; ++i) {
        for (int j = 0; j < n_x_max; ++j) {
            for (int k = 0; k < n_z_max; ++k) {
                hamiltonian(i * n_x_max * n_z_max + j * n_z_max + k) = H(i, i) + params.omega_x * j + params.omega_z * k;
            }
        }
    }
    return hamiltonian;
}

ComplexMat define_V_plus(const States& states,
                         const Params& params,
                         int f,
                         int l,
                         int I_index,
                         int D_index)
{
    double Isat = parameters::pi * parameters::plancks_constant * parameters::speed_of_light * states.G_tot /
                  (3.0 * pow(states.transition_lambda, 3));
    int size_ground = states.n_ground_states * params.n_x_max * params.n_z_max;
    int size_excited = states.n_excited_states * params.n_x_max * params.n_z_max;
    Vector3d wavevector = params.k.row(f) * 2.0 * parameters::pi / states.transition_lambda;
    double sat_param = sqrt(params.I(I_index, f) / (2.0 * Isat));
    MatrixXcd V_plus(size_excited, size_ground);

    for (int i = 0; i < size_excited; ++i) {
        auto [e, xi, zi] = from_number_to_tuple(i, params.n_x_max, params.n_z_max);
        auto [g, xj, zj] = from_number_to_tuple(l, params.n_x_max, params.n_z_max);
        complex<double> contrib = 0.0;
        for (int pol_dir = 0; pol_dir < 3; ++pol_dir) {
            complex<double> Gijk = states.G[pol_dir](e, g);
            Vector3i from = Vector3i{xj, xi, zj};
            Vector3i to   = Vector3i{xi, xi, zi};
            chance_of_jump::Vec3 wavevector_ = { wavevector[0], wavevector[1], wavevector[2] };
            std::array<int, 3> from_ = { from[0], from[1], from[2] };
            std::array<int, 3> to_ = { to[0], to[1], to[2] };
            
            complex<double> ov = chance_of_jump::calculate_overlap_from_k(
                wavevector_, params.mass,
                params.omega_x, params.omega_z,
                from_, to_
            );
            ov = conj(ov);
            contrib += Gijk * sat_param * params.s.row(f)[pol_dir] * ov;
        }
        V_plus(i, l) = contrib / Complex(2, 0);
    }
    return V_plus;
}

ComplexMat define_V_minus(const States& states,
                          const Params& params,
                          int f,
                          int I_index,
                          int D_index)
{
    double Isat = parameters::pi * parameters::plancks_constant * parameters::speed_of_light * states.G_tot /
                  (3.0 * pow(states.transition_lambda, 3));
    int size_ground = states.n_ground_states * params.n_x_max * params.n_z_max;
    int size_excited = states.n_excited_states * params.n_x_max * params.n_z_max;
    Vector3d wavevector = params.k.row(f) * 2.0 * parameters::pi / states.transition_lambda;
    double sat_param = sqrt(params.I(I_index, f) / (2.0 * Isat));
    MatrixXcd V_minus(size_ground, size_excited);
    
    for (int i = 0; i < size_ground; ++i) {
        auto [g, xi, zi] = from_number_to_tuple(i, params.n_x_max, params.n_z_max);
        for (int j = 0; j < size_excited; ++j) {
            auto [e, xj, zj] = from_number_to_tuple(j, params.n_x_max, params.n_z_max);
            complex<double> contrib = 0.0;
            for (int pol_dir = 0; pol_dir < 3; ++pol_dir) {
                complex<double> Gijk = states.G[pol_dir](e, g);
                Vector3i from = Vector3i{xi, xj, zi};
                Vector3i to   = Vector3i{xj, xj, zj};
                chance_of_jump::Vec3 wavevector_ = { wavevector[0], wavevector[1], wavevector[2] };
                std::array<int, 3> from_ = { from[0], from[1], from[2] };
                std::array<int, 3> to_ = { to[0], to[1], to[2] };
                
                complex<double> ov = chance_of_jump::calculate_overlap_from_k(
                    wavevector_, params.mass,
                    params.omega_x, params.omega_z,
                    from_, to_
                );
                contrib += Gijk * sat_param * params.s.row(f)[pol_dir] * ov;
            }
            V_minus(i, j) = contrib / Complex(2, 0);
        }
    }
    return V_minus;
}

ComplexMat V_minus(const States& states,
                   const Params& params,
                   int I_index,
                   int D_index)
{
    int n_excited_states = states.n_excited_states;
    int n_ground_states = states.n_ground_states;
    ComplexMat V_minus_total = ComplexMat::Zero(n_ground_states * params.n_x_max * params.n_z_max,
                                                n_excited_states * params.n_x_max * params.n_z_max);
    
    for (int f = 0; f < params.n_beams; ++f) {
        ComplexMat V_minus_f = define_V_minus(states, params, f, I_index, D_index);
        V_minus_total += V_minus_f;
    }
    
    return V_minus_total;
}

ComplexMat build_W(const States& states,
                   const Params& params,
                   int I_index,
                   int D_index)
{
    // W = \sum_{f=0}^{n_beams-1} \sum_{l=0}^{n_ground_states * params.n_x_max * params.n_z_max - 1}
    // H_NH(f, l)^{-1} * V_plus(f, l)
    // Where H_NH(f, l)^{-1}(e, e) = 1 / (H_excited(e, e) - (i/2) G_tot - H_ground(l, l) - Delta(f))
    int n_excited_states = states.n_excited_states;
    int n_ground_states = states.n_ground_states;
    ComplexMat W = ComplexMat::Zero(n_excited_states * params.n_x_max * params.n_z_max,
                                    n_ground_states * params.n_x_max * params.n_z_max);
    ComplexVec H_ground_diag = define_partition_hamiltonian(states.H_ground, params.n_x_max, params.n_z_max, params);
    ComplexVec H_excited_diag = define_partition_hamiltonian(states.H_excited, params.n_x_max, params.n_z_max, params);
    for (int f = 0; f < params.n_beams; ++f) {
        for (int l = 0; l < n_ground_states * params.n_x_max * params.n_z_max; ++l) {
            ComplexMat V_plus_fl = define_V_plus(states, params, f, l, I_index, D_index);
            ComplexMat H_NH(n_excited_states * params.n_x_max * params.n_z_max,
                            n_excited_states * params.n_x_max * params.n_z_max);
            for (int e = 0; e < n_excited_states * params.n_x_max * params.n_z_max; ++e) {
                    H_NH(e, e) = Complex(1,0) / (H_excited_diag(e) - Complex(0, 0.5 * states.G_tot) -
                                            H_ground_diag(l) - params.D(D_index, f)); // blue detuning
            }
            W += H_NH * V_plus_fl;
        }
    }
    return W;
}

std::vector<std::tuple<int, int, double>> build_L(
    const std::vector<Eigen::MatrixXd>& G,
    int n_x, int n_z,
    int num_nonzero,
    double mass,
    double omega_x,
    double omega_z,
    double wavelength,
    const std::array<double, 3>& B_direction,
    unsigned int seed) 
{
    const int np = static_cast<int>(G.size());
    const int ne = static_cast<int>(G[0].rows());
    const int ng = static_cast<int>(G[0].cols());

    int full_size = num_nonzero * (n_x * n_x * n_x) * (n_z * n_z);
    int full_full_size = np * ne * ng * (n_x * n_x * n_x) * (n_z * n_z);

    std::vector<std::tuple<int, int, double>> Lt;

    auto from_tuple_to_number = [&](int s, int x, int z) {
        return (s * n_x + x) * n_z + z;
    };

    auto bound = [&](int x) {
        return std::clamp(x, 0, n_x - 1);
    };

    for (int i = 0; i < full_full_size; ++i) {
        // Unpack i to np_, ne_, ng_, nxf, nyf, nzf, nxi, nzi
        int rem = i;
        int nzi = rem % n_z; rem /= n_z;
        int nxi = rem % n_x; rem /= n_x;
        int nzf = rem % n_z; rem /= n_z;
        int nyf = rem % n_x; rem /= n_x;
        int nxf = rem % n_x; rem /= n_x;
        int ng_ = rem % ng;  rem /= ng;
        int ne_ = rem % ne;  rem /= ne;
        int np_ = rem % np;

        if (G[np_](ne_, ng_) == 0.0) continue;

        int excited_state = from_tuple_to_number(ne_, nxi, nzi);
        int ground_state = from_tuple_to_number(ng_, bound(nxf + nyf - nxi), nzf);

        std::array<int, 3> start = {nxi, nxi, nzi};
        std::array<int, 3> end   = {nxf, nyf, nzf};

        std::vector<double> chance = chance_of_jump::calculate_chance_of_jump(
            seed + i, // ensure variability if needed
            np_, mass, omega_x, omega_z, wavelength,
            B_direction, start, end, 500
        );

        double weighted_sum = 0.0;
        for (double val : chance) {
            weighted_sum += std::abs(val);
        }

        double average = weighted_sum / static_cast<double>(chance.size());
        if (average <= 1e-3) continue;
        average *= G[np_](ne_, ng_);
        double weight_of_jump = std::sqrt(average);

        Lt.emplace_back(ground_state, excited_state, weight_of_jump);
    }

    return Lt;
}

ComplexMat build_H(const States& states,
                   const Params& params,
                   int I_index,
                   int D_index,
                   ComplexMat& W)
{
    ComplexVec H_ground_diag = define_partition_hamiltonian(states.H_ground, params.n_x_max, params.n_z_max, params);
    // H_eff = H_ground_diag - 0.5 * [V_minus * W + h.c.]
    ComplexMat H_eff = ComplexMat::Zero(states.n_ground_states * params.n_x_max * params.n_z_max,
                                        states.n_ground_states * params.n_x_max * params.n_z_max);
    ComplexMat V_minus_total = V_minus(states, params, I_index, D_index);
    // first add H_ground_diag
    for (int i = 0; i < states.n_ground_states * params.n_x_max * params.n_z_max; ++i) {
        H_eff(i, i) = H_ground_diag(i);
    }
    ComplexMat V_minus_W = V_minus_total * W;
    // now add -0.5 * (V_minus * W + h.c.)
    H_eff -= 0.5 * (V_minus_W + V_minus_W.adjoint());
    return H_eff;
}
} // namespace hamiltonian
