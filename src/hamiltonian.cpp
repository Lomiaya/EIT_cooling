#include "hamiltonian.hpp"
#include "constants.hpp"
#include "chance_of_jump.hpp"

#include <iostream>

namespace hamiltonian {
tuple<int, int, int> from_number_to_tuple(int n, int n_x_max, int n_z_max) {
    int z = n % n_z_max;
    int xs = n / n_z_max;
    int x = xs % n_x_max;
    int s = xs / n_x_max;
    return {s, x, z};
}
tuple<int, int, int, int> from_number_to_tuple(int n, int n_x_max, int n_y_max, int n_z_max) {
    int z = n % n_z_max;
    int yx = n / n_z_max;
    int y = yx % n_y_max;
    int xs = yx / n_y_max;
    int x = xs % n_x_max;
    int s = xs / n_x_max;
    return {s, x, y, z};
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

/**
 * @brief Add Stark shift to Hamiltonian
 */
DoubleVec define_partition_hamiltonian(const ComplexMat& H, const ComplexMat& H_stark, int n_x_max, int n_z_max, const Params& params) {
    int n_states = H.rows();
    DoubleVec hamiltonian = DoubleVec::Zero(n_states * n_x_max * n_z_max);
    for (int i = 0; i < n_states; ++i) {
        for (int j = 0; j < n_x_max; ++j) {
            for (int k = 0; k < n_z_max; ++k) {
                hamiltonian(i * n_x_max * n_z_max + j * n_z_max + k) = H(i, i).real() + (params.omega_x * j + params.omega_z * k) * sqrt(H_stark(i, i).real()) + (- params.trap_depth) * H_stark(i, i).real();
            }
        }
    }
    return hamiltonian;
}
DoubleVec define_partition_hamiltonian(const ComplexMat& H, const ComplexMat& H_stark, int n_x_max, int n_y_max, int n_z_max, const Params& params) {
    int n_states = H.rows();
    DoubleVec hamiltonian = DoubleVec::Zero(n_states * n_x_max * n_y_max * n_z_max);
    for (int i = 0; i < n_states; ++i) {
        for (int j = 0; j < n_x_max; ++j) {
            for (int k = 0; k < n_y_max; ++k) {
                for (int l = 0; l < n_z_max; ++l) {
                    hamiltonian(i * n_x_max * n_y_max * n_z_max + j * n_y_max * n_z_max + k * n_z_max + l) = H(i, i).real() + (params.omega_x * j + params.omega_y * k + params.omega_z * l) * sqrt(H_stark(i, i).real()) + (- params.trap_depth) * H_stark(i, i).real();
                }
            }
        }
    }
    return hamiltonian;
}
SpectrumMatrix multiply(const SpectrumMatrix& A,
                        const SpectrumMatrix& B)
{
    SpectrumMatrix result;
    for (const auto& [A_mat, A_freq] : A) {
        for (const auto& [B_mat, B_freq] : B) {
            result.emplace_back(A_mat * B_mat, A_freq + B_freq);
        }
    }
    return result;
}
SpectrumMatrix multiply(const SpectrumMatrix& A,
                        const SpectrumMatrix& B,
                        double threshold)
{
    SpectrumMatrix result;
    for (const auto& [A_mat, A_freq] : A) {
        for (const auto& [B_mat, B_freq] : B) {
            if (std::abs(A_freq + B_freq) < threshold) {
                bool found = false;
                for (auto& [res_mat, res_freq] : result) {
                    if (std::abs(B_freq + A_freq - res_freq) < 1) {
                        res_mat += A_mat * B_mat;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    result.emplace_back(A_mat * B_mat, A_freq + B_freq);
                }
            }
        }
    }
    return result;
}
SpectrumMatrix multiply(const ComplexMat& A,
                        const SpectrumMatrix& B)
{
    SpectrumMatrix result;
    for (const auto& [B_mat, B_freq] : B) {
        result.emplace_back(A * B_mat, B_freq);
    }
    return result;
}
SpectrumMatrix multiply(const SpectrumMatrix& A,
                        const ComplexMat& B)
{
    SpectrumMatrix result;
    for (const auto& [A_mat, A_freq] : A) {
        result.emplace_back(A_mat * B, A_freq);
    }
    return result;
}
SpectrumMatrix multiply(double c, const SpectrumMatrix& A)
{
    SpectrumMatrix result;
    for (const auto& [A_mat, A_freq] : A) {
        result.emplace_back(c * A_mat, A_freq);
    }
    return result;
}
SpectrumMatrix multiply(Complex c, const SpectrumMatrix& A)
{
    SpectrumMatrix result;
    for (const auto& [A_mat, A_freq] : A) {
        result.emplace_back(c * A_mat, A_freq);
    }
    return result;
}
SpectrumMatrix adjoint(const SpectrumMatrix& A)
{
    SpectrumMatrix result;
    for (const auto& [A_mat, A_freq] : A) {
        result.emplace_back(A_mat.adjoint(), -A_freq); // Frequency is negated for adjoint operation
    }
    return result;
}
SpectrumMatrix addition(const SpectrumMatrix& A,
                        const SpectrumMatrix& B)
{
    SpectrumMatrix result = A;
    for (const auto& [B_mat, B_freq] : B) {
        bool found = false;
        for (auto& [res_mat, res_freq] : result) {
            if (std::abs(res_freq - B_freq) < 1) {
                res_mat += B_mat;
                found = true;
                break;
            }
        }
        if (!found) {
            result.emplace_back(B_mat, B_freq);
        }
    }
    return result;
}
ComplexMat evaluate(const SpectrumMatrix& A, double t)
{
    ComplexMat result = ComplexMat::Zero(A[0].first.rows(), A[0].first.cols());
    for (const auto& [A_mat, A_freq] : A) {
        result += A_mat * std::exp(Complex(0, A_freq * t));
    }
    return result;
}
SpectrumMatrix low_pass_filter(const SpectrumMatrix& A, double threshold)
{
    SpectrumMatrix result;
    for (const auto& [A_mat, A_freq] : A) {
        if (abs(A_freq) < threshold) {
            result.emplace_back(A_mat, A_freq);
        }
    }
    return result;
}
SpectrumMatrix cleanup(const SpectrumMatrix& A)
{
    auto result = addition({}, A);
    // Remove near-zero matrices
    result.erase(std::remove_if(result.begin(), result.end(),
                                [](const auto& pair) { return pair.first.norm() < 1e-12; }),
                 result.end());
    return result;
}

SpectrumMatrix define_V_plus_2d(const States& states,
                                const Params& params,
                                int f,
                                int l,
                                int I_index,
                                int D_index,
                                DoubleVec& H_ground_diag)
{
    double Isat = parameters::pi * parameters::plancks_constant * parameters::speed_of_light * states.G_tot /
                  (3.0 * pow(states.transition_lambda, 3));
    int size_ground = states.n_ground_states * params.n_x_max * params.n_z_max;
    int size_excited = states.n_excited_states * params.n_x_max * params.n_z_max;
    Vector3d wavevector = params.k.row(f) * 2.0 * parameters::pi / states.transition_lambda;
    double sat_param = sqrt(params.I(I_index, f) / (2.0 * Isat));
    MatrixXcd V_plus = MatrixXcd::Zero(size_excited, size_ground);

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
    double freq = - H_ground_diag(l) - params.D(D_index, f); // blue detuning
    return {{V_plus, freq}};
}
SpectrumMatrix define_V_plus_3d(const States& states,
                                const Params& params,
                                int f,
                                int l,
                                int I_index,
                                int D_index,
                                DoubleVec& H_ground_diag)
{
    double Isat = parameters::pi * parameters::plancks_constant * parameters::speed_of_light * states.G_tot /
                (3.0 * pow(states.transition_lambda, 3));
    int size_ground = states.n_ground_states * params.n_x_max * params.n_y_max * params.n_z_max;
    int size_excited = states.n_excited_states * params.n_x_max * params.n_y_max * params.n_z_max;
    Vector3d wavevector = params.k.row(f) * 2.0 * parameters::pi / states.transition_lambda;
    double sat_param = sqrt(params.I(I_index, f) / (2.0 * Isat));
    MatrixXcd V_plus = MatrixXcd::Zero(size_excited, size_ground);

    for (int i = 0; i < size_excited; ++i) {
        auto [e, xi, yi, zi] = from_number_to_tuple(i, params.n_x_max, params.n_y_max, params.n_z_max);
        auto [g, xj, yj, zj] = from_number_to_tuple(l, params.n_x_max, params.n_y_max, params.n_z_max);
        complex<double> contrib = 0.0;
        for (int pol_dir = 0; pol_dir < 3; ++pol_dir) {
            complex<double> Gijk = states.G[pol_dir](e, g);
            Vector3i from = Vector3i{xj, yj, zj};
            Vector3i to   = Vector3i{xi, yi, zi};
            chance_of_jump::Vec3 wavevector_ = { wavevector[0], wavevector[1], wavevector[2] };
            std::array<int, 3> from_ = { from[0], from[1], from[2] };
            std::array<int, 3> to_ = { to[0], to[1], to[2] };
            
            complex<double> ov = chance_of_jump::calculate_overlap_from_k(
                wavevector_, params.mass,
                params.omega_x, params.omega_y, params.omega_z,
                from_, to_
            );
            ov = conj(ov);
            contrib += Gijk * sat_param * params.s.row(f)[pol_dir] * ov;
        }
        V_plus(i, l) = contrib / Complex(2, 0);
    }
    double freq = - H_ground_diag(l) - params.D(D_index, f); // blue detuning
    return {{V_plus, freq}};
}
SpectrumMatrix define_V_plus(const States& states,
                             const Params& params,
                             int f,
                             int l,
                             int I_index,
                             int D_index,
                             DoubleVec& H_ground_diag)
{
    if (params.do_2d_sim) {
        return define_V_plus_2d(states, params, f, l, I_index, D_index, H_ground_diag);
    } else {
        return define_V_plus_3d(states, params, f, l, I_index, D_index, H_ground_diag);
    }
}

SpectrumMatrix define_V_minus_2d(const States& states,
                                 const Params& params,
                                 int f,
                                 int I_index,
                                 int D_index,
                                 DoubleVec& H_ground_diag)
{
    double Isat = parameters::pi * parameters::plancks_constant * parameters::speed_of_light * states.G_tot /
                  (3.0 * pow(states.transition_lambda, 3));
    int size_ground = states.n_ground_states * params.n_x_max * params.n_z_max;
    int size_excited = states.n_excited_states * params.n_x_max * params.n_z_max;
    Vector3d wavevector = params.k.row(f) * 2.0 * parameters::pi / states.transition_lambda;
    double sat_param = sqrt(params.I(I_index, f) / (2.0 * Isat));
    SpectrumMatrix V_minus;
    
    for (int i = 0; i < size_ground; ++i) {
        auto [g, xi, zi] = from_number_to_tuple(i, params.n_x_max, params.n_z_max);
        MatrixXcd V_minus_this = MatrixXcd::Zero(size_ground, size_excited);
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
            V_minus_this(i, j) = contrib / Complex(2, 0);
        }
        double freq = params.D(D_index, f) + H_ground_diag(i); // blue detuning
        V_minus.push_back({V_minus_this, freq});
    }
    return V_minus;
}
SpectrumMatrix define_V_minus_3d(const States& states,
                                 const Params& params,
                                 int f,
                                 int I_index,
                                 int D_index,
                                 DoubleVec& H_ground_diag)
{
    double Isat = parameters::pi * parameters::plancks_constant * parameters::speed_of_light * states.G_tot /
                (3.0 * pow(states.transition_lambda, 3));
    int size_ground = states.n_ground_states * params.n_x_max * params.n_y_max * params.n_z_max;
    int size_excited = states.n_excited_states * params.n_x_max * params.n_y_max * params.n_z_max;
    Vector3d wavevector = params.k.row(f) * 2.0 * parameters::pi / states.transition_lambda;
    double sat_param = sqrt(params.I(I_index, f) / (2.0 * Isat));
    SpectrumMatrix V_minus;
    
    for (int i = 0; i < size_ground; ++i) {
        auto [g, xi, yi, zi] = from_number_to_tuple(i, params.n_x_max, params.n_y_max, params.n_z_max);
        MatrixXcd V_minus_this = MatrixXcd::Zero(size_ground, size_excited);
        for (int j = 0; j < size_excited; ++j) {
            auto [e, xj, yj, zj] = from_number_to_tuple(j, params.n_x_max, params.n_y_max, params.n_z_max);
            complex<double> contrib = 0.0;
            for (int pol_dir = 0; pol_dir < 3; ++pol_dir) {
                complex<double> Gijk = states.G[pol_dir](e, g);
                Vector3i from = Vector3i{xi, yi, zi};
                Vector3i to   = Vector3i{xj, yj, zj};
                chance_of_jump::Vec3 wavevector_ = { wavevector[0], wavevector[1], wavevector[2] };
                std::array<int, 3> from_ = { from[0], from[1], from[2] };
                std::array<int, 3> to_ = { to[0], to[1], to[2] };
                
                complex<double> ov = chance_of_jump::calculate_overlap_from_k(
                    wavevector_, params.mass,
                    params.omega_x, params.omega_y, params.omega_z,
                    from_, to_
                );
                contrib += Gijk * sat_param * params.s.row(f)[pol_dir] * ov;
            }
            V_minus_this(i, j) = contrib / Complex(2, 0);
        }
        double freq = params.D(D_index, f) + H_ground_diag(i); // blue detuning
        V_minus.push_back({V_minus_this, freq});
    }
    return V_minus;
}
SpectrumMatrix define_V_minus(const States& states,
                              const Params& params,
                              int f,
                              int I_index,
                              int D_index,
                              DoubleVec& H_ground_diag)
{
    if (params.do_2d_sim) {
        return define_V_minus_2d(states, params, f, I_index, D_index, H_ground_diag);
    } else {
        return define_V_minus_3d(states, params, f, I_index, D_index, H_ground_diag);
    }
}

SpectrumMatrix V_minus(const States& states,
                       const Params& params,
                       int I_index,
                       int D_index,
                       DoubleVec& H_ground_diag)
{
    int n_excited_states = states.n_excited_states;
    int n_ground_states = states.n_ground_states;
    SpectrumMatrix V_minus_total;
    
    for (int f = 0; f < params.n_beams; ++f) {
        SpectrumMatrix V_minus_f = define_V_minus(states, params, f, I_index, D_index, H_ground_diag);
        V_minus_total = addition(V_minus_total, V_minus_f);
    }
    
    return V_minus_total;
}


SpectrumMatrix build_W_2d(const States& states,
                       const Params& params,
                       int I_index,
                       int D_index)
{
    int n_excited_states = states.n_excited_states;
    int n_ground_states = states.n_ground_states;
    SpectrumMatrix W;
    DoubleVec H_ground_diag = define_partition_hamiltonian(states.H_ground, states.H_ground_stark, params.n_x_max, params.n_z_max, params);
    DoubleVec H_excited_diag = define_partition_hamiltonian(states.H_excited, states.H_excited_stark, params.n_x_max, params.n_z_max, params);
    for (int f = 0; f < params.n_beams; ++f) {
        for (int l = 0; l < n_ground_states * params.n_x_max * params.n_z_max; ++l) {
            SpectrumMatrix V_plus_fl = define_V_plus(states, params, f, l, I_index, D_index, H_ground_diag);
            ComplexMat H_NH = ComplexMat::Zero(n_excited_states * params.n_x_max * params.n_z_max,
                                               n_excited_states * params.n_x_max * params.n_z_max);
            for (int e = 0; e < n_excited_states * params.n_x_max * params.n_z_max; ++e) {
                    H_NH(e, e) = Complex(1,0) / (H_excited_diag(e) - Complex(0, 0.5 * states.G_tot) -
                                                 H_ground_diag(l) - params.D(D_index, f)); // blue detuning
            }
            W = addition(W, multiply(H_NH, V_plus_fl));
        }
    }
    return W;
}
SpectrumMatrix build_W_3d(const States& states,
                       const Params& params,
                       int I_index,
                       int D_index)
{
    int n_excited_states = states.n_excited_states;
    int n_ground_states = states.n_ground_states;
    SpectrumMatrix W;
    DoubleVec H_ground_diag = define_partition_hamiltonian(states.H_ground, states.H_ground_stark, params.n_x_max, params.n_y_max, params.n_z_max, params);
    DoubleVec H_excited_diag = define_partition_hamiltonian(states.H_excited, states.H_excited_stark, params.n_x_max, params.n_y_max, params.n_z_max, params);
    for (int f = 0; f < params.n_beams; ++f) {
        for (int l = 0; l < n_ground_states * params.n_x_max * params.n_y_max * params.n_z_max; ++l) {
            SpectrumMatrix V_plus_fl = define_V_plus(states, params, f, l, I_index, D_index, H_ground_diag);
            ComplexMat H_NH = ComplexMat::Zero(n_excited_states * params.n_x_max * params.n_y_max * params.n_z_max,
                                               n_excited_states * params.n_x_max * params.n_y_max * params.n_z_max);
            for (int e = 0; e < n_excited_states * params.n_x_max * params.n_y_max * params.n_z_max; ++e) {
                    H_NH(e, e) = Complex(1,0) / (H_excited_diag(e) - Complex(0, 0.5 * states.G_tot) -
                                                 H_ground_diag(l) - params.D(D_index, f)); // blue detuning
            }
            W = addition(W, multiply(H_NH, V_plus_fl));
        }
    }
    return W;
}
/**
 * @brief Build auxiliary term W = \sum_{f=0}^{n_beams-1} \sum_{l=0}^{n_ground_states * params.n_x_max * params.n_z_max - 1}
 * H_NH(f, l)^{-1} * V_plus(f, l)
 * Where H_NH(f, l)^{-1}(e, e) = 1 / (H_excited(e, e) - (i/2) G_tot - H_ground(l, l) - Delta(f))
 */
SpectrumMatrix build_W(const States& states,
                       const Params& params,
                       int I_index,
                       int D_index)
{
    // W = \sum_{f=0}^{n_beams-1} \sum_{l=0}^{n_ground_states * params.n_x_max * params.n_z_max - 1}
    // H_NH(f, l)^{-1} * V_plus(f, l)
    // Where H_NH(f, l)^{-1}(e, e) = 1 / (H_excited(e, e) - (i/2) G_tot - H_ground(l, l) - Delta(f))
    if (params.do_2d_sim) {
        return build_W_2d(states, params, I_index, D_index);
    } else {
        return build_W_3d(states, params, I_index, D_index);
    }
}

/**
 * @brief Build relaxation operator
 * 
 * Returns (decay to which ground_state, from which excited_state, martrix element)
 */
std::vector<std::tuple<int, int, double>> build_L_2d(
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
std::vector<std::tuple<int, int, double>> build_L_3d(
    const std::vector<Eigen::MatrixXd>& G,
    int n_x, int n_y, int n_z,
    int num_nonzero,
    double mass,
    double omega_x,
    double omega_y,
    double omega_z,
    double wavelength,
    const std::array<double, 3>& B_direction,
    unsigned int seed) 
{
    const int np = static_cast<int>(G.size());
    const int ne = static_cast<int>(G[0].rows());
    const int ng = static_cast<int>(G[0].cols());

    int full_size = num_nonzero * (n_x * n_x) * (n_y * n_y) * (n_z * n_z);
    int full_full_size = np * ne * ng * (n_x * n_x) * (n_y * n_y) * (n_z * n_z);

    std::vector<std::tuple<int, int, double>> Lt;

    auto from_tuple_to_number = [&](int s, int x, int y, int z) {
        return ((s * n_x + x) * n_y + y) * n_z + z;
    };

    auto bound = [&](int x) {
        return std::clamp(x, 0, n_x - 1);
    };

    for (int i = 0; i < full_full_size; ++i) {
        // Unpack i to np_, ne_, ng_, nxf, nyf, nzf, nxi, nyi, nzi
        int rem = i;
        int nzi = rem % n_z; rem /= n_z;
        int nyi = rem % n_y; rem /= n_y;
        int nxi = rem % n_x; rem /= n_x;
        int nzf = rem % n_z; rem /= n_z;
        int nyf = rem % n_y; rem /= n_y;
        int nxf = rem % n_x; rem /= n_x;
        int ng_ = rem % ng;  rem /= ng;
        int ne_ = rem % ne;  rem /= ne;
        int np_ = rem % np;

        if (G[np_](ne_, ng_) == 0.0) continue;

        int excited_state = from_tuple_to_number(ne_, nxi, nyi, nzi);
        int ground_state = from_tuple_to_number(ng_, nxf, nyf, nzf);

        std::array<int, 3> start = {nxi, nyi, nzi};
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

/**
 * @brief Build Hamiltonian.
 */
SpectrumMatrix build_H(const States& states,
                       const Params& params,
                       int I_index,
                       int D_index,
                       SpectrumMatrix& W,
                       double threshold)
{
    // H_eff = - V_minus * W. Memory-inefficient, but it works so I don't care.
    DoubleVec H_ground_diag;
    if (params.do_2d_sim)
        H_ground_diag = define_partition_hamiltonian(states.H_ground, states.H_ground_stark, params.n_x_max, params.n_z_max, params);
    else
        H_ground_diag = define_partition_hamiltonian(states.H_ground, states.H_ground_stark, params.n_x_max, params.n_y_max, params.n_z_max, params);
    SpectrumMatrix V_minus_total = V_minus(states, params, I_index, D_index, H_ground_diag);
    SpectrumMatrix V_minus_W = multiply(V_minus_total, W, threshold);
    // SpectrumMatrix H_eff = multiply(-0.5, (addition(V_minus_W, adjoint(V_minus_W))));
    SpectrumMatrix H_eff = multiply(-1, V_minus_W);
    return H_eff;
}
} // namespace hamiltonian
