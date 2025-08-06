#include "hamiltonian.hpp"
#include "constants.hpp"
#include "chance_of_jump.hpp"
#include "define_params.hpp"
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <iostream>
#include <complex>
#include <cmath>

namespace hamiltonian {

using namespace Eigen;
using namespace std;
using namespace define_params;

tuple<int, int, int> from_number_to_tuple(int n, int n_x_max, int n_z_max) {
    int z = n % n_z_max;
    int xs = n / n_z_max;
    int x = xs % n_x_max;
    int s = xs / n_x_max;
    return {s, x, z};
}

MatrixXcd build_H_zero_freq(const MatrixXcd& H_ground,
                            const MatrixXcd& H_excited,
                            int n_ground_states,
                            int n_excited_states,
                            int n_x_max,
                            int n_z_max,
                            const Params& params) {

    int ng = n_ground_states;
    int ne = n_excited_states;
    MatrixXcd top(ng, ng + ne);
    MatrixXcd bottom(ne, ng + ne);
    top << H_ground, MatrixXcd::Zero(ng, ne);
    bottom << MatrixXcd::Zero(ne, ng), H_excited;
    MatrixXcd H_states(ng + ne, ng + ne);
    H_states << top, bottom;

    int size = (ng + ne) * n_x_max * n_z_max;
    MatrixXcd H0 = MatrixXcd::Zero(size, size);

    for (int i = 0; i < size; ++i) {
        auto [si, xi, zi] = from_number_to_tuple(i, n_x_max, n_z_max);
        for (int j = 0; j < size; ++j) {
            auto [sj, xj, zj] = from_number_to_tuple(j, n_x_max, n_z_max);
            if (xi == xj && zi == zj)
                H0(i, j) = H_states(si, sj);
        }
    }

    for (int i = 0; i < size; ++i) {
        auto [s, x, z] = from_number_to_tuple(i, n_x_max, n_z_max);
        double deltaE = (x * params.omega_x + z * params.omega_z);
        H0(i, i) += complex<double>(deltaE, 0);
    }

    return H0;
}

MatrixXcd build_H_light_transition(double Iii,
                                   const Vector3cd& pol,
                                   const Vector3d& vec,
                                   const vector<MatrixXd>& G,
                                   double Isat,
                                   double mass,
                                   double wavelength,
                                   int n_ground_states,
                                   int n_excited_states,
                                   int n_x_max,
                                   int n_z_max,
                                   const Params& params,
                                   bool excite) {
    int size = (n_ground_states + n_excited_states) * n_x_max * n_z_max;
    Vector3d wavevector = vec * 2.0 * parameters::pi / wavelength;
    double sat_param = sqrt(Iii / (2.0 * Isat));
    MatrixXcd Hex = MatrixXcd::Zero(size, size);

    for (int i = 0; i < size; ++i) {
        auto [si, xi, zi] = from_number_to_tuple(i, n_x_max, n_z_max);
        for (int j = 0; j < size; ++j) {
            auto [sj, xj, zj] = from_number_to_tuple(j, n_x_max, n_z_max);

            bool cond = excite ? (si >= n_ground_states && sj < n_ground_states)
                               : (sj >= n_ground_states && si < n_ground_states);
            if (!cond) continue;

            int e = excite ? si - n_ground_states : sj - n_ground_states;
            int g = excite ? sj : si;

            complex<double> contrib = 0.0;
            for (int pol_dir = 0; pol_dir < 3; ++pol_dir) {
                complex<double> Gijk = G[pol_dir](e, g);
                Vector3i from = excite ? Vector3i{xj, xi, zj} : Vector3i{xi, xj, zi};
                Vector3i to   = excite ? Vector3i{xi, xi, zi} : Vector3i{xj, xj, zj};
                chance_of_jump::Vec3 wavevector_ = { wavevector[0], wavevector[1], wavevector[2] };
                std::array<int, 3> from_ = { from[0], from[1], from[2] };
                std::array<int, 3> to_ = { to[0], to[1], to[2] };
                
                complex<double> ov = chance_of_jump::calculate_overlap_from_k(
                    wavevector_, mass,
                    params.omega_x, params.omega_z,
                    from_, to_
                );
                if (excite) ov = conj(ov);
                contrib += Gijk * sat_param * pol[pol_dir] * ov;
            }
            Hex(i, j) = contrib / Complex(2, 0);
        }
    }
    return Hex;
}

MatrixXcd build_H_light_transition_excite(double Iii,
                                          const Vector3cd& pol,
                                          const Vector3d& vec,
                                          const vector<MatrixXd>& G,
                                          double Isat,
                                          double mass,
                                          double wavelength,
                                          int n_ground_states,
                                          int n_excited_states,
                                          int n_x_max,
                                          int n_z_max,
                                          const Params& params) {
    return build_H_light_transition(Iii, pol, vec, G, Isat, mass, wavelength,
                                    n_ground_states, n_excited_states, n_x_max,
                                    n_z_max, params, true);
}

MatrixXcd build_H_light_transition_deexci(double Iii,
                                          const Vector3cd& pol,
                                          const Vector3d& vec,
                                          const vector<MatrixXd>& G,
                                          double Isat,
                                          double mass,
                                          double wavelength,
                                          int n_ground_states,
                                          int n_excited_states,
                                          int n_x_max,
                                          int n_z_max,
                                          const Params& params) {
    return build_H_light_transition(Iii, pol, vec, G, Isat, mass, wavelength,
                                    n_ground_states, n_excited_states, n_x_max,
                                    n_z_max, params, false);
}

} // namespace hamiltonian
