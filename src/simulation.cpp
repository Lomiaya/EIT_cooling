#include <iostream>
#include <vector>
#include <complex>
#include <chrono>
#include <Eigen/Dense>
#include <map>
#include <tuple>
#include "constants.hpp"
#include "hamiltonian.hpp"
#include "simulation.hpp"
#include "define_params.hpp"
#include "stochastic_schrodinger_spin_motion.hpp"

namespace simulation {

using namespace std;
using namespace Eigen;
using namespace ss_spin;
using namespace define_params;

int from_tuple_to_number(int s, int x, int z, const Params& p) {
    return (s * p.n_x_max + x) * p.n_z_max + z;
}

tuple<int, int, int> from_number_to_tuple(int n, const Params& p) {
    int z = n % p.n_z_max;
    int xs = n / p.n_z_max;
    int x = xs % p.n_x_max;
    int s = xs / p.n_x_max;
    return make_tuple(s, x, z);
}

std::tuple<MatrixXd, MatrixXd, MatrixXd> simulate(const Params& params, const States& def_states, const int N, const double t_0, const int num_keys) {
    double Isat = parameters::pi * parameters::plancks_constant * parameters::speed_of_light * def_states.G_tot /
                  (3.0 * pow(def_states.transition_lambda, 3));
    cout << "Isat: " << Isat << endl;

    const auto& Ds = params.D;
    const auto& Is = params.I;
    const auto& polarizations = params.s;
    const auto& k_vectors = params.k;
    const int n_beams = params.n_beams;

    MatrixXd tot_jumps(Ds.rows(), Is.rows());
    MatrixXd avg_tempsx(Ds.rows(), Is.rows());
    MatrixXd avg_tempsz(Ds.rows(), Is.rows());

    for (size_t D_index = 0; D_index < Ds.rows(); ++D_index) {
        for (size_t I_index = 0; I_index < Is.rows(); ++I_index) {

            cout << D_index << ", " << I_index << endl;

            auto t0 = chrono::high_resolution_clock::now();

            int n_states = def_states.n_ground_states + def_states.n_excited_states;
            int size = n_states * params.n_x_max * params.n_z_max;

            MatrixXcd H0 = hamiltonian::build_H_zero_freq(def_states.H_ground, def_states.H_excited,
                def_states.n_ground_states, def_states.n_excited_states,
                params.n_x_max, params.n_z_max,
                params);

            vector<pair<MatrixXcd, double>> H = {{H0, 0.0}};

            for (int beam_n = 0; beam_n < n_beams; ++beam_n) {
                auto H_excite = hamiltonian::build_H_light_transition_excite(
                    Is(I_index,beam_n),
                    polarizations.row(beam_n),
                    k_vectors.row(beam_n),
                    def_states.G,
                    Isat,
                    params.mass,
                    def_states.transition_lambda,
                    def_states.n_ground_states,
                    def_states.n_excited_states,
                    params.n_x_max,
                    params.n_z_max,
                    params);

                auto H_deexci = hamiltonian::build_H_light_transition_deexci(
                    Is(I_index,beam_n),
                    polarizations.row(beam_n),
                    k_vectors.row(beam_n),
                    def_states.G,
                    Isat,
                    params.mass,
                    def_states.transition_lambda,
                    def_states.n_ground_states,
                    def_states.n_excited_states,
                    params.n_x_max,
                    params.n_z_max,
                    params);

                H.emplace_back(H_excite, -Ds(D_index,beam_n));
                H.emplace_back(H_deexci, Ds(D_index,beam_n));
            }

            cout << "1. Finished Building H!" << endl;

            unsigned int seed = 42;
            int num_G_nonzero_entries = 0;
            for (const auto& g : def_states.G) num_G_nonzero_entries += (g.array() != 0.0).count();

            auto L = ss_spin::build_L(def_states.G, n_states, params.n_x_max, params.n_z_max,
                                      num_G_nonzero_entries, params.mass,
                                      params.omega_x, params.omega_z,
                                      def_states.transition_lambda, def_states.B_direction, seed);

            cout << "2. Finished Building L!" << endl;

            cout << "Size of L: " << L.size() << endl;

            VectorXd t = VectorXd::LinSpaced(N, 0.0, t_0);
            cout << "dt: " << (t[1] - t[0]) << endl;

            VectorXcd psi0 = VectorXcd::Zero(size);
            int idx0 = from_tuple_to_number(2, params.n_x_init, params.n_z_init, params);
            psi0[idx0] = 1.0;

            auto [psi_final, jumps, nx_over_t, nz_over_t] =
                ss_spin::solve(t, psi0, H, L, def_states.G_tot, def_states.n_ground_states, n_states, params.n_x_max, params.n_z_max, num_keys);

            cout << "3. Function returned!" << endl;

            tot_jumps(D_index, I_index) = jumps;
            avg_tempsx(D_index, I_index) = nx_over_t[nx_over_t.size() - 1];
            avg_tempsz(D_index, I_index) = nz_over_t[nz_over_t.size() - 1];

            cout << "Total jumps: " << jumps << endl;
            cout << "Resulting heat x: " << nx_over_t[nx_over_t.size() - 1] << endl;
            cout << "Resulting heat z: " << nz_over_t[nz_over_t.size() - 1] << endl;

            auto t1 = chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed = t1 - t0;
            cout << "Elapsed time (s) per 1e3 timestamps: " << (elapsed.count() * 1e3 / (N - 1)) << endl;
        }
    }

    cout << tot_jumps << endl;
    cout << avg_tempsx << endl;
    cout << avg_tempsz << endl;
    return {tot_jumps, avg_tempsx, avg_tempsz};
}

} // namespace simulation
