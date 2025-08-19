#include "simulation.hpp"
#include "stochastic_schrodinger_spin_motion.hpp"
#include "hamiltonian.hpp"
#include "define_params.hpp"
#include "constants.hpp"

#include <iostream>

namespace simulation {

using namespace ss_spin;
using namespace hamiltonian;
using namespace define_params;
using namespace std;

std::tuple<MatrixXd, MatrixXd, MatrixXd> simulate(const Params& params, const States& def_states, const int N, const double t_0, const int num_keys, double low_pass_threshold)
{
    auto Isat = parameters::pi * parameters::plancks_constant * parameters::speed_of_light * def_states.G_tot /
                  (3.0 * pow(def_states.transition_lambda, 3));
    cout << "Isat: " << Isat << endl;

    const auto& Ds = params.D;
    const auto& Is = params.I;
    auto states = diagonalize_hamiltonian(def_states);
    // std::cout << states.G[0] << std::endl;
    // std::cout << states.G[1] << std::endl;
    // std::cout << states.G[2] << std::endl;
    // std::cout << states.H_ground << std::endl;
    // std::cout << states.H_excited << std::endl;
    MatrixXd tot_jumps(Ds.rows(), Is.rows());
    MatrixXd avg_tempsx(Ds.rows(), Is.rows());
    MatrixXd avg_tempsz(Ds.rows(), Is.rows());
    int n_expanded_ground_states = states.n_ground_states * params.n_x_max * params.n_z_max;
    double dt = t_0 / (N - 1);
    cout << "dt: " << dt << endl;
    for (size_t D_index = 0; D_index < Ds.rows(); ++D_index) {
        for (size_t I_index = 0; I_index < Is.rows(); ++I_index) {

            cout << D_index << ", " << I_index << endl;

            auto t0 = chrono::high_resolution_clock::now();
            
            unsigned int seed = 42;
            int num_G_nonzero_entries = 0;
            for (const auto& g : states.G) num_G_nonzero_entries += (g.array() != 0.0).count();

            auto L = build_L(states.G, params.n_x_max, params.n_z_max,
                             num_G_nonzero_entries, params.mass,
                             params.omega_x, params.omega_z,
                             states.transition_lambda, states.B_direction, seed);
            std::cout << "Finished building L!" << std::endl;
            auto W = build_W(states, params, I_index, D_index);
            W = cleanup(W);
            std::cout << "Finished building W!" << std::endl;
            auto H = build_H(states, params, I_index, D_index, W, low_pass_threshold);
            H = cleanup(H);
            std::cout << "Finished building H!" << std::endl;
            std::cout << "Size of W: " << W.size() << std::endl;
            std::cout << "Size of H: " << H.size() << std::endl;

            // for (const auto& [g, e, v] : L) {
            //     cout << "L: " << g << ", " << e << ", " << v << endl;
            // }

            // for (const auto& [H_mat, H_freq] : H) {
            //     cout << "H: " << H_freq << ", " << H_mat << endl;
            // }

            // for (const auto& [W_mat, W_freq] : W) {
            //     cout << "W: " << W_freq << ", " << W_mat << endl;
            // }

            VectorXcd psi0 = VectorXc::Zero(n_expanded_ground_states);
            int idx0 = params.n_x_init * params.n_z_max + params.n_z_init; // initial state index
            psi0[idx0] = 1.0;

            auto [psi_final, jumps, nx_over_t, nz_over_t] =
                ss_spin::solve(dt, N, N / 10, psi0, H, W, L, states.G_tot, params.n_x_max, params.n_z_max, num_keys, low_pass_threshold);

            cout << "3. Function returned!" << endl;

            tot_jumps(D_index, I_index) = jumps;
            avg_tempsx(D_index, I_index) = nx_over_t[nx_over_t.size() - 1];
            avg_tempsz(D_index, I_index) = nz_over_t[nz_over_t.size() - 1];

            cout << "Total jumps: " << jumps << endl;
            cout << "Resulting heat x: " << nx_over_t[nx_over_t.size() - 1] << endl;
            cout << "Resulting heat z: " << nz_over_t[nz_over_t.size() - 1] << endl;
            cout << "Resulting psi:" << psi_final << endl;

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

}