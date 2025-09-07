#include "simulation.hpp"
#include "stochastic_schrodinger_spin_motion.hpp"
#include "hamiltonian.hpp"
#include "define_params.hpp"
#include "constants.hpp"

#include <iostream>
#include <fstream>

namespace simulation {

using namespace ss_spin;
using namespace hamiltonian;
using namespace define_params;
using namespace std;

void write_to_file(std::string file_name, std::string label, const VecD& data) {
    // Open file for writing
    std::ofstream file(file_name, std::ios::out | std::ios::app);
    if (!file.is_open()) {
        std::cerr << "Error: could not open file for writing. " << file_name;
    }

    // Optional: format without brackets, with space separation
    Eigen::IOFormat format(Eigen::FullPrecision, Eigen::DontAlignCols, " ", "\n");

    // Write matrix to file
    file << label << " " << data.format(format) << std::endl;
}

/**
 * @brief Runs quantum simulation over grids of detunings and intensities, and computes heating and jump statistics.
 *
 * This function simulates the dynamics of a quantum system under a time-dependent Hamiltonian and associated 
 * Lindblad operators, over a range of detunings (D) and intensities (I). It diagonalizes the initial state 
 * Hamiltonian, constructs Liouvillian, Hamiltonian, and jump operators, and evolves the system using a solver.
 * It outputs total jump counts and average final occupation numbers along two spatial axes.
 *
 * @param params         Struct containing simulation parameters
 * @param def_states     Struct describing the default quantum states
 * @param N              Number of time steps in the simulation.
 * @param t_0            Total simulation time (final time).
 * @param num_keys       Number of keys/frequencies used for Hamiltonian/W modulation.
 * @param low_pass_threshold  Threshold to filter high-frequency components in Hamiltonian and W terms.
 *
 * @return std::tuple<MatrixXd, MatrixXd, MatrixXd> 
 *         A tuple containing:
 *         - tot_jumps:     Matrix of total quantum jumps per (D, I) pair.
 *         - avg_tempsx:    Matrix of average final x-direction phonon occupation per (D, I) pair.
 *         - avg_tempsz:    Matrix of average final z-direction phonon occupation per (D, I) pair.
 *
 * ### Details:
 * - **Isat Calculation**: Uses physical constants and decay rates to compute saturation intensity.
 * - **Diagonalization**: Diagonalizes the Hamiltonian to extract energy eigenstates and decay matrix G.
 * - **Operator Construction**: Builds time-dependent Hamiltonian (H), jump operators (W), and Liouvillian (L).
 * - **Initial State**: Initializes wavefunction in ground state with given spatial mode indices.
 * - **Solver**: Evolves system using a solver (`ss_spin::solve`) and collects statistics on quantum jumps and phonon heating.
 * - **Output Logging**: Intermediate results are logged to console and to files `heat_x.txt` and `heat_z.txt`.
 * - **Performance**: Reports per-timestep runtime performance for each parameter set.
 *
 * ### Notes:
 * - The matrices returned have dimensions `(Ds.rows(), Is.rows())`.
 * - Assumes fixed initial random seed (42) for reproducibility.
 * - Hamiltonian and W undergo "cleanup" filtering to reduce high-frequency contributions.
 *
 * ### Example Usage:
 * ```cpp
 * auto [jumps, heat_x, heat_z] = simulate(params, def_states, 500, 1e-3, 100, 0.01);
 * ```
 */
std::tuple<MatrixXd, MatrixXd, MatrixXd, MatrixXd> simulate(const Params& params, const States& def_states, const int N, const double t_0, const int num_keys, double low_pass_threshold)
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
    MatrixXd avg_tempsy(Ds.rows(), Is.rows());
    MatrixXd avg_tempsz(Ds.rows(), Is.rows());
    double dt = t_0 / (N - 1);
    cout << "dt: " << dt << endl;
    for (size_t D_index = 0; D_index < Ds.rows(); ++D_index) {
        for (size_t I_index = 0; I_index < Is.rows(); ++I_index) {

            cout << D_index << ", " << I_index << endl;

            auto t0 = chrono::high_resolution_clock::now();
            
            unsigned int seed = 42;
            int num_G_nonzero_entries = 0;
            for (const auto& g : states.G) num_G_nonzero_entries += (g.array() != 0.0).count();

            std::vector<std::tuple<int, int, double>> L;
            VectorXcd psi0;
            if (params.do_2d_sim) {
                cout << "2D simulation!" << endl;
                L = build_L_2d(states.G, params.n_x_max, params.n_z_max,
                                num_G_nonzero_entries, params.mass,
                                params.omega_x, params.omega_z,
                                states.transition_lambda, states.B_direction, seed);
                int n_expanded_ground_states = states.n_ground_states * params.n_x_max * params.n_z_max;
                psi0 = VectorXc::Zero(n_expanded_ground_states);
                int idx0 = params.n_x_init * params.n_z_max + params.n_z_init; // initial state index
                psi0[idx0] = 1.0;
            } else {
                cout << "3D simulation!" << endl;
                L = build_L_3d(states.G, params.n_x_max, params.n_y_max, params.n_z_max,
                                num_G_nonzero_entries, params.mass,
                                params.omega_x, params.omega_y, params.omega_z,
                                states.transition_lambda, states.B_direction, seed);
                int n_expanded_ground_states = states.n_ground_states * params.n_x_max * params.n_y_max * params.n_z_max;
                psi0 = VectorXc::Zero(n_expanded_ground_states);
                int idx0 = ((params.n_x_init * params.n_y_max) + params.n_y_init) * params.n_z_max + params.n_z_init; // initial state index
                psi0[idx0] = 1.0;
            }
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
            //     cout << "H: " << H_freq << ", " << H_mat.norm()
            //     // << ", " << H_mat << endl
            //     ;
            // }

            // for (const auto& [W_mat, W_freq] : W) {
            //     cout << "W: " << W_freq << ", " << W_mat.norm()
            //     // << ", " << W_mat << endl
            //     ;
            // }

            auto [psi_final, jumps, nx_over_t, ny_over_t, nz_over_t] =
                ss_spin::solve(dt, N, N / 10, psi0, H, W, L, states.G_tot, 
                               params.n_x_max, params.n_y_max, params.n_z_max, 
                               params.do_2d_sim, num_keys, low_pass_threshold);

            cout << "3. Function returned!" << endl;

            tot_jumps(D_index, I_index) = jumps;
            avg_tempsx(D_index, I_index) = nx_over_t[nx_over_t.size() - 1];
            avg_tempsy(D_index, I_index) = ny_over_t[nx_over_t.size() - 1];
            avg_tempsz(D_index, I_index) = nz_over_t[nz_over_t.size() - 1];

            cout << "Total jumps: " << jumps << endl;
            cout << "Resulting heat x: " << nx_over_t[nx_over_t.size() - 1] << endl;
            cout << "Resulting heat y: " << ny_over_t[ny_over_t.size() - 1] << endl;
            cout << "Resulting heat z: " << nz_over_t[nz_over_t.size() - 1] << endl;
            cout << "Resulting psi:" << psi_final << endl;

            write_to_file("./heat_x.txt", "nx over t: ", nx_over_t);
            write_to_file("./heat_y.txt", "ny over t: ", ny_over_t);
            write_to_file("./heat_z.txt", "nz over t: ", nz_over_t);

            auto t1 = chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed = t1 - t0;
            cout << "Elapsed time (s) per 1e3 timestamps: " << (elapsed.count() * 1e3 / (N - 1)) << endl;
        }
    }

    cout << tot_jumps << endl;
    cout << avg_tempsx << endl;
    cout << avg_tempsy << endl;
    cout << avg_tempsz << endl;
    return {tot_jumps, avg_tempsx, avg_tempsy, avg_tempsz};
}

}