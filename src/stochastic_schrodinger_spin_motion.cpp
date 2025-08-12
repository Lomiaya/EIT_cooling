#include "stochastic_schrodinger_spin_motion.hpp"
#include <thread>
#include <future>
#include <unsupported/Eigen/MatrixFunctions>

namespace ss_spin {

MatrixXc propagator(const MatrixXc& H, double dt) {
    MatrixXc exp_H = (H * dt).exp();
    return exp_H;
}

std::pair<double,double> compute_avg_x_z(const VectorXc& psi,
                                         int n_x, int n_z) {
    int n_tot = psi.size();
    VecD probs = psi.cwiseAbs2().real();
    VecD idx = VecD::LinSpaced(n_tot, 0, n_tot - 1);

    // compute 3D indices: s, x, z from linear index
    VecD s = (idx / (n_x * n_z)).array().floor();
    VecD x = ((idx - s * n_x * n_z) / n_z).array().floor();
    VecD z = (idx - s * n_x * n_z - x * n_z).array();

    // ignore spin s for average
    double avg_x = (probs.array() * x.array()).sum();
    double avg_z = (probs.array() * z.array()).sum();
    return {avg_x, avg_z};
}

std::pair<StateCarry, void*> step(const StateCarry& carry)
{
    auto [psi, scatter_count, ii, step_count, per_log_step, nx, nz,
          Lt, G_tot, U, W, n_x, n_z, dt, rng] = carry;

    // Propagate the state
    psi = U * psi;
    psi.normalize();
    
    // The excited state population
    VectorXc excited = W * psi;

    double jump_prob = excited.squaredNorm() * G_tot * dt;
    std::vector<double> probs_jump = {jump_prob, 1 - jump_prob};
    std::discrete_distribution<> dist_jump(probs_jump.begin(), probs_jump.end());
    int idx_jump = dist_jump(rng);

    if (idx_jump == 0) {
        // Compute probabilities
        std::vector<double> probs;
        probs.reserve(Lt.size());

        for (const auto& [i, j, val] : Lt) {
            Complex Li_psi = val * excited(j);
            Complex Li_dag_Li_psi_j = std::conj(val) * Li_psi;
            Complex dot = std::conj(excited(j)) * Li_dag_Li_psi_j;
            probs.push_back(dt * dot.real());
        }

        std::discrete_distribution<> dist(probs.begin(), probs.end());
        int idx = dist(rng);

        auto [i, j, val] = Lt[idx];
        VectorXc psi_post = VectorXc::Zero(psi.size());
        psi_post(i) = 1;
        psi = psi_post;
        scatter_count += 1;
    }

    // Record observable
    if (step_count % per_log_step == 0) {
        auto [nxi, nzi] = compute_avg_x_z(psi, n_x, n_z);
        nx(ii) = nxi;
        nz(ii) = nzi;
        ii += 1;
    }
    step_count += 1;
    
    StateCarry new_state = std::make_tuple(psi, scatter_count, ii, step_count, per_log_step, nx, nz,
                                           Lt, G_tot, U, W, n_x, n_z, dt, std::move(rng));
    // Return the updated state carry
    return {new_state, nullptr};
}

std::tuple<VectorXc,double,VecD,VecD>
solve(const double time_step,
      const long long num_steps,
      const long long per_log_step,
      const VectorXc& psi0,
      const MatrixXc& H,
      const MatrixXc& W,
      const std::vector<std::tuple<int, int, double>>& Lt,
      const double G_tot,
      int n_x, int n_z,
      int num_keys)
{
    const size_t N = num_keys;   // Number of consumer threads & loops per thread
    const size_t M = 10;   // Number of producer threads
    const size_t stash_capacity = 100;   // Stash capacity

    std::atomic<size_t> produce_count{0};
    std::atomic<size_t> consume_count{0};
    std::atomic<size_t> next_expected_idx{0}; // controls sequential emplacement

    const size_t total = num_steps; // Total number of elements to produce/consume

    // Consumer threads
    auto make_consumer = [&](StateCarry init_state, std::promise<StateCarry> result_promise) {
        return [&, state = std::move(init_state), result = std::move(result_promise)]() mutable {
            for (size_t i = 0; i < total; ++i) {
                state = step(state).first;
            }

            // Return the final state
            result.set_value(state);
        };
    };

    // Launch consumers
    std::vector<std::thread> consumers;
    std::vector<std::future<StateCarry>> consumer_results;
    for (size_t i = 0; i < N; ++i) {
        const auto& seed = i;
        std::mt19937 rng(seed);
        size_t num_log = (num_steps - 1) / per_log_step + 1;
        VecD nx = VecD::Zero(num_log);
        VecD nz = VecD::Zero(num_log);
        StateCarry init_state = std::make_tuple(
            psi0,
            0, // scatter_count
            0, // ii
            0, // step_count
            per_log_step,
            nx,
            nz,
            Lt,
            G_tot,
            propagator(H, time_step), // U
            W,
            n_x, n_z,
            time_step,
            rng
        );
        std::promise<StateCarry> p;
        consumer_results.emplace_back(p.get_future());
        consumers.emplace_back(make_consumer(init_state, std::move(p)));
    }

    // Join all
    for (auto &t : consumers) t.join();

    // Collect results
    std::vector<StateCarry> final_states;
    for (auto& fut : consumer_results) {
        final_states.push_back(fut.get());
    }
    
    // Rearrange collected results
    std::vector<VecD> all_nx;
    std::vector<VecD> all_nz;
    size_t num_log = (num_steps - 1) / per_log_step + 1;
    VecD avg_nx = VecD::Zero(num_log);
    VecD avg_nz = VecD::Zero(num_log);
    std::vector<int> all_jumps;
    std::vector<VectorXc> all_final_psis;
    for (auto& final_state: final_states) {
        auto& [psi, scatter_count, ii, step_count, per_log_step, nx, nz, Lt, G_tot, n_g, n_s, n_x, n_z, dt, rng] = final_state;
        all_final_psis.push_back(psi);
        all_jumps.push_back(scatter_count);
        all_nx.push_back(nx);
        all_nz.push_back(nz);
    }

    // Average nx and nz
    for (int i = 0; i < num_log; ++i) {
        for (int k = 0; k < all_nx.size(); ++k) {
            avg_nx(i) += all_nx[k](i);
            avg_nz(i) += all_nz[k](i);
        }
        avg_nx(i) /= static_cast<double>(all_nx.size());
        avg_nz(i) /= static_cast<double>(all_nz.size());
    }

    // Average number of jumps
    double avg_jumps = std::accumulate(all_jumps.begin(), all_jumps.end(), 0.0) /
                       static_cast<double>(all_jumps.size());

    // Return final psi of last run (just like in Python)
    return {all_final_psis.back(), avg_jumps, avg_nx, avg_nz};
}

}