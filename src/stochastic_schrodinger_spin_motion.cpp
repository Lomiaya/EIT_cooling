// stochastic_schrodinger_spin_motion.cpp
#include "stochastic_schrodinger_spin_motion.hpp"
#include "chance_of_jump.hpp"
#include <unsupported/Eigen/MatrixFunctions> // for expm
#include <cmath>
#include <numeric>
#include <omp.h>

#include <thread>
#include <vector>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <future>

//
#include <iostream>

namespace ss_spin {

/**
 * @brief Calculate evolution matrix
 */
SparseMat brute_force(const SpectrumMatrix& Hi_omegas,
                     double t0, double t1,
                     int n_steps) {
    double dt = (t1 - t0) / n_steps;
    int d = Hi_omegas[0].first.rows();
    SparseMat U = SparseMat(d,d);
    U.setIdentity();
    SparseMat Id = U; // ??? deep copy ???
    for(int k=0; k<n_steps; ++k) {
        double t = t0 + k*dt;
        SparseMat drive = SparseMat(d,d);
        for(auto& p: Hi_omegas) drive += p.first * std::exp(Complex(0,1)*p.second*t);
        U = (Id - Complex(0,1)*drive*dt) * U;
    }
    return U;
}

// MatrixXc magnus2_analytic(const SpectrumMatrix& Hi_omegas,
//                           double t0, double t1) {
//     auto K = [&](double w){
//         if(std::abs(w) < 1e-12) return Complex(0,1)*(t0-t1);
//         return (std::exp(Complex(0,1)*w*t0) - std::exp(Complex(0,1)*w*t1))/w;
//     };
//     MatrixXc Omega1 = -Complex(0,1)*(MatrixXc::Zero(Hi_omegas[0].first.rows(),Hi_omegas[0].first.cols()));
//     for(auto& p: Hi_omegas) Omega1 += p.first * K(p.second);
//     // Omega2
//     MatrixXc Omega2 = -Complex(0,1)*(MatrixXc::Zero(Hi_omegas[0].first.rows(),Hi_omegas[0].first.cols()));
//     // I integral
//     auto Iij = [&](double w_i,double w_j){
//         if(std::abs(w_i)<1e-12 && std::abs(w_j)<1e-12) return Complex(1,0)*0.5*(t1-t0)*(t1-t0);
//         if(std::abs(w_j)<1e-12) return std::exp(Complex(0,1)*w_i*t1)*(Complex(1,0)/(w_i*w_i)+(t1-t0)/(Complex(0,1)*w_i))
//                                         - std::exp(Complex(0,1)*w_i*t0)*(Complex(1,0)/(w_i*w_i));
//         if(std::abs(w_i)<1e-12) return (std::exp(Complex(0,1)*(w_i+w_j)*t1)-std::exp(Complex(0,1)*(w_i+w_j)*t0))/(-w_j*(w_i+w_j)) 
//                                         - std::exp(Complex(0,1)*w_j*t0)*(t1-t0)/(Complex(0,1)*w_j);
//         if(std::abs(w_i+w_j)<1e-12) return ((t1-t0))/(Complex(0,1)*w_j)
//                                             - std::exp(Complex(0,1)*w_j*t0)*(std::exp(Complex(0,1)*w_i*t1)-std::exp(Complex(0,1)*w_i*t0))/(-w_i*w_j);
//         return (std::exp(Complex(0,1)*(w_i+w_j)*t1)-std::exp(Complex(0,1)*(w_i+w_j)*t0))/(-w_j*(w_i+w_j))
//                 - std::exp(Complex(0,1)*w_j*t0)*(std::exp(Complex(0,1)*w_i*t1)-std::exp(Complex(0,1)*w_i*t0))/(-w_i*w_j);
//     };
//     // drive-drive
//     for(auto& pi: Hi_omegas) for(auto& pj: Hi_omegas) {
//         double wi=pi.second, wj=pj.second;
//         Omega2 += -0.5*(pi.first*pj.first - pj.first*pi.first)*Iij(pi.second,pj.second);
//     }
//     MatrixXc Omega = Omega1 + Omega2;
//     return Omega.exp();
// }

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

std::tuple<double,double,double> compute_avg_x_y_z(const VectorXc& psi, int n_x, int n_y, int n_z) {
    int n_tot = psi.size();
    VecD probs = psi.cwiseAbs2().real();
    VecD idx = VecD::LinSpaced(n_tot, 0, n_tot - 1);

    // compute 3D indices: s, x, z from linear index
    VecD s = (idx / (n_x * n_y * n_z)).array().floor();
    VecD x = ((idx - s * n_x * n_y * n_z) / (n_y * n_z)).array().floor();
    VecD y = ((idx - s * n_x * n_y * n_z - x * n_y * n_z) / n_z).array().floor();
    VecD z = ((idx - s * n_x * n_y * n_z - x * n_y * n_z - y * n_z)).array();

    // ignore spin s for average
    double avg_x = (probs.array() * x.array()).sum();
    double avg_y = (probs.array() * y.array()).sum();
    double avg_z = (probs.array() * z.array()).sum();
    return {avg_x, avg_y, avg_z};
}

Eigen::VectorXi expand(int idx,
                       const std::array<int,5>& size_list) {
    int np = size_list[0];
    int ne = size_list[1];
    int ng = size_list[2];
    int n_x = size_list[3];
    int n_z = size_list[4];
    int nzi = idx % n_z;
    idx = idx / n_z;
    int nxi = idx % n_x;
    idx = idx / n_x;
    int nzf = idx % n_z;
    idx = idx / n_z;
    int nyf = idx % n_z;
    idx = idx / n_z;
    int nxf = idx % n_x;
    idx = idx / n_x;
    int ng_ = idx % ng;
    idx = idx / ng;
    int ne_ = idx % ne;
    idx = idx / ne;
    int np_ = idx;
    Eigen::VectorXi out(7);
    out << np_, ne_, ng_, nxf, nyf, nzf, nxi, nzi;
    return out;
}

std::vector<std::tuple<int, int, double>> build_L(
    const std::vector<Eigen::MatrixXd>& G,
    int n_s, int n_x, int n_z,
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

        int excited_state = from_tuple_to_number(ne_ + ng, nxi, nzi);
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
        double chance_of_jump = std::sqrt(average);

        Lt.emplace_back(ground_state, excited_state, chance_of_jump);
    }

    return Lt;
}

// Kernel function that consumes Data
std::pair<StateCarry, void*> step(const StateCarry& carry,
                                  SparseMat& expH,
                                  SparseMat& Wt) {
    auto [psi, scatter_count, ii, step_count, per_log_step, nx, ny, nz, Lt, G_tot, n_x, n_y, n_z, is_2d_sim, dt, rng] = carry;

    // Non-stochastic evolution using Magnus expansion
    psi = expH * psi;
    psi.normalize();

    VectorXc excited = Wt * psi;
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
        psi_post(i) = 1.0;
        psi = psi_post;
        scatter_count += 1;
    }

    // Record observable
    if (step_count % per_log_step == 0) {
        if (is_2d_sim) {
            auto [nxi, nzi] = compute_avg_x_z(psi, n_x, n_z);
            nx(ii) = nxi;
            ny(ii) = nxi;
            nz(ii) = nzi;
        } else {
            auto [nxi, nyi, nzi] = compute_avg_x_y_z(psi, n_x, n_y, n_z);
            nx(ii) = nxi;
            ny(ii) = nyi;
            nz(ii) = nzi;
        }
        ii += 1;
    }
    step_count += 1;

    StateCarry new_state = std::make_tuple(psi, scatter_count, ii, step_count, per_log_step, nx, ny, nz,
                                           Lt, G_tot, n_x, n_y, n_z, is_2d_sim, dt, std::move(rng));
    return {new_state, nullptr};
}

/**
 * @brief A parallelized simulation of a quantum system's evolution over time, using producer-consumer threading to manage data processing. 
 * It returns a tuple containing the final quantum state, the average number of quantum jumps, and averaged statistical data about the system's evolution.
 * @param time_step Time step dt
 * @param num_steps Number of the steps, N=t/dt
 * @param per_long_step Log data every per_long_step steps
 * @param psi0 Initial quantum state
 * @param H Hamiltonian
 * @param W Auxilary matirx
 * @param Lt Relaxation operator
 * @param G_tot Total scattering rate
 * @param n_x, n_y, n_z Number of motional state
 * @param is_2d_sim Is the simulation 2-D
 * @param num_keys Number of consumer threads & loops per thread
 * @param low_pass_threshold Threshold to filter high-frequency components in Hamiltonian and W terms.
 */
std::tuple<VectorXc,double,VecD,VecD,VecD>
solve(const double time_step,
      const long long num_steps,
      const long long per_log_step,
      const VectorXc& psi0,
      const SpectrumMatrix& H,
      const SpectrumMatrix& W,
      const std::vector<std::tuple<int, int, double>>& Lt,
      const double G_tot,
      int n_x, int n_y, int n_z,
      bool is_2d_sim,
      int num_keys,
      double low_pass_threshold)
{
    const size_t N = num_keys;   // Number of consumer threads & loops per thread
    const size_t M = 32;   // Number of producer threads
    const size_t stash_capacity = 100;   // Stash capacity

    std::queue<std::tuple<SparseMat, SparseMat, size_t, size_t>> stash;
    std::mutex mtx;
    std::condition_variable cv_produce;
    std::condition_variable cv_consume;

    std::atomic<size_t> produce_count{0};
    std::atomic<size_t> consume_count{0};
    std::atomic<size_t> next_expected_idx{0}; // controls sequential emplacement
    
    const size_t total = num_steps; // Total number of elements to produce/consume

    SpectrumMatrix H_eff = H;

    // Producer threads
    auto producer = [&]() {
        while (true) {
            size_t idx = produce_count.fetch_add(1);
            if (idx >= total) break;

            // Wait until it's this producer's turn to emplace
            {
                // Produce data
                double time0 = time_step * idx;
                double time1 = time0 + time_step;
                // MatrixXc d = magnus2_analytic(H_eff, time0, time1);
                SparseMat d = brute_force(H_eff, time0, time1, 1);
                SparseMat w = evaluate(W, time0);
                std::unique_lock<std::mutex> lock(mtx);
                cv_produce.wait(lock, [&]() { return idx == next_expected_idx && stash.size() < stash_capacity; });
                stash.emplace(d, w, 0, idx);  // stash maintains insertion order
                ++next_expected_idx;
                lock.unlock();

                if (idx % per_log_step == 0) {
                    std::cout << "Produced element " << idx << " by producer "
                            << std::this_thread::get_id() << std::endl;
                }
            }

            // Notify all consumers and producers (one may be waiting to insert next)
            cv_consume.notify_all();
            cv_produce.notify_all();
        }
    };

    // Consumer threads
    auto make_consumer = [&](StateCarry init_state, std::promise<StateCarry> result_promise) {
        return [&, state = std::move(init_state), result = std::move(result_promise)]() mutable {
            for (size_t i = 0; i < total; ++i) {
                {
                    std::unique_lock<std::mutex> lock(mtx);
                    cv_consume.wait(lock, [&]() { return !stash.empty() && std::get<3>(stash.front()) == i; });

                    auto& [d, w, count, idx] = stash.front();

                    // std::cout << "Consumed element " << idx << ", " << i << " by consumer "
                    //         << std::this_thread::get_id() << std::endl;

                    lock.unlock();

                    // Process the data
                    state = step(state, d, w).first; // actively process data since d might be freed after stash.pop().
                }

                {
                    std::unique_lock<std::mutex> lock(mtx);
                    cv_consume.wait(lock, [&]() { return !stash.empty() && std::get<3>(stash.front()) == i; });

                    auto& [d, w, count, idx] = stash.front();
                    ++count;

                    if (count == N) {
                        stash.pop();
                        lock.unlock();
                        cv_produce.notify_all(); // Need to be all on some hardware.
                        cv_consume.notify_all();
                    } else {
                        lock.unlock(); // Don't notify producer yet
                    }
                }
            }

            // Return the final state
            result.set_value(state);
        };
    };


    // Launch producers
    std::vector<std::thread> producers;
    for (size_t i = 0; i < M; ++i) {
        producers.emplace_back(producer);
    }

    // Launch consumers
    std::vector<std::thread> consumers;
    std::vector<std::future<StateCarry>> consumer_results;
    for (size_t i = 0; i < N; ++i) {
        const auto& seed = i;
        std::mt19937 rng(seed);
        size_t num_log = (num_steps - 1) / per_log_step + 1;
        VecD nx = VecD::Zero(num_log);
        VecD ny = VecD::Zero(num_log);
        VecD nz = VecD::Zero(num_log);
        StateCarry init_state = std::make_tuple(
            psi0,
            0, // scatter_count
            0, // ii
            0, // step_count
            per_log_step,
            nx,
            ny,
            nz,
            Lt,
            G_tot,
            n_x, n_y, n_z,
            is_2d_sim,
            time_step,
            rng
        );
        std::promise<StateCarry> p;
        consumer_results.emplace_back(p.get_future());
        consumers.emplace_back(make_consumer(init_state, std::move(p)));
    }


    // Join all
    for (auto &t : producers) t.join();
    for (auto &t : consumers) t.join();

    // Collect results
    std::vector<StateCarry> final_states;
    for (auto& fut : consumer_results) {
        final_states.push_back(fut.get());
    }
    
    // Rearrange collected results
    std::vector<VecD> all_nx;
    std::vector<VecD> all_ny;
    std::vector<VecD> all_nz;
    size_t num_log = (num_steps - 1) / per_log_step + 1;
    VecD avg_nx = VecD::Zero(num_log);
    VecD avg_ny = VecD::Zero(num_log);
    VecD avg_nz = VecD::Zero(num_log);
    std::vector<int> all_jumps;
    std::vector<VectorXc> all_final_psis;
    for (auto& final_state: final_states) {
        auto& [psi, scatter_count, ii, step_count, per_log_step, nx, ny, nz, Lt, G_tot, n_x, n_y, n_z, is2dsim, dt, rng] = final_state;
        all_final_psis.push_back(psi);
        all_jumps.push_back(scatter_count);
        all_nx.push_back(nx);
        all_ny.push_back(ny);
        all_nz.push_back(nz);
    }

    // Average nx and nz
    for (int i = 0; i < num_log; ++i) {
        for (int k = 0; k < all_nx.size(); ++k) {
            avg_nx(i) += all_nx[k](i);
            avg_ny(i) += all_ny[k](i);
            avg_nz(i) += all_nz[k](i);
        }
        avg_nx(i) /= static_cast<double>(all_nx.size());
        avg_ny(i) /= static_cast<double>(all_ny.size());
        avg_nz(i) /= static_cast<double>(all_nz.size());
    }

    // Average number of jumps
    double avg_jumps = std::accumulate(all_jumps.begin(), all_jumps.end(), 0.0) /
                       static_cast<double>(all_jumps.size());

    // Return final psi of last run (just like in Python)
    return {all_final_psis.back(), avg_jumps, avg_nx, avg_ny, avg_nz};
}

} // namespace ss_spin