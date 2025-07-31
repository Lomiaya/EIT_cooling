#include "stochastic_schrodinger.hpp"
#include <unsupported/Eigen/MatrixFunctions> // for matrix exponential
#include <numeric>
#include <iostream>

namespace stochastic_schrodinger {

// i = imaginary unit
static const Complex I(0,1);

MatrixXc expmt(const MatrixXc& M, const MatrixXc& Mu, double dt) {
    MatrixXc A = -I * M * dt - Mu * dt;
    return A.exp();
}

std::tuple<VectorXc, int> step(
    const VectorXc& psi_in,
    int scatter_count,
    const std::vector<MatrixXc>& Lt,
    const MatrixXc& expH,
    const double dt,
    std::mt19937& rng)
{
    // non-stochastic evolution
    VectorXc psi = expH * psi_in;
    psi /= psi.norm();

    // calculate jump probabilities
    std::vector<double> probs;
    for (const auto& Li : Lt) {
        VectorXc temp = Li * psi;
        double p = (temp.adjoint() * temp)(0,0).real() * dt;
        probs.push_back(p);
    }

    // probability of no jump
    double sum_probs = std::accumulate(probs.begin(), probs.end(), 0.0);
    double p_no_jump = std::max(0.0, 1.0 - sum_probs);

    // create discrete distribution including no jump
    std::discrete_distribution<int> dist({probs.begin(), probs.end()});
    // add no jump index at end:
    std::vector<double> extended_probs = probs;
    extended_probs.push_back(p_no_jump);
    std::discrete_distribution<int> dist_full(extended_probs.begin(), extended_probs.end());

    int idx = dist_full(rng);

    if (idx < static_cast<int>(Lt.size())) {
        // jump
        psi = Lt[idx] * psi;
        scatter_count += 1;
    }
    // else no jump

    psi /= psi.norm();

    return {psi, scatter_count};
}

std::tuple<VectorXc, int, MatrixXc> solve(
    const std::vector<double>& time,
    const VectorXc& psi0,
    const MatrixXc& Ht,
    const std::vector<MatrixXc>& Lt,
    unsigned int seed)
{
    std::mt19937 rng(seed);

    // Compute sum_LdaggerL = sum Li^dagger * Li
    MatrixXc sum_LdaggerL = MatrixXc::Zero(Ht.rows(), Ht.cols());
    for (const auto& Li : Lt) {
        sum_LdaggerL += Li.adjoint() * Li;
    }

    double dt = time[1] - time[0];
    MatrixXc expH = expmt(Ht, sum_LdaggerL, dt);

    VectorXc psi = psi0;
    int scatter_count = 0;

    for (size_t i = 0; i < time.size() - 1; ++i) {
        std::tie(psi, scatter_count) = step(psi, scatter_count, Lt, expH, dt, rng);
    }

    return {psi, scatter_count, expH};
}

} // namespace stochastic_schrodinger
