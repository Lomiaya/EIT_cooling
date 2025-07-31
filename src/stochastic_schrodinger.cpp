#include "stochastic_schrodinger.hpp"
#include <unsupported/Eigen/MatrixFunctions> // for matrix exponential
#include <numeric>
#include <iostream>
#include <random>

namespace stochastic_schrodinger {

// imaginary unit
static const Complex I(0,1);

MatrixXc expmt(const MatrixXc& M, const MatrixXc& Mu, double dt) {
    return (-I * M * dt - Mu * dt).exp();
}

std::tuple<VectorXc, int> step(
    const VectorXc& psi_in,
    int scatter_count,
    const std::vector<MatrixXc>& Lt,
    const MatrixXc& expH,
    double dt,
    std::mt19937& rng)
{
    // Non-stochastic evolution step
    VectorXc psi = expH * psi_in;
    psi.normalize();

    // Calculate jump probabilities
    std::vector<double> probs;
    probs.reserve(Lt.size());
    double sum_probs = 0.0;

    for (const auto& Li : Lt) {
        VectorXc temp = Li * psi;
        double p = temp.squaredNorm() * dt;
        probs.push_back(p);
        sum_probs += p;
    }

    // Probability of no jump
    double p_no_jump = std::max(0.0, 1.0 - sum_probs);

    // Prepare probabilities vector including no jump
    probs.push_back(p_no_jump);

    // Create discrete distribution once with all probabilities (jumps + no jump)
    std::discrete_distribution<int> dist(probs.begin(), probs.end());
    int idx = dist(rng);

    if (idx < static_cast<int>(Lt.size())) {
        // jump occurs
        psi = Lt[idx] * psi;
        psi.normalize();
        ++scatter_count;
    }
    // else no jump, psi already normalized

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
        sum_LdaggerL.noalias() += Li.adjoint() * Li;
    }

    double dt = time[1] - time[0];
    MatrixXc expH = expmt(Ht, sum_LdaggerL, dt);

    VectorXc psi = psi0;
    int scatter_count = 0;

    for (size_t i = 0; i + 1 < time.size(); ++i) {
        std::tie(psi, scatter_count) = step(psi, scatter_count, Lt, expH, dt, rng);
    }

    return {psi, scatter_count, expH};
}

} // namespace stochastic_schrodinger
