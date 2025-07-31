// stochastic_schrodinger.cpp
#include "stochastic_schrodinger.hpp"
#include <Eigen/Dense>
#include <random>
#include <numeric>
#include <iostream>
#include <unsupported/Eigen/MatrixFunctions> // for matrix exponential

namespace stochastic {

using namespace Eigen;
using namespace std::complex_literals;

ComplexMatrix expmt(const ComplexMatrix& H, const ComplexMatrix& Mu, double dt) {
    ComplexMatrix A = -1i * H * dt - Mu * dt;
    return A.exp();
}

StepResult step(
    const ComplexVector& psi_in,
    int scatter_count,
    const std::vector<ComplexMatrix>& L,
    const ComplexMatrix& expH,
    std::mt19937& gen
) {
    ComplexVector psi = expH * psi_in;
    psi.normalize();

    std::vector<double> probs;
    for (const auto& Li : L) {
        Complex val = psi.adjoint() * (Li.adjoint() * Li) * psi;
        probs.push_back(std::real(val));
    }

    double sum_probs = std::accumulate(probs.begin(), probs.end(), 0.0);
    probs.push_back(std::max(1.0 - sum_probs, 0.0)); // No-jump probability

    std::discrete_distribution<> dist(probs.begin(), probs.end());
    int idx = dist(gen);

    if (idx < static_cast<int>(L.size())) {
        ComplexVector psi_post = L[idx] * psi;
        psi_post.normalize();
        return {psi_post, scatter_count + 1};
    } else {
        return {psi, scatter_count};
    }
}

SolveResult solve(
    const std::vector<double>& time,
    const ComplexVector& psi0,
    const ComplexMatrix& H,
    const std::vector<ComplexMatrix>& L,
    unsigned seed
) {
    std::mt19937 gen(seed);
    double dt = time[1] - time[0];

    ComplexMatrix Mu = ComplexMatrix::Zero(H.rows(), H.cols());
    for (const auto& Li : L) {
        Mu += Li.adjoint() * Li;
    }

    ComplexMatrix expH = expmt(H, Mu, dt);
    ComplexVector psi = psi0;
    int jumps = 0;

    for (size_t i = 0; i < time.size() - 1; ++i) {
        StepResult res = step(psi, jumps, L, expH, gen);
        psi = res.psi;
        jumps = res.jumps;
    }

    return {psi, jumps, expH};
}

} // namespace stochastic
