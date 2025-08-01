// stochastic_schrodinger.hpp
#ifndef STOCHASTIC_SCHRODINGER_HPP
#define STOCHASTIC_SCHRODINGER_HPP

#include <Eigen/Dense>
#include <complex>
#include <vector>
#include <random>
#include <tuple>

namespace stochastic_schrodinger {

using Complex = std::complex<double>;
using MatrixXc = Eigen::MatrixXcd;
using VectorXc = Eigen::VectorXcd;

// Matrix exponential for small dt using Eigen's matrix exponential (requires unsupported/Eigen/MatrixFunctions)
MatrixXc expmt(const MatrixXc& M, const MatrixXc& Mu, double dt);

// Single time step of the stochastic evolution
std::tuple<VectorXc, int> step(
    const VectorXc& psi,
    int scatter_count,
    const std::vector<MatrixXc>& Lt,
    const MatrixXc& expH,
    std::mt19937& rng);

// Solve the stochastic Schrodinger equation over given time points
// Returns final psi, total jumps count, and expH matrix used
std::tuple<VectorXc, int, MatrixXc> solve(
    const std::vector<double>& time,
    const VectorXc& psi0,
    const MatrixXc& Ht,
    const std::vector<MatrixXc>& Lt,
    unsigned int seed);

} // namespace stochastic_schrodinger

#endif // STOCHASTIC_SCHRODINGER_HPP
