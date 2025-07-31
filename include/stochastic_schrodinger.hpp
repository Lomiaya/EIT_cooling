// stochastic_schrodinger.hpp
#ifndef STOCHASTIC_SCHRODINGER_HPP
#define STOCHASTIC_SCHRODINGER_HPP

#include <complex>
#include <vector>
#include <random>
#include <array>

namespace stochastic_schrodinger {

using Complex = std::complex<double>;
using Matrix = std::vector<std::vector<Complex>>;
using Vector = std::vector<Complex>;

// Matrix exponential of -i*H*dt - Mu*dt
Matrix expmt(const Matrix& H, const Matrix& Mu, double dt);

// Time evolution of quantum state with possible jump
std::tuple<Vector, int, Matrix> solve(
    const std::vector<double>& time,
    const Vector& psi0,
    const Matrix& Ht,
    const std::vector<Matrix>& Lt,
    std::mt19937& rng);

} // namespace stochastic_schrodinger

#endif // STOCHASTIC_SCHRODINGER_HPP
