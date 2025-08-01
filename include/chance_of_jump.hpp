// chance_of_jump.hpp
#ifndef CHANCE_OF_JUMP_HPP
#define CHANCE_OF_JUMP_HPP

#include <array>
#include <complex>
#include <vector>

namespace chance_of_jump {

using Complex = std::complex<double>;
using Vec3 = std::array<double, 3>;

// Computes associated Laguerre polynomial L_n^a(x)
double assoc_laguerre(int n, int a, double x);

// Harmonic oscillator matrix element in 1D
Complex ho_1d(int n_i, int n_f, double k_alpha, double x0_alpha);

// Solve for cos(theta) given uniform sample u and transition type np
double solve_cos_theta(double u, int np);

// Deterministic matrix element for <nf|e^{i k.r}|ni>
Complex calculate_overlap_from_k(
    const Vec3& k_vec,
    double mass,
    double omega_x,
    double omega_z,
    const std::array<int, 3>& start,
    const std::array<int, 3>& end);

// Sample matrix elements for random k-directions
std::vector<double> calculate_chance_of_jump(
    int seed,
    int np,
    double mass,
    double omega_x,
    double omega_z,
    double wavelength,
    const Vec3& B_direction,
    const std::array<int, 3>& start,
    const std::array<int, 3>& end,
    int num_samples = 500);

} // namespace chance_of_jump

#endif // CHANCE_OF_JUMP_HPP
