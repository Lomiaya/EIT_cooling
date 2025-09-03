// chance_of_jump.cpp
#include "chance_of_jump.hpp"

#include <cmath>
#include <random>
#include <algorithm>
#include "constants.hpp"

//
#include <iostream>

namespace chance_of_jump {

// Factorial function (non-optimized, for small ints)
static double factorial(int n) {
    double res = 1.0;
    for (int i = 2; i <= n; ++i) res *= i;
    return res;
}

// Associated Laguerre polynomial L_n^a(x)
double assoc_laguerre(int n, int a, double x) {
    double acc = 0.0;
    for (int m = 0; m <= n; ++m) {
        double num = factorial(n + a);
        double denom = factorial(n - m) * factorial(a + m);
        double coeff = num / denom;
        acc += coeff * std::pow(-x, m) / factorial(m);
    }
    return acc;
}

// 1D harmonic oscillator matrix element
Complex ho_1d(int n_i, int n_f, double k_alpha, double x0_alpha) {
    int delta_n = n_f - n_i;
    double eta = k_alpha * x0_alpha / std::sqrt(2.0);
    int n_min = std::min(n_i, n_f);
    int n_max = std::max(n_i, n_f);

    double pref_abs = std::exp(-eta*eta/2) * std::sqrt(factorial(n_min)/factorial(n_max));
    Complex phase = std::pow(Complex(0,1)*eta, std::abs(delta_n));

    double L = assoc_laguerre(n_min, std::abs(delta_n), eta*eta);

    return phase * pref_abs * L;
}

// Newton iteration to solve inverse CDF for cos(theta)
double solve_cos_theta(double u, int np) {
    double c = 2.0*u - 1.0;
    for (int i = 0; i < 6; ++i) {
        double F, dF;
        if (np == 1) {
            F = 1 - 1.5*c + 0.5*c*c*c;
            dF = -1.5 + 1.5*c*c;
        } else {
            F = 1 - 0.75*c - 0.25*c*c*c;
            dF = -0.75 - 0.75*c*c;
        }
        c -= (F - u)/(dF + 1e-16);
    }
    if (c > 1.0) c = 1.0;
    else if (c < -1.0) c = -1.0;
    return c;
}

// Dot product helper
static double dot(const Vec3& a, const Vec3& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// Cross product helper
static Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3{
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    };
}

// Normalize vector
static Vec3 normalize(const Vec3& v) {
    double norm = std::sqrt(dot(v,v));
    if (norm < 1e-16) return Vec3{0,0,0};
    return Vec3{v[0]/norm, v[1]/norm, v[2]/norm};
}

// Calculate deterministic overlap <nf|e^{i k·r}|ni>
Complex calculate_overlap_from_k(
    const Vec3& k_vec,
    double mass,
    double omega_x,
    double omega_z,
    const std::array<int, 3>& start,
    const std::array<int, 3>& end)
{
    return calculate_overlap_from_k(k_vec, mass, omega_x, omega_x, omega_z, start, end);
}
Complex calculate_overlap_from_k(
    const Vec3& k_vec,
    double mass,
    double omega_x,
    double omega_y,
    double omega_z,
    const std::array<int, 3>& start,
    const std::array<int, 3>& end)
{
    // reduced Planck constant, from your constants.hpp or just hardcode here
    constexpr double hbar = 1.054571817e-34;

    double x0 = std::sqrt(hbar / (mass * omega_x));
    double y0 = std::sqrt(hbar / (mass * omega_y));
    double z0 = std::sqrt(hbar / (mass * omega_z));

    int n_xi = start[0], n_yi = start[1], n_zi = start[2];
    int n_xf = end[0], n_yf = end[1], n_zf = end[2];

    double kx = k_vec[0], ky = k_vec[1], kz = k_vec[2];

    Complex Mx = ho_1d(n_xi, n_xf, kx, x0);
    Complex My = ho_1d(n_yi, n_yf, ky, y0);
    Complex Mz = ho_1d(n_zi, n_zf, kz, z0);

    return Mx * My * Mz;
}

// Stochastic sampling of overlaps for random k directions
std::vector<double> calculate_chance_of_jump(
    int seed,
    int np,
    double mass,
    double omega_x,
    double omega_z,
    double wavelength,
    const Vec3& B_direction,
    const std::array<int,3>& start,
    const std::array<int,3>& end,
    int num_samples)
{
    return calculate_chance_of_jump(seed, np, mass, omega_x, omega_x, omega_z, wavelength, B_direction, start, end, num_samples);
}
std::vector<double> calculate_chance_of_jump(
    int seed,
    int np,
    double mass,
    double omega_x,
    double omega_y,
    double omega_z,
    double wavelength,
    const Vec3& B_direction,
    const std::array<int, 3>& start,
    const std::array<int, 3>& end,
    int num_samples)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> uniform_0_1(0.0, 1.0);
    std::uniform_real_distribution<double> uniform_0_2pi(0.0, 2*parameters::pi);

    Vec3 bz = normalize(B_direction);
    Vec3 helper{1.0, 0.0, 0.0};
    if (std::abs(dot(bz, helper)) > 1.0 - 1e-3) {
        helper = Vec3{0.0, 1.0, 0.0};
    }
    Vec3 bx = normalize(cross(bz, helper));
    Vec3 by = cross(bz, bx);

    double k_mag = 2 * parameters::pi / wavelength;

    std::vector<double> results;
    results.reserve(num_samples);

    for (int i=0; i<num_samples; ++i) {
        double u_theta = uniform_0_1(rng);
        double phi = uniform_0_2pi(rng);
        double c = solve_cos_theta(u_theta, np);
        double s = std::sqrt(std::max(0.0, 1.0 - c*c));

        Vec3 k_hat {
            s * std::cos(phi) * bx[0] + s * std::sin(phi) * by[0] + c * bz[0],
            s * std::cos(phi) * bx[1] + s * std::sin(phi) * by[1] + c * bz[1],
            s * std::cos(phi) * bx[2] + s * std::sin(phi) * by[2] + c * bz[2]
        };

        Vec3 k_vec {
            k_mag * k_hat[0],
            k_mag * k_hat[1],
            k_mag * k_hat[2]
        };

        Complex overlap = calculate_overlap_from_k(k_vec, mass, omega_x, omega_y, omega_z, start, end);

        results.push_back((overlap * std::conj(overlap)).real());
    }

    return results;
}

} // namespace chance_of_jump
