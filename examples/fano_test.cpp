// fano_test.cpp
#include <iostream>
#include <vector>
#include <complex>
#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include "stochastic_schrodinger.hpp"

// Placeholder for your solve function, should be implemented elsewhere
// std::tuple<Eigen::VectorXcd, int, Eigen::MatrixXcd> solve(...);

using namespace std::complex_literals;
using namespace Eigen;

using cdouble = std::complex<double>;

MatrixXcd add_light_coupling_component(
    const MatrixXcd& H, int g, int e, cdouble D, cdouble O) 
{
    MatrixXcd H_new = H;
    H_new(g, e) += std::conj(O);
    H_new(e, g) += O;
    H_new(g, g) += D;
    return H_new;
}

MatrixXcd decay_rate(int size, int g, int e, cdouble G) {
    MatrixXcd L = MatrixXcd::Zero(size, size);
    L(g, e) = std::sqrt(G);
    return L;
}

int main() {
    // Constants
    constexpr double pi = 3.14159265358979323846;
    constexpr double MHz = 1e6;

    // Define parameters
    const int size = 3;

    std::vector<cdouble> D1s;
    // Fill D1s with linspace -100 to 0 (20 points) and 0 to 1 (20 points), * 2*pi*1e6 + 2*pi*94.5e6
    // Approximate linspace:
    for (int i = 0; i < 20; ++i) {
        double val = -100.0 + i * (100.0 / 19.0);  // -100 to 0
        D1s.push_back(2 * pi * 94.5e6 + 2 * pi * val * 1e6);
    }
    for (int i = 0; i < 20; ++i) {
        double val = i * (1.0 / 19.0); // 0 to 1
        D1s.push_back(2 * pi * 94.5e6 + 2 * pi * val * 1e6);
    }

    cdouble D2 = 2 * pi * 94.5e6;
    cdouble O1 = 2 * pi * 2e6;
    cdouble O2 = 2 * pi * 5.2e6;
    cdouble G1 = 2 * pi * 3e6;
    cdouble G2 = 2 * pi * 3e6;
    cdouble G = G1 + G2;

    std::cout << "Check ratios:\n";
    std::cout << (D2 * (O2 * O2 - O1 * O1) / (4.0 * D2 * D2 + G * G) / (2.0 * pi * 1e6)) << "\n";
    std::cout << (G * (O2 * O2 + O1 * O1) / (4.0 * D2 * D2 + G * G) / (2.0 * pi * 1e6)) << "\n";

    std::vector<double> tot_jumps(D1s.size(), 0.0);

    for (size_t i = 0; i < D1s.size(); ++i) {
        auto start = std::chrono::steady_clock::now();

        MatrixXcd H = MatrixXcd::Zero(size, size);
        H = add_light_coupling_component(H, 1, 0, D1s[i], O1);
        H = add_light_coupling_component(H, 2, 0, D2, O2);

        MatrixXcd L1 = decay_rate(size, 1, 0, G1);
        MatrixXcd L2 = decay_rate(size, 2, 0, G2);

        std::vector<MatrixXcd> L = {L1, L2};

        // Prepare time vector (0 to 10000 us, 200001 points)
        int num_points = 20001;
        VectorXd t = VectorXd::LinSpaced(num_points, 0.0, 1000e-6);
        double dt = t(1) - t(0);
        std::vector<double> t_double_vector(t.data(), t.data() + t.size());
        std::cout << "dt: " << dt << "\n";

        VectorXcd psi0 = VectorXcd::Zero(size);
        psi0(1) = 1.0;

        // Here you need your RNG key or equivalent
        unsigned int key = 42;

        // Call solve function - placeholder
        auto [psi, jumps, expH] = stochastic_schrodinger::solve(t_double_vector, psi0, H, L, key);

        // For now, just simulate jumps = 0 and psi = psi0 (replace with real solve)
        // int jumps = 0;
        // VectorXcd psi = psi0;

        tot_jumps[i] = jumps;

        std::cout << "Total jumps: " << jumps << "\n";
        std::cout << "Psi: " << psi.transpose() << "\n";

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Elapsed time (s) per 1e3 timestamps: " << elapsed.count() * 1e3 / (num_points - 1) << "\n";
    }

    // Output results to console (no plot)
    std::cout << "D1 detuning (MHz), Tot jumps:\n";
    for (size_t i = 0; i < D1s.size(); ++i) {
        double detuning_MHz = (std::real(D1s[i]) - std::real(D2)) / (2 * pi * 1e6);
        std::cout << detuning_MHz << ", " << tot_jumps[i] << "\n";
    }

    return 0;
}
