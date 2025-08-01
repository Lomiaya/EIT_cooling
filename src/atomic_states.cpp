// atomic_states.cpp
#include "atomic_states.hpp"
#include "constants.hpp"

#include <cmath>

namespace atomic_states {

const std::array<double, 3> B_direction = {0.0, 0.0, 1.0};

const double G1 = 0.4 * 2 * parameters::pi * 6e6;
const double G2 = 0.6 * 2 * parameters::pi * 6e6;

const double G_tot = 2 * parameters::pi * 6e6;

// Transition amplitudes: 3 × n_excited × n_ground
const std::vector<Eigen::MatrixXd> G = {
    // σ⁻
    (Eigen::MatrixXd(1, 2) << 0.0, G2).finished(),
    // π
    (Eigen::MatrixXd(1, 2) << G1, 0.0).finished(),
    // σ⁺
    (Eigen::MatrixXd(1, 2) << 0.0, 0.0).finished()
};

// Ground state Hamiltonian: 2 × 2
const Matrix H_ground = (Matrix(2, 2) << 
    Complex(0.0), Complex(0.0),
    Complex(0.0), Complex(2 * parameters::pi * 5e6)
).finished();

// Excited state Hamiltonian: 1 × 1
const Matrix H_excited = (Matrix(1, 1) << 
    Complex(0.0)
).finished();

const double transition_lambda = 780e-9;

} // namespace atomic_states
