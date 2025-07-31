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
const std::vector<Matrix> G = {
    // σ⁻
    Matrix{
        {Complex(0.0), Complex(G2)}
    },
    // π
    Matrix{
        {Complex(G1), Complex(0.0)}
    },
    // σ⁺
    Matrix{
        {Complex(0.0), Complex(0.0)}
    }
};

const Matrix H_ground = {
    {Complex(0.0), Complex(0.0)},
    {Complex(0.0), Complex(2 * parameters::pi * 500e6)}
};

const Matrix H_excited = {
    {Complex(0.0)}
};

const double transition_lambda = 780e-9;

} // namespace atomic_states
