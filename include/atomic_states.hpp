// atomic_states.hpp
#ifndef ATOMIC_STATES_HPP
#define ATOMIC_STATES_HPP

#include <array>
#include <vector>
#include <complex>
#include <Eigen/Dense>

namespace atomic_states {

constexpr int n_ground_states = 2;
constexpr int n_excited_states = 1;

using Complex = std::complex<double>;
using Matrix = Eigen::MatrixXcd;

extern const std::array<double, 3> B_direction;

extern const double G1;
extern const double G2;
extern const double G_tot;

extern const std::vector<Eigen::MatrixXd> G; // shape: [3][n_excited_states][n_ground_states]

extern const Matrix H_ground;  // n_ground_states x n_ground_states
extern const Matrix H_excited; // n_excited_states x n_excited_states

extern const double transition_lambda;

} // namespace atomic_states

#endif // ATOMIC_STATES_HPP
