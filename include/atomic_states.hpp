// atomic_states.hpp
#ifndef ATOMIC_STATES_HPP
#define ATOMIC_STATES_HPP

#include <array>
#include <map>
#include <vector>
#include <complex>
#include <Eigen/Dense>
#include "define_params.hpp"

namespace atomic_states {
using Complex = std::complex<double>;
using Matrix = Eigen::MatrixXcd;
using namespace define_params;

States define_states();
Params create_params(const States& states);
} // namespace atomic_states

#endif // ATOMIC_STATES_HPP
