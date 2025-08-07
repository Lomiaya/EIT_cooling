// atomic_states.hpp
#ifndef CAF_STATES_HPP
#define CAF_STATES_HPP

#include <array>
#include <map>
#include <vector>
#include <complex>
#include <Eigen/Dense>
#include "define_params.hpp"
#include "hunds_case_b.hpp"

namespace caf_states {
using Complex = std::complex<double>;
using Matrix = Eigen::MatrixXcd;
using namespace define_params;

std::vector<hunds_case_b::HundsCaseB_Rot> X_basis_states();
std::vector<hunds_case_b::HundsCaseB_Rot> A_basis_states();
States define_states();
Params create_params(const States& states);
} // namespace caf_states

#endif // CAF_STATES_HPP
