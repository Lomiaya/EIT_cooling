// simulation.hpp
#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "define_params.hpp"
#include <Eigen/Dense>
#include <tuple>
using namespace Eigen;

namespace simulation {

using namespace define_params;
using namespace std;

// Actual simulation function wrapper
std::tuple<MatrixXd, MatrixXd, MatrixXd> simulate(const Params& params, const States& def_states, const int N, const double t_0, const int num_keys);

} // namespace simulation

#endif
