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

int from_tuple_to_number(int s, int x, int z, const Params& p);
tuple<int, int, int> from_number_to_tuple(int n, const Params& p);

std::tuple<MatrixXd, MatrixXd, MatrixXd> simulate(const Params& params, const States& def_states);

} // namespace simulation

#endif
