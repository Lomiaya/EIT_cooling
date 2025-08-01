#ifndef HAMILTONIAN_HPP
#define HAMILTONIAN_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <complex>
#include <vector>
#include <tuple>
#include <map>
#include <string>
#include "define_params.hpp"

namespace hamiltonian {

using namespace Eigen;
using namespace std;
using namespace define_params;

MatrixXcd build_H_zero_freq(const MatrixXcd& H_ground,
                            const MatrixXcd& H_excited,
                            int n_ground_states,
                            int n_excited_states,
                            int n_x_max,
                            int n_z_max,
                            const Params& params);

tuple<int, int, int> from_number_to_tuple(int n, int n_x_max, int n_z_max);

MatrixXcd build_H_light_transition_excite(double Iii,
                                          const Vector3cd& pol,
                                          const Vector3d& vec,
                                          const vector<MatrixXd>& G,
                                          double Isat,
                                          double mass,
                                          double wavelength,
                                          int n_ground_states,
                                          int n_excited_states,
                                          int n_x_max,
                                          int n_z_max,
                                          const Params& params);

MatrixXcd build_H_light_transition_deexci(double Iii,
                                          const Vector3cd& pol,
                                          const Vector3d& vec,
                                          const vector<MatrixXd>& G,
                                          double Isat,
                                          double mass,
                                          double wavelength,
                                          int n_ground_states,
                                          int n_excited_states,
                                          int n_x_max,
                                          int n_z_max,
                                          const Params& params);

MatrixXcd build_H_light_transition(double Iii,
                                   const Vector3cd& pol,
                                   const Vector3d& vec,
                                   const vector<MatrixXd>& G,
                                   double Isat,
                                   double mass,
                                   double wavelength,
                                   int n_ground_states,
                                   int n_excited_states,
                                   int n_x_max,
                                   int n_z_max,
                                   const Params& params,
                                   bool excite);

} // namespace hamiltonian

#endif // HAMILTONIAN_HPP
