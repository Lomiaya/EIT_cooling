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
using SpectrumMatrix = std::vector<std::pair<ComplexMat, double>>;

tuple<int, int, int> from_number_to_tuple(int n, int n_x_max, int n_z_max);
// Diagonalizes the Hamiltonian and scattering matrices in the States struct
States diagonalize_hamiltonian(const States& states);
// Defines the partition-state Hamiltonian vector
DoubleVec define_partition_hamiltonian(const ComplexMat& H, int n_x_max, int n_z_max, const Params& params);
SpectrumMatrix multiply(const SpectrumMatrix& A,
                        const SpectrumMatrix& B);
SpectrumMatrix multiply(const SpectrumMatrix& A,
                        const SpectrumMatrix& B,
                        double threshold);
SpectrumMatrix multiply(const ComplexMat& A,
                        const SpectrumMatrix& B);
SpectrumMatrix multiply(const SpectrumMatrix& A,
                        const ComplexMat& B);
SpectrumMatrix multiply(double c, const SpectrumMatrix& A);
SpectrumMatrix multiply(Complex c, const SpectrumMatrix& A);
SpectrumMatrix adjoint(const SpectrumMatrix& A);
SpectrumMatrix addition(const SpectrumMatrix& A,
                        const SpectrumMatrix& B);
SpectrumMatrix cleanup(const SpectrumMatrix& A);
ComplexMat evaluate(const SpectrumMatrix& A, double t);
SpectrumMatrix low_pass_filter(const SpectrumMatrix& A, double threshold);
// Defines V_+(f,l) matrix, of shape n'_excited_states x n'_ground_states
SpectrumMatrix define_V_plus(const States& states,
                             const Params& params,
                             int f,
                             int l,
                             int I_index,
                             int D_index,
                             DoubleVec& H_ground_diag);
// Defines V_-(f) matrix, of shape n'_ground_states x n'_excited_states
SpectrumMatrix define_V_minus(const States& states,
                              const Params& params,
                              int f,
                              int I_index,
                              int D_index,
                              DoubleVec& H_ground_diag);
// Sums over the V_-(f) matrices, wrapper, of shape n'_ground_states x n'_excited_states
SpectrumMatrix V_minus(const States& states,
                       const Params& params,
                       int I_index,
                       int D_index,
                       DoubleVec& H_ground_diag);
// Build W matrix, of shape n'_excited_states x n'_ground_states
SpectrumMatrix build_W(const States& states,
                       const Params& params,
                       int I_index,
                       int D_index);
// Build L (each entry: ground state, excited state, value)
// Do NOT build ground state L, it is not efficient.
std::vector<std::tuple<int, int, double>> build_L(const std::vector<Eigen::MatrixXd>& G,
                                                  int n_x, int n_z,
                                                  int num_nonzero,
                                                  double mass,
                                                  double omega_x,
                                                  double omega_z,
                                                  double wavelength,
                                                  const std::array<double, 3>& B_direction,
                                                  unsigned int seed);
// Build the Hamiltonian matrix, of shape n'_ground_states x n'_ground_states
SpectrumMatrix build_H(const States& states,
                       const Params& params,
                       int I_index,
                       int D_index,
                       SpectrumMatrix& W);

} // namespace hamiltonian

#endif // HAMILTONIAN_HPP
