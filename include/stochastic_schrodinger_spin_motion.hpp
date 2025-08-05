// stochastic_schrodinger_spin_motion.hpp
#ifndef STOCHASTIC_SCHRODINGER_SPIN_MOTION_HPP
#define STOCHASTIC_SCHRODINGER_SPIN_MOTION_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>
#include <tuple>
#include <complex>

namespace ss_spin {

using Complex = std::complex<double>;
using MatrixXc = Eigen::MatrixXcd;
using VectorXc = Eigen::VectorXcd;
using VecD = Eigen::VectorXd;

// Brute‐force first‐order propagator
MatrixXc brute_force(const std::vector<std::pair<MatrixXc,double>>& Hi_omegas,
                     const MatrixXc& G,
                     double t0, double t1,
                     int n_steps = 1000);

// Magnus second‐order analytic integrator
MatrixXc magnus2_analytic(const std::vector<std::pair<MatrixXc,double>>& Hi_omegas,
                           const MatrixXc& G,
                           double t0, double t1);

// Average position computation
std::pair<double,double> compute_avg_x_z(const VectorXc& psi,
                                          int n_s, int n_x, int n_z);

// Builders for Li entries and translate to full operator
Eigen::VectorXi expand(int idx,
                       const std::array<int,5>& size_list);

// Top‐level build_L
std::vector<std::tuple<int, int, double>> build_L(const std::vector<Eigen::MatrixXd>& G,
                                                  int n_s, int n_x, int n_z,
                                                  int num_nonzero,
                                                  double mass,
                                                  double omega_x,
                                                  double omega_z,
                                                  double wavelength,
                                                  const std::array<double, 3>& B_direction,
                                                  unsigned int seed);

// Time‐step for spin‐motion evolution
using StateCarry = std::tuple<VectorXc,int,int,VecD,VecD, // psi, scatter_count, ii, nx, nz
                              std::vector<std::tuple<int, int, double>>, // Lt entries: i,j,val
                              double, // G_tot
                              int, int, int, int, // sizes: n_g,n_s,n_x,n_z
                              double, std::mt19937>; // dt, rng seed
std::pair<StateCarry, void*> step(const StateCarry& carry,
                                  const std::pair<double,double>& t_pair);

// Full solve: returns (final psi, mean jumps, avg_x, avg_z)
std::tuple<VectorXc,double,VecD,VecD>
solve(const VecD& time,
      const VectorXc& psi0,
      const std::vector<std::pair<MatrixXc,double>>& Ht,
      const std::vector<std::tuple<int, int, double>>& Lt,
      const double G_tot, int n_g,
      int n_s, int n_x, int n_z,
      int num_keys);

} // namespace ss_spin

#endif // STOCHASTIC_SCHRODINGER_SPIN_MOTION_HPP