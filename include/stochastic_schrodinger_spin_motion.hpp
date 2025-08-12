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

MatrixXc propagator(const MatrixXc& H, double dt);

// compute average x and z coordinates
std::pair<double,double> compute_avg_x_z(const VectorXc& psi,
                                         int n_x, int n_z);

using StateCarry = std::tuple<VectorXc,int,long long,long long,long long,VecD,VecD, // psi, scatter_count, ii, step_count, per_log_step, nx, nz
                              std::vector<std::tuple<int, int, double>>, // Lt entries: i,j,val
                              double, MatrixXc, MatrixXc, // G_tot, U, W
                              int, int, // sizes: n_x,n_z
                              double, std::mt19937>; // dt, rng seed
std::pair<StateCarry, void*> step(const StateCarry& carry);

// Full solve: returns (final psi, mean jumps, avg_x, avg_z)
std::tuple<VectorXc,double,VecD,VecD>
solve(const double time_step,
      const long long num_steps,
      const long long per_log_step,
      const VectorXc& psi0,
      const MatrixXc& H,
      const MatrixXc& W,
      const std::vector<std::tuple<int, int, double>>& Lt,
      const double G_tot,
      int n_x, int n_z,
      int num_keys);
} // namespace ss_spin

#endif // STOCHASTIC_SCHRODINGER_SPIN_MOTION_HPP