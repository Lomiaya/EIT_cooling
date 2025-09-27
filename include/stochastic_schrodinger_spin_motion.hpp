// stochastic_schrodinger_spin_motion.hpp
#ifndef STOCHASTIC_SCHRODINGER_SPIN_MOTION_HPP
#define STOCHASTIC_SCHRODINGER_SPIN_MOTION_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>
#include <tuple>
#include <complex>
#include <hamiltonian.hpp>

namespace ss_spin {

using namespace hamiltonian;
using Complex = std::complex<double>;
using MatrixXc = Eigen::MatrixXcd;
using VectorXc = Eigen::VectorXcd;
using SparseMat = SparseMatrix<complex<double>>;
using VecD = Eigen::VectorXd;

SparseMat brute_force(const SpectrumMatrix& Hi_omegas,
                     double t0, double t1,
                     int n_steps = 1000);

// // Magnus second‐order analytic integrator
// MatrixXc magnus2_analytic(const SpectrumMatrix& Hi_omegas,
//                           double t0, double t1);

// compute average x and z coordinates
std::pair<double,double> compute_avg_x_z(const VectorXc& psi,
                                         int n_x, int n_z);

std::tuple<double,double,double> compute_avg_x_y_z(const VectorXc& psi, int n_x, int n_y, int n_z);

using StateCarry = std::tuple<VectorXc,int,long long,long long,long long,VecD,VecD,VecD, // psi, scatter_count, ii, step_count, per_log_step, nx, ny, nz
                              std::vector<std::tuple<int, int, double>>, // Lt entries: i,j,val
                              double, // G_tot
                              int, int, int, bool, // sizes: n_x,n_y,n_z; is_2d_sim
                              double, std::mt19937>; // dt, rng seed
std::pair<StateCarry, void*> step(const StateCarry& carry,
                                  const SparseMat& expH,
                                  const SparseMat& Wt);

// Full solve: returns (final psi, mean jumps, avg_x, avg_z)
std::tuple<VectorXc,double,VecD,VecD,VecD>
solve(const double time_step,
      const long long num_steps,
      const long long per_log_step,
      const VectorXc& psi0,
      const SpectrumMatrix& H,
      const SpectrumMatrix& W,
      const std::vector<std::tuple<int, int, double>>& Lt,
      const double G_tot,
      int n_x, int n_y, int n_z,
      bool is_2d_sim,
      int num_keys,
      double low_pass_threshold);
} // namespace ss_spin

#endif // STOCHASTIC_SCHRODINGER_SPIN_MOTION_HPP