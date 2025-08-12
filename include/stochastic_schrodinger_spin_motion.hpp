// stochastic_schrodinger_spin_motion.hpp
#ifndef STOCHASTIC_SCHRODINGER_SPIN_MOTION_HPP
#define STOCHASTIC_SCHRODINGER_SPIN_MOTION_HPP

#include <Eigen/Dense>
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
using VecD = Eigen::VectorXd;

MatrixXc brute_force(const std::vector<std::pair<MatrixXc,double>>& Hi_omegas,
                     double t0, double t1,
                     int n_steps = 1000);

// Magnus second‐order analytic integrator
MatrixXc magnus2_analytic(const SpectrumMatrix& Hi_omegas,
                          double t0, double t1);

// compute average x and z coordinates
std::pair<double,double> compute_avg_x_z(const VectorXc& psi,
                                         int n_x, int n_z);

using StateCarry = std::tuple<VectorXc,int,long long,long long,long long,VecD,VecD, // psi, scatter_count, ii, step_count, per_log_step, nx, nz
                              std::vector<std::tuple<int, int, double>>, // Lt entries: i,j,val
                              double, // G_tot
                              int, int, // sizes: n_x,n_z
                              double, std::mt19937>; // dt, rng seed
std::pair<StateCarry, void*> step(const StateCarry& carry,
                                  const MatrixXc& expH,
                                  const MatrixXc& Wt);

// Full solve: returns (final psi, mean jumps, avg_x, avg_z)
std::tuple<VectorXc,double,VecD,VecD>
solve(const double time_step,
      const long long num_steps,
      const long long per_log_step,
      const VectorXc& psi0,
      const SpectrumMatrix& H,
      const SpectrumMatrix& W,
      const std::vector<std::tuple<int, int, double>>& Lt,
      const double G_tot,
      int n_x, int n_z,
      int num_keys);
} // namespace ss_spin

#endif // STOCHASTIC_SCHRODINGER_SPIN_MOTION_HPP