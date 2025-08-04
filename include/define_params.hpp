// define_params.hpp
#ifndef DEFINE_PARAMS_HPP
#define DEFINE_PARAMS_HPP

#include <array>
#include <complex>
#include <vector>
#include <Eigen/Dense>

namespace define_params {

using Complex = std::complex<double>;
using ComplexVec = Eigen::VectorXcd;
using ComplexMat = Eigen::MatrixXcd;
using DoubleVec = Eigen::VectorXd;
using DoubleMat = Eigen::MatrixXd;

struct Params {
    int n_beams;

    DoubleMat D;  // [n_scans_D][n_beams]
    DoubleMat I;   // [n_scans_I][n_beams]
    ComplexMat s;   // [n_beams][3]
    DoubleMat k;   // [n_beams][3]

    double omega_x;
    double omega_z;

    int n_x_max;
    int n_z_max;

    int n_x_init;
    int n_z_init;

    double mass;
};

struct States {
    int n_ground_states;
    int n_excited_states;

    std::array<double, 3> B_direction;

    double G_tot;

    std::vector<Eigen::MatrixXd> G; // shape: [3][n_excited_states][n_ground_states]

    ComplexMat H_ground;  // n_ground_states x n_ground_states
    ComplexMat H_excited; // n_excited_states x n_excited_states

    double transition_lambda;
};

} // namespace define_params

#endif // DEFINE_PARAMS_HPP
