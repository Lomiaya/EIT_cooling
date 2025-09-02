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

/**
 * @brief Define simulation parameters
 * @param n_beams int, Number of beams
 * @param D Detuning, [n_scans_D][n_beams] in rad/s
 * @param I Intensity, [n_scans_I][n_beams], W/m^2
 * @param s Polarization, [n_beams][3], normalized [sigma-, pi, sigma+] for each beam
 * @param k K vector, [n_beams][3], normalized [x, y, z], can't have y component
 * @param omega_x Trap frequency in x direction, rad/s
 * @param omega_y Trap frequency in y direction, rad/s
 * @param omega_z Trap frequency in z direction, rad/s
 * @param trap_depth Trap depth in rad/s, used to calculate differential trap shift.
 * @param n_x_max Max motional number in x direction considered in the siumuation.
 * @param n_y_max Max motional number in y direction considered in the siumuation.
 * @param n_z_max Max motional number in z direction considered in the siumuation.
 * @param n_x_init Initial motion number in x direction.
 * @param n_y_init Initial motion number in y direction.
 * @param n_z_init Initial motion number in z direction.
 * @param do_2d_sim Whether to do 2D simulation, if true, only x and z direction is considered.
 * @param mass Mass in kg.
 *  **/
struct Params {
    int n_beams;

    DoubleMat D;  // [n_scans_D][n_beams], rad/s
    DoubleMat I;   // [n_scans_I][n_beams], W/m^2
    ComplexMat s;   // [n_beams][3], [sigma-, pi, sigma+] for each beam
    DoubleMat k;   // [n_beams][3], in the unit of lambda/2*pi, can't have y component

    double omega_x; // Trap frequency, rad/s
    double omega_y = 0;
    double omega_z;
    double trap_depth; // in rad/s, used to calculate differential trap shift.

    int n_x_max; // Max motional number
    int n_y_max = 1;
    int n_z_max;

    int n_x_init;
    int n_y_init = 0;
    int n_z_init;

    bool do_2d_sim = true;

    double mass;    // kg
};

/**
 * @brief Define quantum states
 * @param n_ground_states Number of ground states
 * @param n_excited_states Number of excited states
 * @param B_direction Unit vectore representing the direction of the B field [x, y, z]
 * @param G_tot Excited state scatering rate in rad/s
 * @param G Dipole matrix in the shape: [3][n_excited_states][n_ground_states], in order of [sigma-, pi, sigma+]
 * @param H_ground Ground state Hamiltonian, [n_ground_states][n_ground_states]
 * @param H_excited Excited state Hamiltonian, [n_excites_states][n_excites_states]
 * @param H_ground_stark Ground state Stark shift, n_ground_states x n_ground_states, relative, set to 1 for specified wx, wz, off-diagonal terms will be neglected
 * @param H_excited_stark Excites state Stark shift, n_excited_states x n_excited_states, relative, set to 1 for specified wx, wz, off-diagonal terms will be neglected
 * @param transition_lambda wavelength of the transition light, in m
 */
struct States {
    int n_ground_states;
    int n_excited_states;

    std::array<double, 3> B_direction; // [x, y, z] unit vector

    double G_tot; // Gamma, in rad/s

    std::vector<Eigen::MatrixXd> G; // shape: [3][n_excited_states][n_ground_states]

    ComplexMat H_ground;  // n_ground_states x n_ground_states
    ComplexMat H_excited; // n_excited_states x n_excited_states
    ComplexMat H_ground_stark; // n_ground_states x n_ground_states, relative, set to 1 for specified wx, wz, off-diagonal terms will be neglected
    ComplexMat H_excited_stark; // n_excited_states x n_excited_states, relative, set to 1 for specified wx, wz, off-diagonal terms will be neglected

    double transition_lambda; // wavelength of the transition light, in m
};

} // namespace define_params

#endif // DEFINE_PARAMS_HPP
