// atomic_states.cpp
#include "atomic_states.hpp"
#include "constants.hpp"
#include "define_params.hpp"

#include <cmath>

namespace atomic_states {

using namespace define_params;

States define_states() {
    int n_excited_states = 1;
    int n_ground_states = 2;

    std::array<double, 3> B_direction = {0.0, 0.0, 1.0};

    double G1 = 0.6 * 2 * parameters::pi * 6e6;
    double G2 = 0.4 * 2 * parameters::pi * 6e6;

    double G_tot = 2 * parameters::pi * 6e6;

    // Transition amplitudes: 3 × n_excited × n_ground
    std::vector<Eigen::MatrixXd> G = {
        // σ⁻
        (Eigen::MatrixXd(1, 2) << 0.0, G2).finished(),
        // π
        (Eigen::MatrixXd(1, 2) << G1, 0.0).finished(),
        // σ⁺
        (Eigen::MatrixXd(1, 2) << 0.0, 0.0).finished()
    };

    // Ground state Hamiltonian: 2 × 2
    Matrix H_ground = (Matrix(2, 2) << 
        Complex(0.0), Complex(0.0),
        Complex(0.0), Complex(2 * parameters::pi * 5e6)
    ).finished();

    // Excited state Hamiltonian: 1 × 1
    Matrix H_excited = (Matrix(1, 1) << 
        Complex(0.0)
    ).finished();

    double transition_lambda = 780e-9;

    States states;
    states.B_direction = B_direction;
    states.G = G;
    states.G_tot = G_tot;
    states.H_excited = H_excited;
    states.H_ground = H_ground;
    states.n_excited_states = n_excited_states;
    states.n_ground_states = n_ground_states;
    states.transition_lambda = transition_lambda;
    return states;
}

Params create_params(const States& states) {
    Params params;

    params.n_beams = 2;

    int size_D = 10;

    // Detunings
    DoubleVec D1s(size_D);
    for (int i = 0; i < size_D; ++i) {
        double delta = (static_cast<double>(i) - 0.0) / 10.0 * 0.5e6;
        double val = - 2.0 * parameters::pi * 94.5e6 - 2.0 * parameters::pi * delta;
        D1s(i) = val;
    }

    double D2 = - 2.0 * parameters::pi * 94.5e6 - states.H_ground(1,1).real(); // this is blue detuning!
    DoubleMat D(size_D,params.n_beams);
    int i = 0;
    for (const auto& D1 : D1s) {
        D(i,0) = D1;
        D(i,1) = D2;
        i += 1;
    }
    params.D = D; // this is blue detuning!

    // Intensities
    params.I = DoubleMat(1,params.n_beams);
    params.I << 25.0, 400.0;

    // Polarizations
    params.s = ComplexMat(params.n_beams,3);
    params.s << 0.0, 1.0, 0.0,  1.0, 0.0, 0.0;

    // Wave vectors
    params.k = DoubleMat(params.n_beams,3);
    params.k << 1.0, 0.0, 0.0,  0.0, 0.0, 1.0;

    params.omega_x = 2 * parameters::pi * 73e3;
    params.omega_z = 2 * parameters::pi * 10e3;

    params.n_x_max = 5;
    params.n_z_max = 1;

    params.n_x_init = 2;
    params.n_z_init = 0;

    params.mass = 87 * parameters::atomic_unit_weight;

    return params;
}

} // namespace atomic_states
