// atomic_states.cpp
#include "caf_states.hpp"
#include "constants.hpp"
#include "define_params.hpp"
#include "half_integer.hpp"
#include "state_hamiltonian.hpp"
#include "hunds_case_b.hpp"

#include <cmath>

namespace caf_states {

using namespace define_params;
using namespace state_hamiltonian;
using namespace half_integer;
using namespace hunds_case_b;

std::vector<HundsCaseB_Rot> X_basis_states() {
    std::vector<HundsCaseB_Rot> basis;
    std::string label = "X";
    int v = 0;
    HalfInteger S = HalfInteger::from_twice(1);
    HalfInteger I = HalfInteger::from_twice(1);
    HalfInteger L = HalfInteger::from_twice(0);
    HalfInteger N = HalfInteger::from_twice(2);
    for (int M_ = -1; M_ <= 1; ++M_) {
        HalfInteger J = HalfInteger::from_twice(1);
        HalfInteger F = HalfInteger(1);
        HundsCaseB_Rot state = {label, v, S, I, L, N, J, F, HalfInteger(M_)};
        basis.push_back(state);
    }
    for (int M_ = -0; M_ <= 0; ++M_) {
        HalfInteger J = HalfInteger::from_twice(1);
        HalfInteger F = HalfInteger(0);
        HundsCaseB_Rot state = {label, v, S, I, L, N, J, F, HalfInteger(M_)};
        basis.push_back(state);
    }
    for (int M_ = -1; M_ <= 1; ++M_) {
        HalfInteger J = HalfInteger::from_twice(3);
        HalfInteger F = HalfInteger(1);
        HundsCaseB_Rot state = {label, v, S, I, L, N, J, F, HalfInteger(M_)};
        basis.push_back(state);
    }
    for (int M_ = -2; M_ <= 2; ++M_) {
        HalfInteger J = HalfInteger::from_twice(3);
        HalfInteger F = HalfInteger(2);
        HundsCaseB_Rot state = {label, v, S, I, L, N, J, F, HalfInteger(M_)};
        basis.push_back(state);
    }
    return basis;
}

std::vector<HundsCaseB_Rot> A_basis_states() {
    std::vector<HundsCaseB_Rot> basis;
    std::string label = "A";
    int v = 0;
    HalfInteger S = HalfInteger::from_twice(1);
    HalfInteger I = HalfInteger::from_twice(1);
    HalfInteger L = HalfInteger::from_twice(2); // only half of the correct parity state...
    HalfInteger N = HalfInteger::from_twice(0);
    for (int M_ = -1; M_ <= 1; ++M_) {
        HalfInteger J = HalfInteger::from_twice(1);
        HalfInteger F = HalfInteger(1);
        HundsCaseB_Rot state = {label, v, S, I, L, N, J, F, HalfInteger(M_)};
        basis.push_back(state);
    }
    for (int M_ = -0; M_ <= 0; ++M_) {
        HalfInteger J = HalfInteger::from_twice(1);
        HalfInteger F = HalfInteger(0);
        HundsCaseB_Rot state = {label, v, S, I, L, N, J, F, HalfInteger(M_)};
        basis.push_back(state);
    }
    return basis;
}

States define_states() {
    int n_excited_states = 4;
    int n_ground_states = 12;

    double G_tot = 2 * parameters::pi * 8.3e6;

    // ground state
    auto X_state_basis = X_basis_states();
    auto X_hamiltonian = StateHamiltonian(X_state_basis);

    // parameters; B, D, gamma, bF, c
    std::vector<double> molecular_params = {10303.988 * 1e6, 0.014060 * 1e6, 39.65891 * 1e6, 122.5569 * 1e6, 40.1190 * 1e6};
    std::vector<double (*)(const HundsCaseB_Rot&, const HundsCaseB_Rot&)> ops = 
    {Rotation, RotationDistortion, SpinRotation, Hyperfine_IS, Hyperfine_Dipolar};

    X_hamiltonian.add_operators_to_matrix(molecular_params, ops);

    Matrix H_ground = X_hamiltonian.get_matrix();

    // excited state
    auto A_state_basis = A_basis_states();
    auto A_hamiltonian = StateHamiltonian(A_state_basis);

    Matrix H_excited = (Matrix(4, 4) << 
        Complex(0.0), Complex(0.0), Complex(0.0), Complex(0.0),
        Complex(0.0), Complex(0.0), Complex(0.0), Complex(0.0),
        Complex(0.0), Complex(0.0), Complex(0.0), Complex(0.0),
        Complex(0.0), Complex(0.0), Complex(0.0), Complex(2 * parameters::pi * 1e6)
    ).finished();

    auto G_pre = calculate_TDM(X_hamiltonian, A_hamiltonian, TDM);

    std::vector<Eigen::MatrixXd> G = { 2 * G_tot * G_pre[0], 2 * G_tot * G_pre[1], 2 * G_tot * G_pre[2] };

    std::array<double, 3> B_direction = {std::sqrt(2.0), 0.0, std::sqrt(2.0)};

    double transition_lambda = 606e-9;

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

// haven't really done params!
Params create_params(const States& states) {
    Params params;

    auto energies = get_eigenvalues(states.H_ground);

    params.n_beams = 2;

    int size_D = 20;

    // Detunings
    DoubleVec D1s(size_D);
    for (int i = 0; i < size_D; ++i) {
        double delta = (static_cast<double>(i) - 0.0) / 20.0 * 0.2e6;
        double val = 2.0 * parameters::pi * 94.5e6 + 2.0 * parameters::pi * delta;
        D1s(i) = val;
    }

    double D2 = 2.0 * parameters::pi * 94.5e6 - states.H_ground(1,1).real(); // this is blue detuning!
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
    params.k << std::sqrt(2.0), 0.0, -std::sqrt(2.0),  std::sqrt(2.0), 0.0, std::sqrt(2.0);

    params.omega_x = 2 * parameters::pi * 73e3;
    params.omega_z = 2 * parameters::pi * 10e3;

    params.n_x_max = 5;
    params.n_z_max = 1;

    params.n_x_init = 2;
    params.n_z_init = 0;

    params.mass = 87 * parameters::atomic_unit_weight;

    return params;
}

} // namespace caf_states
