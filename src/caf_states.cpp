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
    std::vector<double> molecular_params = {10303.988 * 1e6, 0.014060 * 1e6, 39.65891 * 1e6, 122.5569 * 1e6, 40.1190 * 1e6 / 3, 
                        (parameters::electron_spin_g_factor * parameters::bohr_magneton / parameters::plancks_constant) * 0e-4}; // B-field in Gauss

    for (auto& param: molecular_params) {
        param *= 2 * parameters::pi;
    }
    std::vector<double (*)(const HundsCaseB_Rot&, const HundsCaseB_Rot&)> ops = 
    {Rotation, RotationDistortion, SpinRotation, Hyperfine_IS, Hyperfine_Dipolar, Zeeman_z};

    X_hamiltonian.add_operators_to_matrix(molecular_params, ops);

    Matrix H_ground = X_hamiltonian.get_matrix();
    auto energies_lowest = get_eigenvalues(H_ground)[0];
    for (int i = 0; i < n_ground_states; ++i) {
        H_ground(i,i) -= energies_lowest;
    }

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

    std::array<double, 3> B_direction = {0.0, 0.0, 1.0};

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
    auto gamma = states.G_tot;

    params.n_beams = 4;

    int size_D = 20;

    // Detunings

    double D1 = gamma * 4.5 - energies[7]; // 4.5gamma detuning from F=2

    DoubleVec D2s(size_D);
    for (int i = 0; i < size_D; ++i) {
        double delta = (static_cast<double>(i) - 10.0) / 10.0 * 400e3;
        double val = 2.0 * parameters::pi * delta + gamma * 4.5 - energies[4]; // 4.5gamma detuning from F=1+
        D2s(i) = val;
    }
    double D3 = gamma * 4.5 - energies[3]; // 4.5gamma detuning from F=0

    DoubleVec D4s(size_D);
    for (int i = 0; i < size_D; ++i) {
        double delta = (static_cast<double>(i) - 10.0) / 10.0 * 400e3;
        double val = 2.0 * parameters::pi * delta + gamma * 4.5 - energies[7];
        D4s(i) = val;
    }
    DoubleMat D(size_D,params.n_beams);
    for (int i = 0; i < size_D; ++i) {
        D(i,0) = D1;
        D(i,1) = D2s[i];
        D(i,2) = D3;
        D(i,3) = D4s[i];
    }
    params.D = D; // this is blue detuning!

    // Intensities
    params.I = DoubleMat(1,params.n_beams);
    // params.I << 400.0, 25.0, 100.0, 25.0;
    params.I << 1000.0, 50.0, 200.0, 50.0;

    // Polarizations
    params.s = ComplexMat(params.n_beams,3);
    params.s << 1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  0.0, 1.0, 0.0,  0.0, 0.0, 1.0;

    // Wave vectors
    params.k = DoubleMat(params.n_beams,3);
    params.k << 0.0, 0.0, 1.0,  1.0, 0.0, 0.0,  -1.0, 0.0, 0.0,  0.0, 0.0, -1.0;

    params.omega_x = 2 * parameters::pi * 100e3;
    params.omega_z = 2 * parameters::pi * 15e3;

    params.n_x_max = 5;
    params.n_z_max = 1;

    params.n_x_init = 2;
    params.n_z_init = 0;

    params.mass = 59 * parameters::atomic_unit_weight;

    return params;
}

} // namespace caf_states
