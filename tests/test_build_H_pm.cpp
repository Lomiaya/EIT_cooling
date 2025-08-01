#include "hamiltonian.hpp"
#include "constants.hpp"
#include "define_params.hpp"
#include <Eigen/Dense>
#include <iostream>

int main() {
    using namespace Eigen;
    using namespace std;
    using namespace hamiltonian;

    double Iii = 2.0;
    double G_tot = 200.0;
    Vector3d pol(1.0, 0.0, 0.0);
    Vector3d vec(0.0, 0.0, 1.0);

    MatrixXd Gp1b(1, 2);
    Gp1b << 0, G_tot / 2;

    MatrixXd Gp0(1, 2);
    Gp0 << G_tot / 2, 0;

    MatrixXd Gp1(1, 2);
    Gp1 << 0, 0;

    vector<MatrixXd> G = {Gp1b, Gp0, Gp1};

    double mass = 59 * parameters::atomic_unit_weight;
    double wavelength = 606e-9;
    double Isat = 1.0;

    int n_ground_states = 2;
    int n_excited_states = 2;
    int n_x_max = 2;
    int n_z_max = 2;

    define_params::Params params;
    params.omega_x = 2 * parameters::pi * 73e3;
    params.omega_z = 2 * parameters::pi * 10e3;

    auto start_time = std::chrono::high_resolution_clock::now();

    MatrixXcd Hp, Hm;

    for (int i=0; i<200; ++i) {
        Hp = build_H_light_transition_excite(Iii, pol, vec, G, Isat, mass, wavelength,
                                                    n_ground_states, n_excited_states,
                                                    n_x_max, n_z_max, params);

        Hm = build_H_light_transition_deexci(Iii, pol, vec, G, Isat, mass, wavelength,
                                                    n_ground_states, n_excited_states,
                                                    n_x_max, n_z_max, params);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "Time taken by function: " << duration << " milliseconds" << std::endl;

    cout << "Hp =\n" << Hp << endl;
    cout << "Hp - Hm^\u2020 =\n" << (Hp - Hm.adjoint()) << endl;

    return 0;
}
