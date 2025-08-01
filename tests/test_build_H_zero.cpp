#include "hamiltonian.hpp"
#include "constants.hpp"
#include "define_params.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <chrono>

int main() {
    using namespace Eigen;
    using namespace std;

    MatrixXcd H_ground(2, 2);
    H_ground << 1, 0,
                0, 1;

    MatrixXcd H_excited(2, 2);
    H_excited << 0, 1,
                 1, 0;

    int n_ground_states = 2;
    int n_excited_states = 2;

    define_params::Params params;
    params.omega_x = 0.01 / parameters::reduced_plancks_constant;
    params.omega_z = 0.1 / parameters::reduced_plancks_constant;
    params.n_x_max = 2;
    params.n_z_max = 2;

    auto start_time = std::chrono::high_resolution_clock::now();

    MatrixXcd H0;

    for (int i=0; i<200; ++i) {
        H0 = hamiltonian::build_H_zero_freq(H_ground, H_excited,
                                            n_ground_states, n_excited_states,
                                            params.n_x_max, params.n_z_max,
                                            params);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "Time taken by function: " << duration << " milliseconds" << std::endl;

    cout << "H0 =\n" << H0 << endl;
    return 0;
}
