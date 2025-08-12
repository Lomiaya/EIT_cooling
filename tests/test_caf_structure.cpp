#include <iostream>
#include "caf_states.hpp"
#include "define_params.hpp"
#include "state_hamiltonian.hpp"
#include "hunds_case_b.hpp"

int main() {
    define_params::States states = caf_states::define_states(0.0);
    auto X_matrix = states.H_ground;
    auto energies = state_hamiltonian::get_eigenvalues(states.H_ground);
    std::cout << X_matrix << std::endl;
    for (auto& energy: energies) {
        std::cout << energy << std::endl;
    }
    auto G = states.G;
    for (auto& Gi: G) {
        std::cout << Gi / states.G_tot << std::endl;
    }
    auto X_states = caf_states::X_basis_states();
    auto state_a = X_states[3];
    auto state_b = X_states[3];
    std::cout << hunds_case_b::SpinRotation(state_a, state_b) << std::endl;
    return 0;
}
