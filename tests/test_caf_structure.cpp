#include <iostream>
#include "caf_states.hpp"
#include "define_params.hpp"
#include "state_hamiltonian.hpp"

int main() {
    define_params::States states = caf_states::define_states();
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
    return 0;
}
