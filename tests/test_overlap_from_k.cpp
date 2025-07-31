#include <iostream>
#include <chrono>
#include "chance_of_jump.hpp"
#include "constants.hpp"

int main() {
    using namespace chance_of_jump;

    double wavelength = 606e-9;
    Vec3 k{0, 0, 1};
    k[0] *= 2 * parameters::pi / wavelength;
    k[1] *= 2 * parameters::pi / wavelength;
    k[2] *= 2 * parameters::pi / wavelength;

    double atomic_unit_weight = 1.66053873e-27;
    double mass = atomic_unit_weight * 59;
    double omega_x = 100e3 * 2 * parameters::pi;
    double omega_z = 10e3 * 2 * parameters::pi;
    std::array<int,3> start{3,3,0};
    std::array<int,3> end{3,3,0};
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i=0; i<20; ++i) {
        Complex c = calculate_overlap_from_k(k, mass, omega_x, omega_z, start, end);
        std::cout << "Overlap: " << c << "\n";
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "Time taken by function: " << duration << " milliseconds" << std::endl;
}
