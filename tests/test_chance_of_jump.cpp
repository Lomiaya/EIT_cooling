#include <iostream>
#include <chrono>
#include "chance_of_jump.hpp"
#include "constants.hpp"

int main() {
    using namespace chance_of_jump;

    int seed = 42;
    int np = 1;
    double mass = parameters::atomic_unit_weight * 59;
    double omega_x = 100e3 * 2 * parameters::pi;
    double omega_z = 10e3 * 2 * parameters::pi;
    double wavelength = 606e-9;
    Vec3 B_direction{0,0,1};
    std::array<int,3> start{3,3,0};
    std::array<int,3> end{3,3,2};
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i=0; i<20; ++i) {
        int num_repeats = 500;
        auto results = calculate_chance_of_jump(seed, np, mass, omega_x, omega_z, wavelength, B_direction, start, end, num_repeats);
        double res = 0;
        for (const auto& c : results) {
            res += c;
        }
        std::cout << "Overlap: " << res / num_repeats << "\n";
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "Time taken by function: " << duration << " milliseconds" << std::endl;
}
