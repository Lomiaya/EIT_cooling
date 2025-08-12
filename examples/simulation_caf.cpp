#include "simulation.hpp"
#include "caf_states.hpp"
#include <fstream>

using namespace caf_states;

void write_to_file(std::string file_name, const Eigen::MatrixXd& data) {
    // Open file for writing
    std::ofstream file(file_name);
    if (!file.is_open()) {
        std::cerr << "Error: could not open file for writing. " << file_name;
    }

    // Optional: format without brackets, with space separation
    Eigen::IOFormat format(Eigen::FullPrecision, Eigen::DontAlignCols, " ", "\n");

    // Write matrix to file
    file << data.format(format);
}

int main() {
    for (int i = 0; i < 10; ++i) {
        std::cout << "Running simulation, B = " << (i-2.0) / 10.0 << "G..." << std::endl;
        States states = define_states((i-2.0) / 10.0); // Vary B-field from -0.2 to 0.7 Gauss
        Params params = create_params(states);
        auto [tot_jumps, avg_tempsx, avg_tempsz] = simulation::simulate(params, states, 500001, 5000e-6, 10, 1e6);
        write_to_file("./tot_jumps.txt", tot_jumps);
        write_to_file("./avg_tempsx.txt", avg_tempsx);
        write_to_file("./avg_tempsz.txt", avg_tempsz);
    }
    return 0;
}
