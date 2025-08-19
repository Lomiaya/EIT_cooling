#include "simulation.hpp"
#include "caf_states.hpp"
#include <fstream>

using namespace caf_states;

void write_to_file(std::string file_name, std::string label, const Eigen::MatrixXd& data) {
    // Open file for writing
    std::ofstream file(file_name, std::ios::out | std::ios::app);
    if (!file.is_open()) {
        std::cerr << "Error: could not open file for writing. " << file_name;
    }

    // Optional: format without brackets, with space separation
    Eigen::IOFormat format(Eigen::FullPrecision, Eigen::DontAlignCols, " ", "\n");

    // Write matrix to file
    file << label << " " << data.format(format) << std::endl;
}

int main() {
    for (int i = 0; i < 16; ++i) {
        std::cout << "Running simulation, B = " << (i-0.01) / 16.0 << "G..." << std::endl;
        States states = define_states((i-0.01) / 16.0); // Vary B-field from -0 to 1 Gauss
        Params params = create_params(states);
        auto [tot_jumps, avg_tempsx, avg_tempsz] = simulation::simulate(params, states, 800001, 20000e-6, 10, 1e6);
        write_to_file("./tot_jumps.txt", "B = " + std::to_string((i-0.01) / 16.0), tot_jumps);
        write_to_file("./avg_tempsx.txt", "B = " + std::to_string((i-0.01) / 16.0), avg_tempsx);
        write_to_file("./avg_tempsz.txt", "B = " + std::to_string((i-0.01) / 16.0), avg_tempsz);
    }
    return 0;
}
