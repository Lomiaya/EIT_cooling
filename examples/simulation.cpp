#include "simulation.hpp"
#include "atomic_states.hpp"

using namespace atomic_states;

int main() {
    States states = define_states();
    Params params = create_params(states);
    auto _ = simulation::simulate(params, states, 200001, 5000e-6, 10, 1e6);
    return 0;
}
