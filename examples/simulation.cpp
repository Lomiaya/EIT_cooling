#include "simulation.hpp"
#include "atomic_states.hpp"

using namespace atomic_states;

int main() {
    States states = define_states();
    Params params = create_params(states);
    auto _ = simulation::simulate(params, states, 5000001, 5000e-6, 10);
    return 0;
}
