#include "simulation.hpp"
#include "atomic_states.hpp"

using namespace atomic_states;

int main() {
    States states = define_states();
    Params params = create_params(states);
    auto _ = simulation::simulate(params, states, 100001, 2000e-6);
    return 0;
}
