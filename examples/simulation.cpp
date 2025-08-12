#include "simulation.hpp"
#include "atomic_states.hpp"

using namespace atomic_states;

int main() {
    States states = define_states();
    Params params = create_params(states);
    auto _ = simulation::simulate(params, states, 501, 5e-6, 10, 1e6); // we still have bugs: if n_max is not zero, would lead to jump at every timestamp.
    return 0;
}
