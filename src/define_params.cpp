// define_params.cpp
#include "define_params.hpp"
#include "atomic_states.hpp"
#include "constants.hpp"

#include <cmath>

namespace define_params {

Params create_params() {
    Params params;

    params.n_beams = 2;

    // Detunings
    ComplexVec D1s;
    for (int i = 0; i < 20; ++i) {
        double delta = static_cast<double>(i) / 19.0 * 2e6;  // linspace(0,2e6,20)
        Complex val = Complex(2 * parameters::pi * 94.5e6 + 2 * parameters::pi * delta, 0.0);
        D1s.push_back(val);
    }

    Complex D2 = Complex(2 * parameters::pi * 94.5e6, 0.0) + atomic_states::H_ground[1][1];
    ComplexMat D;
    for (const auto& D1 : D1s) {
        D.push_back({D1, D2});
    }
    params.D = D;

    // Intensities
    params.I = {
        {10.0, 50.0}
    };

    // Polarizations
    params.s = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0} // do not let in y-direction; preserved as per note
    };

    // Wave vectors
    params.k = {
        {0.0, 0.0, 1.0},
        {1.0, 0.0, 0.0}
    };

    params.omega_x = 2 * parameters::pi * 73e3;
    params.omega_z = 2 * parameters::pi * 10e3;

    params.n_x_max = 1;
    params.n_z_max = 4;

    params.n_x_init = 0;
    params.n_z_init = 3;

    params.mass = 87 * parameters::atomic_unit_weight;

    return params;
}

} // namespace define_params
