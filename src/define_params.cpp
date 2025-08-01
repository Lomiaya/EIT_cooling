// define_params.cpp
#include "define_params.hpp"
#include "atomic_states.hpp"
#include "constants.hpp"

#include <cmath>

namespace define_params {

Params create_params() {
    Params params;

    params.n_beams = 2;

    int size_D = 40;

    // Detunings
    DoubleVec D1s(size_D);
    for (int i = 0; i < size_D; ++i) {
        double delta = static_cast<double>(i - 20) / 20.0 * 2e6;  // linspace(-2e6,2e6,40)
        double val = 2.0 * parameters::pi * 94.5e6 + 2 * parameters::pi * delta;
        D1s(i) = val;
    }

    double D2 = 2.0 * parameters::pi * 94.5e6 + atomic_states::H_ground(1,1).real();
    DoubleMat D(size_D,params.n_beams);
    int i = 0;
    for (const auto& D1 : D1s) {
        D(i,0) = D1;
        D(i,1) = D2;
        i += 1;
    }
    params.D = D;

    // Intensities
    params.I = DoubleMat(1,params.n_beams);
    params.I << 10.0, 50.0;

    // Polarizations
    params.s = ComplexMat(params.n_beams,3);
    params.s << 1.0, 0.0, 0.0,  0.0, 1.0, 0.0;

    // Wave vectors
    params.k = DoubleMat(params.n_beams,3);
    params.k << 0.0, 0.0, 1.0,  1.0, 0.0, 0.0;

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
