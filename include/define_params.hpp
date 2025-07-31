// define_params.hpp
#ifndef DEFINE_PARAMS_HPP
#define DEFINE_PARAMS_HPP

#include <array>
#include <complex>
#include <vector>

namespace define_params {

using Complex = std::complex<double>;
using ComplexVec = std::vector<Complex>;
using ComplexMat = std::vector<std::vector<Complex>>;
using DoubleVec = std::vector<double>;
using DoubleMat = std::vector<std::vector<double>>;
using IntVec = std::vector<int>;

struct Params {
    int n_beams;

    ComplexMat D;  // [n_scans_D][n_beams]
    DoubleMat I;   // [n_scans_I][n_beams]
    DoubleMat s;   // [n_beams][3]
    DoubleMat k;   // [n_beams][3]

    double omega_x;
    double omega_z;

    int n_x_max;
    int n_z_max;

    int n_x_init;
    int n_z_init;

    double mass;
};

Params create_params();

} // namespace define_params

#endif // DEFINE_PARAMS_HPP
