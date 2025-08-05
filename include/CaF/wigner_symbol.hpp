// wigner_symbol.hpp
#ifndef WIGNER_SYMBOL_HPP
#define WIGNER_SYMBOL_HPP

#include "half_integer.hpp"
#include <cmath>
#include <numeric>
#include <vector>

namespace wigner_symbol {
using namespace half_integer;
// Delta coefficient
double Delta(HalfInteger a, HalfInteger b, HalfInteger c);

// Wigner 3j symbol
double wigner3j(HalfInteger j1, HalfInteger j2, HalfInteger j3,
                HalfInteger m1, HalfInteger m2, HalfInteger m3);

// Wigner 6j symbol
double wigner6j(HalfInteger j1, HalfInteger j2, HalfInteger j3,
                HalfInteger J1, HalfInteger J2, HalfInteger J3);

// Wigner 9j symbol
double wigner9j(HalfInteger j1, HalfInteger j2, HalfInteger j3,
                HalfInteger j4, HalfInteger j5, HalfInteger j6,
                HalfInteger j7, HalfInteger j8, HalfInteger j9);
} // namespace wigner_symbol

#endif // WIGNER_SYMBOL_HPP
