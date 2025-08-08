// wigner_symbol.cpp
#include "wigner_symbol.hpp"
#include <iostream>

namespace wigner_symbol {
using namespace half_integer;
inline double Delta(HalfInteger a, HalfInteger b, HalfInteger c) {
    return factorial(a + b - c)
         * factorial(a - b + c)
         * factorial(-a + b + c)
         / factorial(a + b + c + HalfInteger(1));
}

// helper x function
inline double x_term(HalfInteger t, HalfInteger j1, HalfInteger j2, HalfInteger j3,
                     HalfInteger m1, HalfInteger m2) {
    return factorial(t)
         * factorial(t + j3 - j2 + m1)
         * factorial(t + j3 - j1 - m2)
         * factorial(j1 + j2 - j3 - t)
         * factorial(j1 - m1 - t)
         * factorial(j2 + m2 - t);
}

double wigner3j(HalfInteger j1, HalfInteger j2, HalfInteger j3,
                HalfInteger m1, HalfInteger m2, HalfInteger m3) {
    double value = 0.0;
    if ((-j1 <= m1 && m1 <= j1) &&
        (-j2 <= m2 && m2 <= j2) &&
        (-j3 <= m3 && m3 <= j3) &&
        (m1 + m2 == -m3) &&
        (abs(j1 - j2) <= j3 && j3 <= j1 + j2)) {
        int t_min = ceiling(max(HalfInteger(0), max(j2 - j3 - m1, j1 - j3 + m2)));
        int t_max = floor(min(j1 + j2 - j3, min(j1 - m1, j2 + m2)));
        if (t_min <= t_max) {
            double sum = 0.0;
            for (int ti = t_min; ti <= t_max; ++ti) {
                double term = ((ti % 2) == 0 ? 1 : -1) / x_term(HalfInteger(ti), j1, j2, j3, m1, m2);
                sum += term;
            }
            int parity_ = int(j1 - j2 - m3);
            int parity = ((parity_ % 2) == 0 ? 1 : -1);
            double multiplier = Delta(j1, j2, j3) * factorial(j1+m1) * factorial(j1-m1) * factorial(j2+m2) * factorial(j2-m2) * factorial(j3+m3) * factorial(j3-m3);
            value = sum * parity
                  * sqrt(multiplier);
        }
    }
    return value;
}

// helper f term
inline double f_term(HalfInteger t,
    HalfInteger j1j2j3,
    HalfInteger j1J2J3,
    HalfInteger J1j2J3,
    HalfInteger J1J2j3,
    HalfInteger j1j2J1J2,
    HalfInteger j2j3J2J3,
    HalfInteger j3j1J3J1) {
    return factorial(t - j1j2j3)
         * factorial(t - j1J2J3)
         * factorial(t - J1j2J3)
         * factorial(t - J1J2j3)
         * factorial(j1j2J1J2 - t)
         * factorial(j2j3J2J3 - t)
         * factorial(j3j1J3J1 - t);
}

double wigner6j(HalfInteger j1, HalfInteger j2, HalfInteger j3,
                HalfInteger J1, HalfInteger J2, HalfInteger J3) {
    double value = 0.0;
    // triangular conditions
    if ((abs(j1 - j2) <= j3 && j3 <= j1 + j2) &&
        (abs(j1 - J2) <= J3 && J3 <= j1 + J2) &&
        (abs(J1 - j2) <= J3 && J3 <= J1 + j2) &&
        (abs(J1 - J2) <= j3 && j3 <= J1 + J2)) {
        // sums of js
        HalfInteger j1j2j3 = j1+j2+j3;
        HalfInteger j1J2J3 = j1+J2+J3;
        HalfInteger J1j2J3 = J1+j2+J3;
        HalfInteger J1J2j3 = J1+J2+j3;

        if (j1j2j3.get_twice() % 2 != 0 ||
            j1J2J3.get_twice() % 2 != 0 ||
            J1j2J3.get_twice() % 2 != 0 ||
            J1J2j3.get_twice() % 2 != 0) {
                return 0.0;
        }

        HalfInteger j1j2J1J2 = j1+j2+J1+J2;
        HalfInteger j2j3J2J3 = j2+j3+J2+J3;
        HalfInteger j3j1J3J1 = j3+j1+J3+J1;

        int min_t = ceiling(max(max(j1j2j3, j1J2J3), max(J1j2J3, J1J2j3)));
        int max_t = floor(min(min(j1j2J1J2, j2j3J2J3), j3j1J3J1));
        if (min_t <= max_t) {
            double sum = 0.0;
            for (int ti = min_t; ti <= max_t; ++ti) {
                double term = ((ti % 2 == 1) ? -1 : 1)
                            * factorial(HalfInteger(ti + 1))
                            / f_term(HalfInteger(ti), j1j2j3, j1J2J3, J1j2J3, J1J2j3,
                                     j1j2J1J2, j2j3J2J3, j3j1J3J1);
                sum += term;
            }
            double mult = std::sqrt(
                Delta(j1, j2, j3)
              * Delta(j1, J2, J3)
              * Delta(J1, j2, J3)
              * Delta(J1, J2, j3)
            );
            value = mult * sum;
        }
    }
    return value;
}


double wigner9j(HalfInteger j1, HalfInteger j2, HalfInteger j3,
                HalfInteger j4, HalfInteger j5, HalfInteger j6,
                HalfInteger j7, HalfInteger j8, HalfInteger j9) {
    double value = 0.0;
    int kmin = ceiling(max(max(abs(j1 - j9), abs(j4 - j8)), abs(j2 - j6)));
    int kmax = floor(min(min(abs(j1 + j9), abs(j4 + j8)), abs(j2 + j6)));
    if (kmax >= kmin) {
        double sum = 0.0;
        for (int ki = kmin; ki <= kmax; ++ki) {
            HalfInteger k(ki);
            sum += (2 * ki + 1)
                 * wigner6j(j1, j4, j7, j8, j9, k)
                 * wigner6j(j2, j5, j8, j4, k, j6)
                 * wigner6j(j3, j6, j9, k, j1, j2);
        }
        value = sum;
    }
    return value;
}
} // namespace wigner_symbol
