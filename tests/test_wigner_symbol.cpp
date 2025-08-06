#include <iostream>
#include "wigner_symbol.hpp"
#include "half_integer.hpp"

int main() {
    using namespace half_integer;
    using namespace wigner_symbol;
    HalfInteger zero = HalfInteger(0);
    HalfInteger one = HalfInteger(1);
    HalfInteger two = HalfInteger(2);
    std::cout << wigner3j(one, one, zero, zero, zero, zero) << std::endl;
    std::cout << wigner6j(one, one, one, one, one, one) << std::endl;
    std::cout << wigner9j(one, one, one, two, two, two, one, one, one) << std::endl;
    return 0;
}
