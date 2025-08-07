#include <iostream>
#include "wigner_symbol.hpp"
#include "half_integer.hpp"

int main() {
    using namespace half_integer;
    using namespace wigner_symbol;
    HalfInteger zero = HalfInteger(0);
    HalfInteger one = HalfInteger(1);
    HalfInteger two = HalfInteger(2);
    HalfInteger three = HalfInteger(3);
    HalfInteger minusone = HalfInteger(-1);
    HalfInteger half = HalfInteger::from_twice(1);
    std::cout << wigner3j(half, half, one, half, half, minusone) << std::endl;
    std::cout << wigner6j(one, half, half, half, one, one) << std::endl;
    std::cout << wigner9j(one, one, one, three, three, three, two, two, two) << std::endl;
    return 0;
}
