#ifndef HALF_INTEGER_HPP
#define HALF_INTEGER_HPP

#include <stdexcept>
#include <iostream>

namespace half_integer {
class HalfInteger {
private:
    int twice_value;

public:
    explicit HalfInteger(int numerator = 0) : twice_value(numerator * 2) {}
    static HalfInteger from_twice(int t) { return HalfInteger(t / 2); }

    int get_twice() const { return twice_value; }
    explicit operator double() const { return static_cast<double>(twice_value) / 2.0; }
    explicit operator int() const { return twice_value / 2; }

    HalfInteger operator+(const HalfInteger& other) const { return from_twice(twice_value + other.twice_value); }
    HalfInteger operator-(const HalfInteger& other) const { return from_twice(twice_value - other.twice_value); }
    HalfInteger operator-() const {return from_twice(-twice_value);}
    HalfInteger operator++() const {return from_twice(twice_value + 1);} // actually only adds 0.5.
    bool operator<(const HalfInteger& other) const { return twice_value < other.twice_value; }
    bool operator<=(const HalfInteger& other) const { return twice_value <= other.twice_value; }
    bool operator>(const HalfInteger& other) const { return twice_value > other.twice_value; }
    bool operator>=(const HalfInteger& other) const { return twice_value >= other.twice_value; }
    bool operator==(const HalfInteger& other) const { return twice_value == other.twice_value; }

    friend HalfInteger abs(const HalfInteger& a) {
        return from_twice(std::abs(a.get_twice()));
    }

    friend HalfInteger max(const HalfInteger& a, const HalfInteger& b) {
        return (a.twice_value >= b.twice_value) ? a : b;
    }

    friend HalfInteger min(const HalfInteger& a, const HalfInteger& b) {
        return (a.twice_value <= b.twice_value) ? a : b;
    }

    friend std::ostream& operator<<(std::ostream& os, const HalfInteger& h) {
        os << static_cast<double>(h);
        return os;
    }
};
inline double factorial(HalfInteger x) {
    double v = static_cast<double>(x);
    return std::tgamma(v + 1.0);
}
inline int ceiling(HalfInteger x) {
    return (x.get_twice() + 1) / 2;
}
inline int floor(HalfInteger x) {
    return (x.get_twice()) / 2;
}
} // namespace half_integer

#endif
