// hunds_case_b.hpp
#ifndef HUNDS_CASE_B_HPP
#define HUNDS_CASE_B_HPP

#include "half_integer.hpp"
#include "wigner_symbol.hpp"
#include <array>
#include <complex>
#include <string>

namespace hunds_case_b {
using namespace half_integer;
using namespace wigner_symbol;
// spherical tensor T_Case_B (3×3)
extern const double T_Case_B[3][3];

// abstract base
struct HundsCaseB {
    virtual ~HundsCaseB() = default;
};

// rotational basis struct
struct HundsCaseB_Rot : HundsCaseB {
    std::string  label;
    int          v;
    HalfInteger  S, I, Lambda, N, J, F, M;

    HundsCaseB_Rot(
        std::string label_,
        int v_,
        HalfInteger S_,
        HalfInteger I_,
        HalfInteger Lambda_,
        HalfInteger N_,
        HalfInteger J_,
        HalfInteger F_,
        HalfInteger M_
    ) : label(std::move(label_)), v(v_), S(S_), I(I_), Lambda(Lambda_),
        N(N_), J(J_), F(F_), M(M_) {}
};

// unpack helper
void unpack(const HundsCaseB_Rot& s,
            int& v, HalfInteger& S, HalfInteger& I,
            HalfInteger& L, HalfInteger& N,
            HalfInteger& J, HalfInteger& F, HalfInteger& M);

// delta equality
bool δ(HalfInteger a, HalfInteger b);

// matrix‐elements
double Rotation(const HundsCaseB_Rot& a, const HundsCaseB_Rot& b);
double RotationDistortion(const HundsCaseB_Rot& a, const HundsCaseB_Rot& b);
double SpinRotation(const HundsCaseB_Rot& a, const HundsCaseB_Rot& b);
double Hyperfine_IS(const HundsCaseB_Rot& a, const HundsCaseB_Rot& b);
double Hyperfine_Dipolar(const HundsCaseB_Rot& a, const HundsCaseB_Rot& b);

// Transition Dipole Moment
double TDM(const HundsCaseB_Rot& a, const HundsCaseB_Rot& b, int p);
double TDM(const HundsCaseB_Rot& a, const HundsCaseB_Rot& b);

// Zeeman interactions
double Zeeman(const HundsCaseB_Rot &a, const HundsCaseB_Rot &b, int p);
double Zeeman_z(const HundsCaseB_Rot &a, const HundsCaseB_Rot &b);
} // namespace hunds_case_b
#endif // HUNDS_CASE_B_HPP
