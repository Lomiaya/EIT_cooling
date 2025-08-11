// hunds_case_b.cpp
#include "hunds_case_b.hpp"
#include <cmath>
#include <iostream>

namespace hunds_case_b {
using namespace half_integer;
using namespace wigner_symbol;
// spherical tensor data
extern const double T_Case_B[3][3] = {
    {  0.0,                    0.0,                    0.0                   },
    { -2.0/std::sqrt(3.0),     0.0, -2.0/std::sqrt(6.0) },
    {  0.0,                    0.0,                    0.0                   }
};

void unpack(const HundsCaseB_Rot& s,
            int& v, HalfInteger& S, HalfInteger& I,
            HalfInteger& L, HalfInteger& N,
            HalfInteger& J, HalfInteger& F, HalfInteger& M)
{
    v = s.v;  S = s.S;  I = s.I;  L = s.Lambda;
    N = s.N;  J = s.J;  F = s.F;  M = s.M;
}

bool delta(HalfInteger a, HalfInteger b) { return a == b; }

// Rotation, not tested
double Rotation(const HundsCaseB_Rot& a, const HundsCaseB_Rot& b) {
    int va, vb; HalfInteger Sa,Ia,La,Na,Ja,Fa,Ma;
    HalfInteger Sb,Ib,Lb,Nb,Jb,Fb,Mb;
    unpack(a, va, Sa, Ia, La, Na, Ja, Fa, Ma);
    unpack(b, vb, Sb, Ib, Lb, Nb, Jb, Fb, Mb);
    if (!delta(La,Lb)||!delta(Na,Nb)||!delta(Ja,Jb)||!delta(Fa,Fb)||!delta(Ma,Mb)||va!=vb)
        return 0.0;
    double Na_double = double(Na);
    double La_double = double(La);
    return Na_double*(Na_double+1) - La_double*La_double;
}

// RotationDistortion, not tested
double RotationDistortion(const HundsCaseB_Rot& a, const HundsCaseB_Rot& b) {
    double r = Rotation(a,b);
    return -r * r;
}

// SpinRotation
double SpinRotation(const HundsCaseB_Rot& a, const HundsCaseB_Rot& b) {
    int va, vb; HalfInteger Sa,Ia,La,Na,Ja,Fa,Ma;
    HalfInteger Sb,Ib,Lb,Nb,Jb,Fb,Mb;
    unpack(a, va, Sa, Ia, La, Na, Ja, Fa, Ma);
    unpack(b, vb, Sb, Ib, Lb, Nb, Jb, Fb, Mb);
    double Sad = double(Sa);
    double Iad = double(Ia);
    double Lad = double(La);
    double Nad = double(Na);
    double Jad = double(Ja);
    double Fad = double(Fa);
    double Mad = double(Ma);
    double Sbd = double(Sb);
    double Ibd = double(Ib);
    double Lbd = double(Lb);
    double Nbd = double(Nb);
    double Jbd = double(Jb);
    double Fbd = double(Fb);
    double Mbd = double(Mb);
    if (!delta(Ja,Jb)||!delta(Fa,Fb)||!delta(Ma,Mb)||va!=vb)
        return 0.0;
    int phase_ = (Jb + Sa + Na).get_twice();
    if (phase_ % 2 != 0) {
        std::cout << "Invalid (Jb, Sa, Na) in SpinRotation: " << Jb << " " << Sa << " " << Na << " " << std::endl;
        return 0.0; // should never happen
    }
    phase_ /= 2;
    double phase = (phase_ % 2 == 0? 1.0: -1.0);
    int phase2_ = (Na - La).get_twice();
    if (phase2_ % 2 != 0) {
        std::cout << "Invalid (Na, La) in SpinRotation: " << Na << " " << La << std::endl;
        return 0.0; // should never happen
    }
    phase2_ /= 2;
    double phase2 = (phase2_ % 2 == 0? 1.0: -1.0);
    double pref = 0.5 * sqrt(Sad * (Sad + 1) * (2 * Sad + 1) * (2 * Nad + 1) * (2 * Nbd + 1));
    double wg6j = wigner6j(Nb, Sa, Ja, Sa, Na, HalfInteger(1));
    double sum = 0.0;
    for (int k = 0; k <= 2; ++k) for (int q = -1; q <= 1; ++q) {
        double mult_phase = (k % 2 == 1)? -1.0: 1.0;
        double mult_pref1 = sqrt(Nbd * (Nbd + 1) * (2 * Nbd + 1));
        double mult_pref2 = sqrt(Nad * (Nad + 1) * (2 * Nad + 1));
        sum += sqrt(2 * k + 1) *
        (
            wigner6j(HalfInteger(1), HalfInteger(1), HalfInteger(k), Na, Nb, Nb) * mult_pref1 * mult_phase +
            wigner6j(HalfInteger(1), HalfInteger(1), HalfInteger(k), Nb, Na, Na) * mult_pref2
        ) * 
        wigner3j(Na, HalfInteger(k), Nb, -La, HalfInteger(q), Lb) * T_Case_B[q+1][k];
    }
    return phase * phase2 * pref * wg6j * sum;
}

// Hyperfine_IS, check
double Hyperfine_IS(const HundsCaseB_Rot& a, const HundsCaseB_Rot& b) {
    int va, vb; HalfInteger Sa,Ia,La,Na,Ja,Fa,Ma;
    HalfInteger Sb,Ib,Lb,Nb,Jb,Fb,Mb;
    unpack(a, va, Sa, Ia, La, Na, Ja, Fa, Ma);
    unpack(b, vb, Sb, Ib, Lb, Nb, Jb, Fb, Mb);
    double Sad = double(Sa);
    double Iad = double(Ia);
    double Lad = double(La);
    double Nad = double(Na);
    double Jad = double(Ja);
    double Fad = double(Fa);
    double Mad = double(Ma);
    double Sbd = double(Sb);
    double Ibd = double(Ib);
    double Lbd = double(Lb);
    double Nbd = double(Nb);
    double Jbd = double(Jb);
    double Fbd = double(Fb);
    double Mbd = double(Mb);
    if (!delta(La,Lb)||!delta(Na,Nb)||!delta(Fa,Fb)||!delta(Ma,Mb)||va!=vb)
        return 0.0;
    int phase_ = (Nb + Sa + Ja).get_twice();
    if (phase_ % 2 != 0) {
        std::cout << "Invalid (Nb, Sa, Ja) in Hyperfine_IS: " << Nb << " " << Sa << " " << Ja << std::endl;
        return 0.0; // should never happen
    }
    phase_ /= 2;
    double phase = (phase_ % 2 == 0? 1.0: -1.0);
    int phase2_ = (Jb + Ia + Fb + HalfInteger(1)).get_twice();
    if (phase2_ % 2 != 0) {
        std::cout << "Invalid (Jb, Ia, Fb) in Hyperfine_IS: " << Jb << " " << Ia << " " << Fb << std::endl;
        return 0.0; // should never happen
    }
    phase2_ /= 2;
    double phase2 = (phase2_ % 2 == 0? 1.0: -1.0);
    double pref = sqrt((2 * Jbd + 1) * (2 * Jad + 1) * 
                       Sad * (Sad + 1) * (2 * Sad + 1) *
                       Iad * (Iad + 1) * (2 * Iad + 1));
    return phase * phase2 * pref *
           wigner6j(Ia, Ja, Fb, Jb, Ia, HalfInteger(1)) *
           wigner6j(Sa, Ja, Nb, Jb, Sa, HalfInteger(1));
}

// Hyperfine_Dipolar
double Hyperfine_Dipolar(const HundsCaseB_Rot& a, const HundsCaseB_Rot& b) {
    int va, vb; HalfInteger Sa,Ia,La,Na,Ja,Fa,Ma;
    HalfInteger Sb,Ib,Lb,Nb,Jb,Fb,Mb;
    unpack(a, va, Sa, Ia, La, Na, Ja, Fa, Ma);
    unpack(b, vb, Sb, Ib, Lb, Nb, Jb, Fb, Mb);
    double Sad = double(Sa);
    double Iad = double(Ia);
    double Lad = double(La);
    double Nad = double(Na);
    double Jad = double(Ja);
    double Fad = double(Fa);
    double Mad = double(Ma);
    double Sbd = double(Sb);
    double Ibd = double(Ib);
    double Lbd = double(Lb);
    double Nbd = double(Nb);
    double Jbd = double(Jb);
    double Fbd = double(Fb);
    double Mbd = double(Mb);
    if (!delta(Fa,Fb)||!delta(Ma,Mb)||va!=vb)
        return 0.0;
    int phase_ = (Na - La).get_twice();
    if (phase_ % 2 != 0) {
        std::cout << "Invalid (Na, La) in Hyperfine_IS: " << Na << " " << La << std::endl;
        return 0.0; // should never happen
    }
    phase_ /= 2;
    double phase = (phase_ % 2 == 0? 1.0: -1.0);
    int phase2_ = (Jb + Ia + Fa + HalfInteger(1)).get_twice();
    if (phase2_ % 2 != 0) {
        std::cout << "Invalid (Jb, Ia, Fa) in Hyperfine_IS: " << Jb << " " << Ia << " " << Fa << std::endl;
        return 0.0; // should never happen
    }
    phase2_ /= 2;
    double phase2 = (phase2_ % 2 == 0? 1.0: -1.0);
    double pref = sqrt((2 * Jbd + 1) * (2 * Jad + 1) * 
                       (2 * Nbd + 1) * (2 * Nad + 1) *
                       Sad * (Sad + 1) * (2 * Sad + 1) *
                       Iad * (Iad + 1) * (2 * Iad + 1) *
                       30);
    return phase * phase2 * pref *
           wigner6j(Ia, Ja, Fb, Jb, Ia, HalfInteger(1)) *
           wigner9j(Na, Nb, HalfInteger(2), Sa, Sa, HalfInteger(1), Ja, Jb, HalfInteger(1)) *
           wigner3j(Na, HalfInteger(2), Nb, -La, HalfInteger(0), Lb);
}

// TDM
double TDM(const HundsCaseB_Rot& a, const HundsCaseB_Rot& b, int p) {
    int va, vb; HalfInteger Sa,Ia,La,Na,Ja,Fa,Ma;
    HalfInteger Sb,Ib,Lb,Nb,Jb,Fb,Mb;
    unpack(a, va, Sa, Ia, La, Na, Ja, Fa, Ma);
    unpack(b, vb, Sb, Ib, Lb, Nb, Jb, Fb, Mb);
    double Sad = double(Sa);
    double Iad = double(Ia);
    double Lad = double(La);
    double Nad = double(Na);
    double Jad = double(Ja);
    double Fad = double(Fa);
    double Mad = double(Ma);
    double Sbd = double(Sb);
    double Ibd = double(Ib);
    double Lbd = double(Lb);
    double Nbd = double(Nb);
    double Jbd = double(Jb);
    double Fbd = double(Fb);
    double Mbd = double(Mb);
    if (!delta(Ia,Ib)||!delta(Sa,Sb))
        return 0.0;
    HalfInteger q_ = La - Lb;
    int doubleq = q_.get_twice();
    if (doubleq != -2 && doubleq != 0 && doubleq != 2) {
        std::cout << "Invalid La, Lb in TDM, they are: " << La << " " << Lb << std::endl;
        return 0.0;
    }
    int q = doubleq / 2;

    if (!(Ma + HalfInteger(p) == Mb))
        return 0.0;
    if (!(abs(Fa - HalfInteger(1)) <= Fb && Fb <= Fa + HalfInteger(1)))
        return 0.0;
    if (!(abs(Na - HalfInteger(1)) <= Nb && Nb <= Na + HalfInteger(1)))
        return 0.0;

    // m1
    int phase1_ = (Fa - Ma).get_twice();
    if (phase1_ % 2 != 0) {
        std::cout << "Invalid (Fa, Ma) in Hyperfine_IS: " << Fa << " " << Ma << std::endl;
        return 0.0; // should never happen
    }
    phase1_ /= 2;
    double phase1 = (phase1_ % 2 == 0? 1.0: -1.0);
    double m1 = phase1 *
                wigner3j(Fa, HalfInteger(1), Fb, -Ma, HalfInteger(-p), Mb);
    // m2
    int phase2_ = (Fb + Ja + Ia + HalfInteger(1)).get_twice();
    if (phase2_ % 2 != 0) {
        std::cout << "Invalid (Fb, Ja, Ia) in Hyperfine_IS: " << Fb << " " << Ja << " " << Ia << std::endl;
        return 0.0; // should never happen
    }
    phase2_ /= 2;
    double phase2 = (phase2_ % 2 == 0? 1.0: -1.0);
    double m2 = phase2 *
                sqrt((2 * Fad + 1) * (2 * Fbd + 1)) *
                wigner6j(Jb, Fb, Ia, Fa, Ja, HalfInteger(1));
    // m3
    int phase3_ = (Jb + Na + Sa + HalfInteger(1)).get_twice();
    if (phase3_ % 2 != 0) {
        std::cout << "Invalid (Jb, Na, Sa) in Hyperfine_IS: " << Jb << " " << Na << " " << Sa << std::endl;
        return 0.0; // should never happen
    }
    phase3_ /= 2;
    double phase3 = (phase3_ % 2 == 0? 1.0: -1.0);
    double m3 = phase3 *
                sqrt((2 * Jad + 1) * (2 * Jbd + 1)) *
                wigner6j(Nb, Jb, Sa, Ja, Na, HalfInteger(1));
    // m4
    int phase4_ = (Na - La).get_twice();
    if (phase4_ % 2 != 0) {
        std::cout << "Invalid (Na, La) in Hyperfine_IS: " << Na << " " << La << std::endl;
        return 0.0; // should never happen
    }
    phase4_ /= 2;
    double phase4 = (phase4_ % 2 == 0? 1.0: -1.0);
    double m4 = phase4 *
                sqrt((2 * Nad + 1) * (2 * Nbd + 1)) *
                wigner3j(Na, HalfInteger(1), Nb, -La, HalfInteger(q), Lb);

    return m1 * m2 * m3 * m4;
}

// Zeeman, only for sigma states!
double Zeeman(const HundsCaseB_Rot &a, const HundsCaseB_Rot &b, int p) {
    int va, vb; HalfInteger Sa,Ia,La,Na,Ja,Fa,Ma;
    HalfInteger Sb,Ib,Lb,Nb,Jb,Fb,Mb;
    unpack(a, va, Sa, Ia, La, Na, Ja, Fa, Ma);
    unpack(b, vb, Sb, Ib, Lb, Nb, Jb, Fb, Mb);
    double Sad = double(Sa);
    double Iad = double(Ia);
    double Lad = double(La);
    double Nad = double(Na);
    double Jad = double(Ja);
    double Fad = double(Fa);
    double Mad = double(Ma);
    double Sbd = double(Sb);
    double Ibd = double(Ib);
    double Lbd = double(Lb);
    double Nbd = double(Nb);
    double Jbd = double(Jb);
    double Fbd = double(Fb);
    double Mbd = double(Mb);
    if (!delta(Na,Nb)||va!=vb)
        return 0.0;
    if (!(La == HalfInteger(0) && Lb == HalfInteger(0))) {
        std::cout << "Invalid La or Lb, should be both zero in Zeeman: " << La << " " << Lb << std::endl;
        return 0.0;
    }
    
    int phase1_ = (Fb - Mb).get_twice();
    if (phase1_ % 2 != 0) {
        std::cout << "Invalid (Fb, Mb) in Hyperfine_IS: " << Fb << " " << Mb << std::endl;
        return 0.0; // should never happen
    }
    phase1_ /= 2;
    double phase1 = (phase1_ % 2 == 0? 1.0: -1.0);
    int phase2_ = (Jb + Ia + Fa + HalfInteger(1)).get_twice();
    if (phase2_ % 2 != 0) {
        std::cout << "Invalid (Jb, Ia, Fa) in Hyperfine_IS: " << Jb << " " << Ia << " " << Fa << std::endl;
        return 0.0; // should never happen
    }
    phase2_ /= 2;
    double phase2 = (phase2_ % 2 == 0? 1.0: -1.0);
    int phase3_ = (Na + Sa + Jb + HalfInteger(1)).get_twice();
    if (phase3_ % 2 != 0) {
        std::cout << "Invalid (Na, Sa, Jb) in Hyperfine_IS: " << Na << " " << Sa << " " << Jb << std::endl;
        return 0.0; // should never happen
    }
    phase3_ /= 2;
    double phase3 = (phase3_ % 2 == 0? 1.0: -1.0);

    return (p % 2 == 0? 1.0: -1.0) * phase1 * phase2 * phase3 *
           sqrt(
                (2 * Fad + 1) * (2 * Fbd + 1) *
                (2 * Jad + 1) * (2 * Jbd + 1) *
                Sad * (Sad + 1) * (2 * Sad + 1)
            ) *
            wigner3j(Fb, HalfInteger(1), Fa, -Mb, HalfInteger(-p), Ma) *
            wigner6j(Jb, Fb, Ia, Fa, Ja, HalfInteger(1)) *
            wigner6j(Sa, Jb, Na, Ja, Sa, HalfInteger(1));
}
double Zeeman_z(const HundsCaseB_Rot &a, const HundsCaseB_Rot &b) {
    return Zeeman(a, b, 0);
}
} // namespace hunds_case_b
