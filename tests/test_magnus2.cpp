// Test file: test_magnus2.cpp
#include "stochastic_schrodinger_spin_motion.hpp"
#include <iostream>

int main() {
    using namespace ss_spin;
    MatrixXc H0(2,2), H1(2,2), H1b(2,2), G(2,2);
    H0 << 0,0, 0,1;
    H1 << 0,0, 1,0;
    H1b<< 0,1, 0,0;
    G  << 1,0, 0,0;
    double t0=0, t1=0.5;
    std::vector<std::pair<MatrixXc,double>> Hs = {{H0,0},{H1,1},{H1b,-1}};
    auto U_bf = brute_force(Hs,G,t0,t1);
    auto U_m2 = magnus2_analytic(Hs,G,t0,t1);
    std::cout << "U_bf=\n"<<U_bf<<"\n";
    std::cout << "U_m2=\n"<<U_m2<<"\n";
    std::cout << "U_m2-U_bf=\n"<<(U_m2-U_bf)<<"\n";
    return 0;
}
