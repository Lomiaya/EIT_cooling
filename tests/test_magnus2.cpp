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
    double freq = 0.0;
    std::vector<std::pair<MatrixXc,double>> Hs = {{H0,0},{H1,freq},{H1b,-freq}};
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i=0; i<20; ++i) {
        auto U_bf = brute_force(Hs,G,t0,t1);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "Time taken by function: " << duration << " milliseconds" << std::endl;
    auto start_time2 = std::chrono::high_resolution_clock::now();
    for (int i=0; i<1000; ++i) {
        auto U_m2 = magnus2_analytic(Hs,G,t0,t1);
    }
    auto end_time2 = std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(end_time2 - start_time2).count();
    std::cout << "Time taken by function: " << duration2 << " milliseconds" << std::endl;
    auto U_bf = brute_force(Hs,G,t0,t1);
    auto U_m2 = magnus2_analytic(Hs,G,t0,t1);
    std::cout << "U_bf=\n"<<U_bf<<"\n";
    std::cout << "U_m2=\n"<<U_m2<<"\n";
    std::cout << "U_m2-U_bf=\n"<<(U_m2-U_bf)<<"\n";
    return 0;
}
