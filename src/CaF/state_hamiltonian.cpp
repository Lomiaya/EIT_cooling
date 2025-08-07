// state_hamiltonian.cpp
#include "state_hamiltonian.hpp"
#include <algorithm>

namespace state_hamiltonian {
    std::vector<double> get_eigenvalues(MatrixXc matrix) {
        auto eigenvalues = matrix.eigenvalues();
        std::vector<double> casted_eigenvalues;
        for (auto& eigenvalue: eigenvalues) {
            casted_eigenvalues.push_back(eigenvalue.real());
        }
        std::sort(casted_eigenvalues.begin(), casted_eigenvalues.end());
        return casted_eigenvalues;
    }
    void StateHamiltonian::add_operator_to_matrix(
        double weight,
        double (*op)(const hunds_case_b::HundsCaseB_Rot&, const hunds_case_b::HundsCaseB_Rot&)
    ) {
        int num_states = basis_states.size();
        for (int i = 0; i < num_states; ++i) for (int j = 0; j < num_states; ++j) {
            matrix(i, j) += op(basis_states[i], basis_states[j]);
        }
    }
    
    void StateHamiltonian::add_operators_to_matrix(
        std::vector<double> weights, 
        std::vector<double (*)(const hunds_case_b::HundsCaseB_Rot&, const hunds_case_b::HundsCaseB_Rot&)> ops
    ) {
        int num_states = basis_states.size();
        int num_ops = ops.size();
        // ops should be large in memory, want to loop over it the last.
        for (int k = 0; k < num_states; ++k) for (int i = 0; i < num_states; ++i) for (int j = 0; j < num_states; ++j) {
            matrix(i, j) += ops[k](basis_states[i], basis_states[j]) * weights[k];
        }
    }

    std::vector<Eigen::MatrixXd> calculate_TDM(
        const StateHamiltonian& X,
        const StateHamiltonian& A,
        double (*op)(const hunds_case_b::HundsCaseB_Rot&, const hunds_case_b::HundsCaseB_Rot&, int)
    ) { // shape: [3][n_excited_states][n_ground_states]
        int n_excited_states = A.basis_states.size();
        int n_ground_states = X.basis_states.size();
        Eigen::MatrixXd G1m = Eigen::MatrixXd(n_excited_states, n_ground_states);
        Eigen::MatrixXd G0 = Eigen::MatrixXd(n_excited_states, n_ground_states);
        Eigen::MatrixXd G1 = Eigen::MatrixXd(n_excited_states, n_ground_states);
        std::vector<Eigen::MatrixXd> G = {G1m, G0, G1};
        for (int p = -1; p <= 1; ++p) for (int i = 0; i < n_excited_states; ++i) for (int j = 0; j < n_ground_states; ++j) {
            auto TDM = op(A.basis_states[i], X.basis_states[j], p);
            G[p+1](i, j) = std::abs(TDM) * std::abs(TDM);
        }
        return G;
    }
} // namespace state_hamiltonian
