// state_hamiltonian.hpp
#ifndef STATE_HAMILTONIAN_HPP
#define STATE_HAMILTONIAN_HPP

#include <Eigen/Dense>
#include <complex>
#include <vector>
#include "hunds_case_b.hpp"

namespace state_hamiltonian {

using Complex = std::complex<double>;
using MatrixXc = Eigen::MatrixXcd;
using VectorXc = Eigen::VectorXcd;
using VecD = Eigen::VectorXd;

class StateHamiltonian {
    private:
        MatrixXc matrix;
        std::vector<hunds_case_b::HundsCaseB_Rot> basis_states;
    public:
        explicit StateHamiltonian(std::vector<hunds_case_b::HundsCaseB_Rot> states_init) : basis_states(states_init), matrix(MatrixXc(states_init.size(), states_init.size())) {
            matrix.setZero();
        }
        MatrixXc get_matrix() const { return matrix; }
        std::vector<hunds_case_b::HundsCaseB_Rot> get_states() const { return basis_states; }
        void add_operator_to_matrix(
            double weight,
            double (*op)(const hunds_case_b::HundsCaseB_Rot&, const hunds_case_b::HundsCaseB_Rot&)
        );
        void add_operators_to_matrix(
            std::vector<double> weights, 
            std::vector<double (*)(const hunds_case_b::HundsCaseB_Rot&, const hunds_case_b::HundsCaseB_Rot&)> ops
        );
        friend std::vector<Eigen::MatrixXd> calculate_TDM(
            const StateHamiltonian& X,
            const StateHamiltonian& A,
            double (*op)(const hunds_case_b::HundsCaseB_Rot&, const hunds_case_b::HundsCaseB_Rot&, int)
        ); // shape: [3][n_excited_states][n_ground_states]
};

std::vector<double> get_eigenvalues(MatrixXc matrix);

} // namespace state_hamiltonian

#endif // STATE_HAMILTONIAN_HPP