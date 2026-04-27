#include <cmath>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <random>
#include <numeric>
#include <tuple>

class NewtonsMethod
{
public:
    NewtonsMethod(std::vector<std::vector<long double>> A_in, int maxiter = 1000, long double tol = 1e-15, bool verbose = true) {
        // Check: Matrix A empty
        if (A_in.empty() || A_in[0].empty()) {
            throw std::invalid_argument("Matrix A cannot be empty.");
        }

        size_t A_rows = A_in.size();
        size_t A_cols = A_in[0].size();

        // Check: Matrix A dimension
        if (A_rows != A_cols) {
            throw std::invalid_argument("Matrix A must be square (n x n).");
        }

        // Check: Matrix A symmetry
        for (size_t i = 0; i < A_rows; i++) {
            for (size_t j = i + 1; j < A_cols; j++) {
                if (std::abs(A_in[i][j] - A_in[j][i]) > 1e-15) { 
                    throw std::invalid_argument("Matrix A must be symmetric.");
                }
            }
        }

        this->A = A_in;                     // Matrix, expects symmetric (n x n)               
        this->maxiter = maxiter;            // Maximum number of iterations
        this->tol = tol;                    // Error tolerance
        this->verbose = verbose;            // If true outputs convergence, error, eigenvalue, and eigenvector information    
    }

    std::pair<std::vector<long double>, long double> solver() {
        int A_n = A.size();
        int J_size = A_n + 1;

        bool converged = false;

        std::vector<long double> E;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<long double> dis(1.0, 20.0);

        // Initialize a 'guess' for vector V, a vector with random values
        std::vector<long double> V(A_n);
        for(long double &val : V) val = dis(gen);

        // Initialize a 'guess' for scalar lambda, a random valued scalar
        long double lambda = dis(gen);

        // Newton's Method for Systems
        for (int k = 0; k < maxiter; k++) {
            
            // 1: Compute F vector and Jacobian matrix
            std::vector<long double> F(A_n+1, 0.0); 
            long double sum_sq_v = 0.0;

            std::vector<std::vector<long double>> J(J_size, std::vector<long double>(J_size, 0.0));

            for (int i = 0; i < A_n; i++) {
                long double F_row_sum = 0.0;

                long double V_i = V[i];
                J[i][A_n] = -V_i;
                J[A_n][i] = V_i;
                
                for (int j = 0; j < A_n; j++) {
                    long double A_ij = A[i][j];
                    F_row_sum += A_ij * V[j];
                    
                    long double J_ij;
                    if (i == j) {
                        J_ij = A_ij - lambda;
                    } else {
                        J_ij = A_ij;
                    }
                    J[i][j] = J_ij;
                }
            
                F[i] = F_row_sum - (lambda * V_i);
                sum_sq_v += V_i * V_i; 
            }

            J[A_n][A_n] = 0;
            F[A_n] = (0.5 * sum_sq_v) - 1.0;

            // 2: Convergence check (using infinity norm)
            long double max_F = 0.0;
            for (const auto& val : F) {
                max_F = std::max(max_F, std::abs(val));
            }

            E.push_back(max_F);

            if (max_F < tol) {
                converged = true;
                if (verbose == true) {
                    std::cout << "Converged in " << k << " iterations" << std::endl;
                }
                break;
            }


            // 3: REF using Gaussian elimination with complete pivoting
            std::vector<int> C(J_size);                              // Column index vector
            std::iota(C.begin(), C.end(), 0);                       

            for (int l = 0; l < J_size - 1; l++) {
                // Complete pivoting strategy
                long double max_el = std::abs(J[l][l]);
                int max_el_i = l;
                int max_el_j = l;

                for (int i = l; i < J_size; i++) {
                    for (int j = l; j < J_size; j++){

                        long double abs_J_ij = std::abs(J[i][j]);
                        if (abs_J_ij > max_el) {                    // Find largest element in absolute value
                            max_el = abs_J_ij;
                            max_el_i = i;
                            max_el_j = j;
                        }
                    }
                }
                std::swap(C[l], C[max_el_j]);                        // Keep track of column indices

                std::swap(J[l], J[max_el_i]);                        // Swap rows in J
                std::swap(F[l], F[max_el_i]);                        // Swap rows in F
                for (int row = 0; row < J_size; row++) {             // Swap columns in J
                    std::swap(J[row][l], J[row][max_el_j]);
                }

                for (int i = l + 1; i < J_size; i++) {
                    long double J_il = J[i][l];
                    if (J_il == 0.0){continue;}

                    long double J_ll = J[l][l];
                    if (J_ll == 0.0){
                        throw std::invalid_argument("Matrix A must be invertible and symmetric.");
                    }

                    long double scaling_factor = J_il / J_ll;         // Compute scaling factor for row elimination
                    J[i][l] = 0.0;                                    // Set values below pivot to zero
                    F[i] -= scaling_factor * F[l];                    // Row elimination for F vector

                    for (int j = l + 1; j < J_size; j++){
                        J[i][j] -= scaling_factor * J[l][j];          // Row elimination for J matrix
                    }
                }

            }

            // 4: Backsolve J*D = F for D
            std::vector<long double> D(J_size, 0.0);
            for (int i = J_size - 1; i >= 0; i--) {
                D[i] = F[i];
                for (int j = J_size - 1; j > i; j--) {
                    D[i] -= J[i][j] * D[j];
                }
                D[i] /= J[i][i];
            }

            // 5. Correct D for column permutations
            std::vector<long double> D_fixed(J_size);
            for (int i = 0; i < J_size; i++) {
                D_fixed[C[i]] = D[i];
            }

            // 6: Compute next V and lambda values
            for (int i = 0; i < A_n; i++) {
                V[i] -= D_fixed[i]; 
            }
            lambda -= D_fixed[A_n];

        }

        if (verbose == true && converged == true) {
            std::cout << "Eigenvalue: " << lambda << std::endl;
            
            std::cout << "Eigenvector: (";
            for (size_t i = 0; i < V.size(); ++i) {
                std::cout << V[i];
                if (i < V.size() - 1) {
                    std::cout << ", ";
                }
            }
            std::cout << ")" << std::endl;
            
            std::cout << "Error Sequence: ";
            for (size_t i = 0; i < E.size(); ++i) {
                std::cout << E[i];
                if (i < E.size() - 1) {
                    std::cout << ", ";
                }
            }
    
        }

        if (converged == true) {
            return {V, lambda};
        } else {
            std::cout << "Newton's method failed to converge in " << maxiter << " iterations" << std::endl;
            return {{}, std::numeric_limits<long double>::quiet_NaN()};
        }
    }

private:
    std::vector<std::vector<long double>> A;
    int maxiter;
    long double tol;
    bool verbose;
};