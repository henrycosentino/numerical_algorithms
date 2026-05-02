#include <cmath>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <random>


class PowerMethod
{
public: 
    PowerMethod(std::vector<std::vector<long double>> A_in, int maxiter = 1000, long double tol = 1e-15, bool verbose = true) {
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
        this->verbose = verbose;            // If true outputs convergence, error, and eigenvalue information    
    }

    long double _l2_norm(std::vector<long double> X) {
        // Calculates the L2-Norm of a vector, returning a scalar
        int n = X.size();
        long double norm_factor = 0.0L;

        for (int i = 0; i < n; i++) {
            norm_factor += X[i]*X[i];
        }
        norm_factor = std::pow(norm_factor, 0.5);

        return norm_factor;
    }

    std::vector<long double> _unit_vec(std::vector<long double> X) {
        // Normalizes a vector to unit magnitude using the L2-Norm, returning a vector
        int n = X.size();
        long double norm_factor = _l2_norm(X);

        for (int i = 0; i < n; i++) {
            X[i] /= norm_factor;
        }
        return X;
    }

    long double _dot_product(std::vector<long double> X, std::vector<long double> Y) {
        // Calculates the dot product between two vectors, returning a scalar
        long double res = 0.0L;
        int n = X.size();

        for (int i = 0; i < n; i++) {
            res += X[i] * Y[i];
        }

        return res;
    }

    std::vector<long double> _mult_mat_vec(std::vector<std::vector<long double>> A, std::vector<long double> X) {
        // Calculates matrix-vector multiplication, Ax, returning a vector
        int n = X.size();
        std::vector<long double> V(n);

        for (int i = 0; i < n; i++) {
            V[i] = _dot_product(A[i], X);
        }

        return V;
    }

    long double _calc_residual(std::vector<long double> Y, std::vector<long double> X, long double lambda) {
        int n = Y.size();
        std::vector<long double> R(n);

        for (int i = 0; i < n; i++) {
            R[i] = Y[i] - X[i] * lambda;
        }
        return _l2_norm(R);
    }

    long double solver() {
        int n = A.size();

        bool converged = false;

        std::vector<long double> E;
        long double residual;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<long double> dis(1.0, 20.0);

        long double lambda;

        // Initialize a 'guess' for vector X, a vector with random values
        std::vector<long double> X(n);
        for(long double &val : X) val = dis(gen);

        // Initialize a vector Y
        std::vector<long double> Y(n);

        // Power Method
        for (int k = 0; k < maxiter; k++) {
            Y = _mult_mat_vec(A, X);
            lambda = _dot_product(X, Y);

            residual = _calc_residual(Y, X, lambda);
            E.push_back(residual);

            if (residual < tol) {
                converged = true;
                if (verbose == true) {
                    std::cout << "Converged in " << k << " iterations" << std::endl;
                }
                break;
            }

            X = _unit_vec(Y);
        }
        
        if (verbose == true && converged == true) {
            std::cout << "Eigenvalue: " << lambda << std::endl;

            std::cout << "Error Sequence: ";
            for (size_t i = 0; i < E.size(); ++i) {
                std::cout << E[i];
                if (i < E.size() - 1) {
                    std::cout << ", ";
                }
            }
    
        }

        if (converged == true) {
            return lambda;
        } else {
            std::cout << "Newton's method failed to converge in " << maxiter << " iterations" << std::endl;
            return std::numeric_limits<long double>::quiet_NaN();
        }

    }

private:
    std::vector<std::vector<long double>> A;
    int maxiter;
    long double tol;
    bool verbose;

};