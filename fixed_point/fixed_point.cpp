#include <cmath>
#include <vector>
#include <iostream>
#include <stdexcept>

using namespace std;

class FixedPoint
{
public:
    FixedPoint(double (*f)(double), double (*g)(double), double p, int maxiter = 1000, double tol = 1e-10, bool diagnostics = false) {
        this->f = f;                     // Polynomial function
        this->g = g;                     // Derivative of the polynomial function
        this->p = p;                     // Starting estimate of the fixed point
        this->maxiter = maxiter;         // The maximum amount of iterations for the fixed point algorithm
        this->tol = tol;                 // Error tolerance   
        this->diagnostics = diagnostics; // Convergence diagnostics   
    }

    double fixedPoint() {
        /*
        Method that estimates the fixed point.
        */
        vector<double> p_sequence;

        for (int i = 0; i < maxiter; i++) {

            // Assign current estimate to the sequence
            p_sequence.push_back(p);

            // Calculate the next estimate
            p = p - f(p) / g(p);

            // Early convergence check
            if (abs(p - p_sequence[i]) < tol) {
                if (diagnostics)
                    std::cout << "Fixed point iteration converged early in " << i << " iterations: " << p << "\n";
                return p;
            }

            // Dynamic error analysis
            if (i >= 2) {
                double alpha = log(abs(p - p_sequence[i]) / abs(p_sequence[i] - p_sequence[i-1])) /
                               log(abs(p_sequence[i] - p_sequence[i-1]) / abs(p_sequence[i-1] - p_sequence[i-2]));

                double error;

                if (abs(alpha - 1) < 0.1) {
                    // Case #1: Linear Convergence
                    double lambda = (p_sequence[i] - p_sequence[i-1]) / (p_sequence[i-1] - p_sequence[i-2]);
                    error = lambda / (lambda - 1) * (p_sequence[i] - p_sequence[i-1]);

                    if (abs(error) < tol) {
                        if (diagnostics) {
                            std::cout << "Fixed point iteration converged linearly in " << i << " iterations.\n"; 
                            std::cout << "Fixed point: " << p << "\n";
                            std::cout << "Order of convergence: " << alpha << "\n";
                        }
                        return p;
                    }
                } else {
                    // Case #2: Superlinear Convergence
                    error = p - p_sequence[i];

                    if (abs(error) < tol) {
                        if (diagnostics) {
                            std::cout << "Fixed point iteration converged superlinearly in " << i << " iterations.\n"; 
                            std::cout << "Fixed point: " << p << "\n";
                            std::cout << "Order of convergence: " << alpha << "\n";
                        }
                        return p;
                    }
                }
            }
        }
        
        if (diagnostics)
            std::cout << "Fixed point iteration did not converge in " << maxiter << " iterations.\n";

        return 0;
    }

private:
    double (*f)(double);
    double (*g)(double);
    double p;
    int maxiter;
    double tol;
    bool diagnostics;
};


double f(double x) {
    return x * x * x * x * x * x - 8;
}
double g(double x) {
    return 6 * x * x * x * x * x;
}


int main() {
    FixedPoint FE(f, g, 2, 1000, 1e-10, true);
    FE.fixedPoint();
    return 0;
}