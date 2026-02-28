#include <cmath>
#include <iostream>
#include <stdexcept>

double bisection(double (*f)(double), double a, double b, double tol, int maxiter) {

    if (a >= b || f(a) * f(b) > 0){
        throw std::runtime_error("Root is not bracketed: f(a) and f(b) must have opposite signs and a must be strictly less than b.");
    }

    double error = INFINITY;
    int cnt = 1;

    while (error > tol && cnt <= maxiter){

        double m = (a + b) / 2;
        double f_m = f(m);
        error = std::abs(f_m);

        if (f_m == 0.0 || error <= tol){
            return m;
        }

        if (f(a) * f_m < 0){
            b = m;
        }
        else{
            a = m;
        }

        std::cout << "Iteration: " << cnt << std::endl;
        std::cout << "Error: " << error << std::endl;
        std::cout << "New window: (" << a << " ," << b << ")" << std::endl;
        
        cnt++;


    }

    std::cout << "Bisection Method did not converge in " << maxiter << " iterations." << std::endl;
    return (a + b) / 2;

}

double randomFunc(double x) {
    return x * x - 4.0;
}

int main() {
    double root = bisection(randomFunc, 0.0, 3.0, 1e-7, 1000);
    std::cout << "Root: " << root << std::endl;
    return 0;
}