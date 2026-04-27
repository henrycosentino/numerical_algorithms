#include <iostream>
#include <stdexcept>
#include <vector>
#include "NewtonsMethod.hpp"


std::vector<std::vector<long double>> A = {
    {3, -1},
    {-1, 3}
};

std::vector<std::vector<long double>> Q = {
    // Q matrix for a Markov Process, eigenvalue of 0 is guranteed!
    {-3, 3, 0},
    {3, -5, 2},
    {0, 2, -2}
};

std::vector<std::vector<std::vector<long double>>> matrices = {A, Q};


int main() {

    std::cout << "--- Newton's Method for Nonlinear Systems ---" << std::endl;

    for (int i = 0; i < 2; i++) {
        std::cout << "Matrix " << i+1 << ": " << std::endl;

        NewtonsMethod nm(matrices[i], 1000, 1e-15, true);
        nm.solver();

        std::cout << "" << std::endl;
    }
    return 0;
}