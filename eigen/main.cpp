#include <iostream>
#include <stdexcept>
#include <vector>
#include "NewtonRaphson.hpp"


std::vector<std::vector<long double>> A = {
    {3, -1},
    {-1, 3}
};

std::vector<std::vector<long double>> B = {
    {4, 1, 2},
    {1, 3, 0},
    {2, 0, 5}
};

std::vector<std::vector<std::vector<long double>>> matrices = {A, B};


int main() {

    std::cout << "--- Newton's Method for Nonlinear Systems ---" << std::endl;

    for (int i = 0; i < 2; i++) {
        std::cout << "Matrix " << i+1 << ": " << std::endl;

        NewtonRaphson nr(matrices[i], 10000, 1e-15, true);
        nr.solver();

        std::cout << "" << std::endl;
    }
    return 0;
}