#include <iostream>
#include <stdexcept>
#include "bisection.hpp"

double randomFunc(double x) {
    return x * x - 4.0;
}

int main() {
    double root = bisection(randomFunc, 0.0, 3.0, 1e-7, 1000);
    std::cout << "Root: " << root << std::endl;
    return 0;
}