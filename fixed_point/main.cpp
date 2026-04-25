#include <iostream>
#include <stdexcept>
#include "NewtonRaphson.hpp"

double f(double x) {
    return x * x * x * x * x * x - 8;
}

double g(double x) {
    return 6 * x * x * x * x * x;
}

int main() {
    NewtonRaphson nr(f, g, 10, 1000, 1e-10, true);
    nr.solver();
    return 0;
}