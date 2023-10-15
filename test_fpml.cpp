#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include "fpml.h"

int main() {
    // Example usage of the laguerre function
    int degree = 3;
    std::vector<std::complex<double>> polynomial(degree + 1);
    std::vector<std::complex<double>> roots(degree);
    std::vector<double> berr(degree);
    std::vector<double> cond(degree);
    std::vector<int> conv(degree);

    // Initialize polynomial coefficients with complex values
    polynomial[0] = 1.0;
    polynomial[1] = -6.0;
    polynomial[2] = 11.0;
    polynomial[3] = -6.0;

    // Call the laguerre function to find the roots
    laguerre(polynomial, degree, roots, berr, cond, conv, 100);

    // Output the roots and associated information
    for (int i = 0; i < degree; i++) {
        std::cout << "Root " << i + 1 << ": " << roots[i] << std::endl;
        std::cout << "Convergence Status: " << conv[i] << std::endl;
        std::cout << "Backward Error: " << berr[i] << std::endl;
        std::cout << "Condition Number: " << cond[i] << std::endl;
    }

    return 0;
}