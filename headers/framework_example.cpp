// Александр Шашков cracecrunch@gmail.com
// Дмитрий Балашов dimabalash0v@yandex.ru

// Example of how to use framework

#include <algorithm>
#include <iostream>
#include <limits>
#include <iomanip>
#include <vector>
#include <cmath>

#include "framework.h"

using std::fma;

// You can check recommended precision for ttmath on https://www.ttmath.org/online_calculator/
// lets use big precision
// FOR 64-bit 
// mantissa = 2048 bits / 64 bits per word = 32 words
// exponent = 256 bits / 64 bits per word = 4 words
// FOR 32-bit
// mantissa = 2048 bits / 32 bits per word = 64 words
// exponent = 256 bits / 32 bits per word = 8 words
constexpr int mantissa = 32; // Mantissa for big number notation
constexpr int exponent = 4; // Exponent for big number notation

int main() {
    std::cout << std::numeric_limits<double>::max() << "\n";
    std::cout << std::numeric_limits<long double>::max() << "\n";

    int l = 10;
    std::vector<long double> roots(l, 0.0);
    std::vector<long double> a(l+1, 0.0);
    generate_polynomial<long double, exponent, mantissa>(l, 0, 0, 0,
        1e-5, -10.0, 10.0, roots, a);

std::cout << "\nFRAMEWORK ROOTS:\n";
    for(auto &el : roots){
        std::cout << "(x-(" << el << "))";
    }

std::cout << "\nFRAMEWORK COEFFS:\n";
    for(auto &el : a){
        std::cout<< el << ',';
    }   

    return 0;
}