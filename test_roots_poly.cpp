#include <iostream>
#include "roots.poly.h"

int main() {
    // Example 1: Polynomial with real coefficients
    VecT<std::complex<double>> a1 = {
            std::complex<double>(52, 0),
            std::complex<double>(-112, 0),
            std::complex<double>(24, 0),
            std::complex<double>(8, 0),
            std::complex<double>(1, 0)
    };
    VecT<std::complex<double>> roots1(a1.size() - 1);

    std::cout << "Example 1: Polynomial with real coefficients" << std::endl;
    zroots(a1, roots1, false, &EPS);

    for (size_t i = 0; i < a1.size(); i++) {
        std::cout << "a[" << i << "]= \t" << a1[i] << '\n';
    }

    std::cout << "Polynomial: ";
    for (size_t i = 0; i < a1.size(); i++) {
        std::cout << a1[i] << "x^" << i;
        if (i < a1.size() - 1) {
            std::cout << '+';
        }
    }

    std::cout << '\n';

    std::cout << "Roots:" << std::endl;
    for (size_t i = 0; i < roots1.size(); i++) {
        std::cout << roots1[i] << '\n';
    }
    std::cout << '\n';

    // Example 2: Another polynomial with real coefficients
    VecT<std::complex<double>> a2 = {
            std::complex<double>(-1000.0, 0),
            std::complex<double>(1000002.0, 0),
            std::complex<double>(-2000.001, 0),
            std::complex<double>(1.0, 0)
    };
    VecT<std::complex<double>> roots2(a2.size() - 1);

    std::cout << "Example 2: Another polynomial with real coefficients" << std::endl;
    zroots(a2, roots2, false, &EPS);

    for (size_t i = 0; i < a2.size(); i++) {
        std::cout << "a[" << i << "]= \t" << a2[i] << '\n';
    }

    std::cout << "Polynomial: ";
    for (size_t i = 0; i < a2.size(); i++) {
        std::cout << a2[i] << "x^" << i;
        if (i < a2.size() - 1) {
            std::cout << '+';
        }
    }

    std::cout << '\n';

    std::cout << "Roots:" << std::endl;
    for (size_t i = 0; i < roots2.size(); i++) {
        std::cout << roots2[i] << '\n';
    }
    std::cout << '\n';

    // Example 3: Polynomial with complex coefficients
    VecT<std::complex<double>> a3 = {
            std::complex<double>(-6.855188152137764e-15, 0),
            std::complex<double>(7.500004043464861e-01, 0),
            std::complex<double>(2.668263562685250e+07, 0),
            std::complex<double>(5.993293694242965e+11, 0),
            std::complex<double>(3.621535214588474e+11, 0)
    };
    VecT<std::complex<double>> roots3(a3.size() - 1);

    std::cout << "Example 3: Polynomial with complex coefficients" << std::endl;
    zroots(a3, roots3, false, &EPS);

    for (size_t i = 0; i < a3.size(); i++) {
        std::cout << "a[" << i << "]= \t" << a3[i] << '\n';
    }

    std::cout << "Polynomial: ";
    for (size_t i = 0; i < a3.size(); i++) {
        std::cout << a3[i] << "x^" << i;
        if (i < a3.size() - 1) {
            std::cout << '+';
        }
    }

    std::cout << '\n';

    std::cout << "Roots:" << std::endl;
    for (size_t i = 0; i < roots3.size(); i++) {
        std::cout << roots3[i] << '\n';
    }
    std::cout << '\n';

    return 0;
}
