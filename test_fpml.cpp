#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include "fpml.h"

// Ваша функция здесь

int main() {
    // Тест 1: Полином x^2 - 1 = 0
    std::vector<std::complex<double>> poly1 = {std::complex<double>(1,0), std::complex<double>(1,0), std::complex<double>(1,0),std::complex<double>(1,0), std::complex<double>(1,0),std::complex<double>(1,0), std::complex<double>(1,0),std::complex<double>(1,0) };
    int deg1 = 7;
    std::vector<std::complex<double>> roots1(deg1);
    std::vector<double> berr1(deg1);
    std::vector<double> cond1(deg1);
    std::vector<int> conv1(deg1);
    int itmax1 = 80;
    laguere(poly1, deg1, roots1, berr1, cond1, conv1, itmax1);

    // Проверка результатов
    for(int i = 0; i < deg1; i++) {
        std::cout << "Root " << i+1 << ": " << roots1[i] << ", backward error: " << berr1[i] << ", condition number: " << cond1[i] << ", convergence: " << conv1[i] << std::endl;
    }

    // Добавьте здесь больше тестов...
}
