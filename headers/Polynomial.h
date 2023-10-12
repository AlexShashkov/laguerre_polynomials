#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <complex>
#include <iostream>

#include "Polynomial.h"

namespace Laguerre{

template<typename T>
class Polynomial{
protected:
    std::vector<std::complex<T>> roots;
    std::vector<T> coeffs;
    int degree = 0;
public:
    using complex = std::complex<T>;

    Polynomial(std::vector<T> _coeffs, std::vector<complex> _roots){
        setCoeffs(_coeffs);
        setRoots(_roots);
    }

    Polynomial(std::vector<T> args){
        setCoeffs(args);
    }

    template<typename... Args>
    Polynomial(Args...args){
        setCoeffs(args...);
    }

    void setRoots(std::vector<complex> args){
        size_t count = args.size();
        if(count != degree)
            throw std::invalid_argument("Count of given roots is different from given polynomial degree. (" 
                    + std::to_string(count) + " vs " + std::to_string(degree) + ")\n");
        roots = args;
    }

    template<typename... Args>
    void setRoots(Args...args){
        setRoots(std::vector<complex>{args...});
    }

    void setCoeffs(std::vector<T> args){
        coeffs = args;
        roots.clear();
        int size = coeffs.size();
        degree = size > 0 ? size - 1 : 0;
    }

    template<typename... Args>
    void setCoeffs(Args...args){
        setCoeffs(std::vector<T>{args...});
    }

    Polynomial& operator=(const Polynomial& other) {
        if (this != &other) {
            this->coeffs = other.coeffs;
        }
        return *this;
    }

    Polynomial operator*(const Polynomial& other) {
        int n = coeffs.size();
        int m = other.coeffs.size();
        std::vector<T> result(n + m - 1, 0);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                result[i + j] += coeffs[i] * other.coeffs[j];
            }
        }
        return Polynomial(result);
    }

    Polynomial& operator*=(const Polynomial& other) {
        *this = *this * other;
        this->roots.clear();
        return *this;
    }

    Polynomial diff(int deg=1){
        Polynomial<T> diff_poly = this;
        for(int j = 0; j < deg; ++j){
            for (int i = 1; i < coeffs.size(); ++i) {
                diff_poly.coeffs[i] *= i;
            }
            if (!diff_poly.coeffs.empty()) {
                diff_poly.coeffs.erase(coeffs.begin());
            }
            else break;
        }
    }

    void print() {
        if(!degree){
            std::cout << "This polynomial object doesnt contain any coefficient.\n";
        }
        else{
            for(int i = degree; i >= 0; --i){
                std::cout << "(" << coeffs[i] << ")";
                std::cout << "x^" << i;
                if(i != 0) {
                    std::cout << " + ";
                }
            }
            std::cout << ". Degree: " << degree << '\n';
        }
    }       
};
};

#endif
