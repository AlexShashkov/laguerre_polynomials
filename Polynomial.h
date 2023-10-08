#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <cmath>
#include <numbers>
#include <vector>
#include <complex>
#include <iostream>
#include <algorithm>
#include <iomanip>

using namespace std::complex_literals;

namespace PolyGen {

template<typename T>
class Polynomial{
private:
    std::vector<T> roots;
    std::vector<T> coeffs;
    int degree = 0;

public:
    Polynomial(std::vector<T> args){
        coeffs = args;
        degree = coeffs.size() - 1;
    }

    template<typename... Args>
    Polynomial(Args...args){
        coeffs = {args...};
        degree = coeffs.size() - 1;
    }

    void setRoots(std::vector<T> args){
        roots = args;
    }

    template<typename... Args>
    void setRoots(Args...args){
        if(!sizeof...(Args)){
            throw std::invalid_argument("Polynomial cant have zero roots!\n");
        }
        if(coeffs.size() != 0){
            if(sizeof...(Args) != degree)
                throw std::invalid_argument("Count of given roots is different from given polynomial degree.\n");
        }
    }

    std::vector<T> diff(int deg){

    }

    void print() {
        for(int i = 0; i < coeffs.size(); ++i, --degree) {
            std::cout << "(" << coeffs[i] << ")";
            if(degree > 0){
                std::cout << "x^" << degree;
                if(i != coeffs.size() - 1) {
                    std::cout << " + ";
                }
            }
        }
        std::cout << '\n';
    }
};
}

#endif
