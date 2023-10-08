#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <iostream>

namespace PolyGen {

template<typename T>
class Polynomial{
private:
    std::vector<T> roots;
    std::vector<T> coeffs;
    int degree = 0;
public:
    Polynomial(std::vector<T> args){
        setCoeffs(args);
    }

    template<typename... Args>
    Polynomial(Args...args){
        setCoeffs(args...);
    }

    void setRoots(std::vector<T> args){
        size_t count = args.size();
        if(args.size() != degree)
            throw std::invalid_argument("Count of given roots is different from given polynomial degree. (" 
                    + std::to_string(count) + " vs " + std::to_string(degree) + ")\n");
        roots = args;
        }

    template<typename... Args>
    void setRoots(Args...args){
        setRoots(std::vector<T>{args...});
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

    std::vector<T> diff(int deg){

    }

    void print() {
        if(!degree){
            std::cout << "This polynomial object doesnt contain any coefficient.\n";
            return;
        }
        int deg = degree;
        for(int i = 0; i < coeffs.size(); ++i, --deg) {
            std::cout << "(" << coeffs[i] << ")";
            if(deg > 0){
                std::cout << "x^" << deg;
                if(i != coeffs.size() - 1) {
                    std::cout << " + ";
                }
            }
        }
        std::cout << ". Degree: " << degree << '\n';
    }
};
}

#endif
