#ifndef POLYNOMIALGENERATOR_H
#define POLYNOMIALGENERATOR_H

#include <random>
#include <chrono>

#include <vector>
#include <iostream>

#include "Polynomial.h"
#include "ExtendedFunctions.h"

namespace Laguerre{

template<typename T>
class Generator{
protected:
    // Pointer to a coeff generator function
    typedef Polynomial<T> (*coeffGeneratorPtr)(std::vector<T>);
    static coeffGeneratorPtr coeffGenerator;

    // Pointer to a roots generator function
    typedef std::vector<T> (*rootsGeneratorPtr)(int, T, T, T);
    static rootsGeneratorPtr rootsGenerator;

    static Polynomial<T> basicCoeffsGeneration(std::vector<T> roots){
        int count = roots.size();
        if(count == 1){
            return Polynomial<T>(-roots[0], static_cast<T>(1));
        }
        else{
            int mid = count / 2;
            auto mid_begin = roots.begin() + mid;
            auto end = roots.end();
            Polynomial left_poly = basicCoeffsGeneration(std::vector<T>(mid_begin -mid, mid_begin));
            Polynomial right_poly = basicCoeffsGeneration(std::vector<T>(mid_begin, end));
            return left_poly*right_poly;
        }
    }

    static std::vector<T> basicRootGeneration(int count = 3, T from = 0, T to = 1, T delta=1){
        std::uniform_real_distribution<double> distribution(from, to);
        std::random_device rd;
        std::default_random_engine generator(rd());
        std::vector<T> roots;
        for(int i = 0; i < count; ++i)
            roots.push_back(distribution(generator)*delta);
        return roots;
    }

    // Returns pointer to a function that generates coeffs
    static coeffGeneratorPtr currentCoeffGenerator() {
        if (coeffGenerator == nullptr) {
            coeffGenerator = &basicCoeffsGeneration;
        }
        return coeffGenerator;
    }

    // Returns pointer to a function that generates roots
    static rootsGeneratorPtr currentRootsGenerator() {
        if (rootsGenerator == nullptr) {
            rootsGenerator = &basicRootGeneration;
        }
        return rootsGenerator;
    }

public:
    static void setCoeffGenerator(coeffGeneratorPtr func){
        coeffGenerator = func;
    }

    static void setRootsGenerator(rootsGeneratorPtr func){
        rootsGenerator = func;
    }

    static std::vector<T> generateRoots(int count = 3, T from = 0, T to = 1, T delta=1){
        return currentRootsGenerator()(count, from, to, delta);
    }

    static Polynomial<T> createFromRoots(std::vector<T> args){
        if(!args.size())
            throw std::invalid_argument("Cant create polynomial from zero given roots.");
        auto res = currentCoeffGenerator()(args);
        res.setRoots(args);
        return res;
    }

    template<typename... Args>
    static Polynomial<T> createFromRoots(Args...args){
        return createFromRoots(std::vector<T>{args...});
    }

    static Polynomial<T> gen(int count = 3, T from = 1, T to = 3, T delta=1){
        std::vector<T> roots = currentRootsGenerator()(count, from, to, delta);
        Polynomial<T> res = createFromRoots(roots);
        res.setRoots(roots);
        return res;
    }
};

template<typename T>
typename Laguerre::Generator<T>::coeffGeneratorPtr Laguerre::Generator<T>::coeffGenerator = nullptr;

template<typename T>
typename Laguerre::Generator<T>::rootsGeneratorPtr Laguerre::Generator<T>::rootsGenerator = nullptr;

}

#endif
