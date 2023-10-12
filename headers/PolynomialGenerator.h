#ifndef POLYNOMIALGENERATOR_H
#define POLYNOMIALGENERATOR_H

#include <vector>
#include <iostream>
#include "Polynomial.h"

namespace Laguerre{

template<typename T>
class Generators{
private:
    static Polynomial<T> generate(std::vector<T> roots){
        int count = roots.size();
        if(count == 1){
            return Polynomial<T>(-roots[0], 1);
        }
        else{
            int mid = count / 2;
            auto mid_begin = roots.begin() + mid;
            auto end = roots.end();
            Polynomial left_poly = generate(std::vector<T>(mid_begin -mid, mid_begin));
            Polynomial right_poly = generate(std::vector<T>(mid_begin, end));
            return left_poly*right_poly;
        }
    }
public:
    static Polynomial<T> createFromRoots(std::vector<T> args){
        size_t count = args.size();
        if(count == 0)
            throw std::invalid_argument("Cant create polynomial from zero given roots.");
        return generate(args);
    }

    template<typename... Args>
    static Polynomial<T> createFromRoots(Args...args){
        return createFromRoots(std::vector<T>{args...});
    }

    Generators() = delete;
};

}

#endif
