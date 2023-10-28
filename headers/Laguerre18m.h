#ifndef LAGUERRE18
#define LAGUERRE18

#include <vector>
#include <complex>
#include <iostream>
#include <algorithm>

#include "Polynomial.h"
#include "ExtendedFunctions.h"

namespace Laguerre{
using std::complex;
using std::vector;
using std::fma;

template<typename T>
class ModifiedLaguerre18{
private:
    static constexpr number PI = std::numbers::pi_v<number>;

    /** \brief Find starting approximation
		 * 	\param poly Polynomial object
    */
    inline int orientationAndrew(vector<T> P){
        T val = 
    }

    /** \brief Find starting approximation
		 * 	\param poly Polynomial object
    */
    inline vector<T> startAndrew(Polynomial<T> poly){
        // C = {(i, log |ai|), i = 0, 1, ..., m}
        std::vector<T> C(poly.degree());
        std::transform(C.begin(), C.end(), C.begin(),
           &{
               return std::log(std::fabs(poly[i]));
           }
        );

    }
public:
    ModifiedLaguerre18(){}

    /** \brief Functor for solving polynomials of any degree
		 * 	\param coeffs Polynomial object that stores coefficcients 
		 * 	\param roots Reference on vector, where the result roots will be stored
		 * 	\param converged Reference on boolean to check if Laguerre converged on roots
    */
    void operator()(Polynomial<T> coeffs, vector<complex<T>> &roots, bool &converged){


    }

};

};
#endif
