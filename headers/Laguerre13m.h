// Александр Шашков cracecrunch@gmail.com
// Дмитрий Балашов dimabalash0v@yandex.ru

// On a family of Laguerre methods to find multiple roots of nonlinear equations
// https://www.academia.edu/29512423/On_a_family_of_Laguerre_methods_to_find_multiple_roots_of_nonlinear_equations

#ifndef LAGUERRE13
#define LAGUERRE13

#include <vector>
#include <limits>
#include <numbers>
#include <complex>
#include <iostream>
#include <algorithm>

#include "BaseSolver.h"
#include "ExtendedFunctions.h"

namespace Laguerre{

template<typename T>
class ModifiedLaguerre13 : public BaseSolver<T>{
private:
    /**
     * \brief Laguerre's method to find a root of a polynomial, multiple roots modification.
        * \param a   Polynomial object.
        * \param x   Initial guess for the root.
     */
     inline void laguer13(const vector<complex<T>>& a, complex<T>& x, int& lam){
        complex<T> dx, x1, b, d, f, g, h, sq, gp, gm, g2;
        T err, abx, abp, abm;
        int m = a.size() - 1;

        b = complex<T>(1.0, 0.0);
        err = abs(b);
        d = f = complex<T>(0.0, 0.0);
        abx = abs(x);

        for (int j = m - 1; j >= 0; --j) {
            f = fma(x, f, d); // Calculation of the second derivative of the polynomial at the point x
            d = fma(x, d, b); // Calculation of the first derivative of the polynomial at the point x
            b = fma(x, b, a[j]); // Calculating the value of the polynomial at the point x
            err = fma(err, abx, abs(b));
        }

        if (abs(b) <= err * BaseSolver<T>::eps)
            return;  // We are on the root.

        g = static_cast<T>(lam)* b / d; 
        h = static_cast<T>(2)*f / d;

        if(complexnotfinite(g, BaseSolver<T>::big) || complexnotfinite(h, BaseSolver<T>::big)){
            throw invalid_argument("During error calculation in Laguerre some numbers converged to NaN.");
        }
        sq = std::sqrt(fma(-g, h, static_cast<T>(lam - 1)));
        dx = g/abs(static_cast<T>(1) + sq);   
        if(complexnotfinite(dx, BaseSolver<T>::big)){
            throw invalid_argument("During error calculation in Laguerre some numbers converged to NaN.");
        }  
        x -= dx;
        return;  
            
    }


public:
    ModifiedLaguerre13(){
        // 
    }

    /**
     * \brief Find the roots of a polynomial using Laguerre's method.
         * \param poly Polynomial object.
         * \param roots Vector to store the roots.
         * \param conv Vector to store convergence status of each root.
         * \param itmax Maximum number of iterations.
     */
     void operator()(const vector<T>& poly, vector<complex<T>>& roots, vector<int>& conv, int itmax=80) override{
        complex<T> x, _b, _c;
        int m = poly.size() - 1;
        vector<complex<T>> ad(m + 1);
        vector<complex<T>> ad_v;
        int lambda;
        int m1;

        // Copy of coefficients for successive deflation.
        // ad = poly; no match operator between complex vector and double vector
        for (int i = 0; i <= m; ++i)
             ad[i] = complex<T>(poly[i], static_cast<T>(0));
        
        auto diff = Laguerre::diff(poly, 1);
        vector<T> remainder;
        Laguerre::getRemainderFromBigNum(poly, diff, remainder);

        bool zeros = std::all_of(remainder.begin(), remainder.end(), [](int i) { return i==0; });
        if (zeros)
        {
            lambda = 2*m;
            m1 = 1;
        }
        else
        {
            vector<T> remainder1;
            Laguerre::getRemainderFromBigNum(diff, remainder, remainder1);
            zeros = std::all_of(remainder1.begin(), remainder1.end(), [](int i) { return i==0; });
            if (zeros){
                m1 = m - remainder.size() + 1;
                lambda = 2*m/m1;
                if(anynotfinite(lambda))
                    throw invalid_argument("Lambda value is not finite.");
            } 
            else{
                // If we got further, then this is a wrong type of polynomial, thus
                // this method may fail to converge
                std::cout << "\nFAILED TO DETERMINE DEGREE OF THE ROOT! METHOD MAY FAIL TO CONVERGE!\n";
                lambda = 2*m;
                m1 = 1;
            }
        }

        int j = m - 1;
        for (int i = 0; i < m1; ++i){
            x = complex<T>(0.0, 0.0);

            // Start at zero to favor convergence to the smallest remaining root.
            ad_v = vector<complex<T>>(ad.cbegin(), ad.cbegin() + j + 2);
            laguer13(ad_v, x, lambda);
            if (abs(imag(x)) <= abs(real(x)) * BaseSolver<T>::eps)
                x.imag(static_cast<T>(0));

            for (int k = 0; k < lambda/2; ++k){
                roots[j] = x;
                conv[j] = 1;
                _b = ad[j + 1];
                for (int jj = j; jj >= 0; --jj) {
                    _c = ad[jj];
                    ad[jj] = _b;
                    _b = fma(x, _b, _c);
                }
                --j;    
            }     
        }
 }

};

};
#endif
