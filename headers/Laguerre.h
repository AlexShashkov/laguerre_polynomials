// Александр Шашков cracecrunch@gmail.com
// Дмитрий Балашов dimabalash0v@yandex.ru

// https://ora.ox.ac.uk/objects/uuid:a1cb2db7-b587-4a0f-a442-dc3e06797726
// Based on original code: https://github.com/andresmmera/Laguerre-Root-Solver

#ifndef LAGUERRE
#define LAGUERRE

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
class Original : public BaseSolver<T>{
private:
    /**
     * \brief Laguerre's method to find a root of a polynomial.
        * \param a   Polynomial object.
        * \param x   Initial guess for the root.
     */
    inline void laguer(const vector<complex<T>>& a, complex<T>& x, int& converged, int itmax){        
        converged = 1;
        complex<T> dx, x1, b, d, f, g, h, sq, gp, gm, g2;
        T err, abx, abp, abm;
        int m = a.size() - 1;

        for (int iter = 1; iter <= itmax; ++iter) {
            b = complex<T>(1.0, 0.0); // b will be used to calculate the value of the polynomial at the point x
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

            g = d / b;
            g2 = g * g;
            h = fma(f / b, -static_cast<T>(2), g2); // Calculation of h, which is used in the Lagger formula

            // if(complexnotfinite(g, BaseSolver<T>::big) || complexnotfinite(h, BaseSolver<T>::big)){
            //     throw std::invalid_argument("During error calculation in Laguerre some numbers converged to NaN.");
            //}

            sq = std::sqrt(static_cast<T>(m - 1) * (fma(h, static_cast<T>(m), -g2))); // square root
            // two possible denominator values
            gp = g + sq;
            gm = g - sq;
            abp = abs(gp);
            abm = abs(gm);

            // den (in this case gp) = max{|G + sq|, |G − sq|}
            gp = abp < abm ? gm : gp; // |G + sq| < |G - sq| ? G − sq : G + sq
            // std::cout << "!!! " << gp << "\n"; 
            // dx - change of x on each step
            dx = (std::max(abp, abm) > static_cast<T>(0.0)) ? (static_cast<T>(m) / gp) : std::polar(static_cast<T>(1) + abx, static_cast<T>(iter));
            // if(complexnotfinite(dx, BaseSolver<T>::big)){
            //    throw std::invalid_argument("During error calculation in Laguerre some numbers converged to NaN.");
            //}
            x1 = x - dx;
            if (x == x1)
                return;  // Converged.
            x = x1;      
        }
        converged = -1;
        // throw "Too many iterations!";
    }


public:
    Original(){
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
        int m = roots.size();
        // std::cout << "M SIZE " << m << "\n";
        vector<complex<T>> ad(m + 1);

        // Copy of coefficients for successive deflation.
        // ad = poly; no match operator between complex vector and double vector
        for (int i = 0; i <= m; ++i)
             ad[i] = complex<T>(poly[i], static_cast<T>(0));

        vector<complex<T>> ad_v;
        for (int j = m - 1; j > 0; --j) {
            // std::cout << "j:" << j << "\n";
            x = complex<T>(0.0, 0.0);

            // Start at zero to favor convergence to the smallest remaining root.
            ad_v = vector<complex<T>>(ad.cbegin(), ad.cbegin() + j + 2);
            laguer(ad_v, x, conv[j], itmax);

            x.imag(std::fabs(imag(x)) <= std::fabs(real(x)) * BaseSolver<T>::eps ? static_cast<T>(0) : x.imag());
            roots[j] = x;
            _b = ad[j + 1]; 
            for (int jj = j; jj >= 0; --jj) {
                _c = ad[jj];
                ad[jj] = _b;
                _b = fma(x, _b, _c);
            }
        }
        roots[0] = -ad[0];
        conv[0] = conv[1];
        
        // Polishing
        if (anycomplex<T>(roots)) {
            // return;
            int i;
            for (int j = 1; j < m; ++j) {
                x = roots[j];
                for (i = j - 1; i>=0; --i) {
                    if (real(roots[i]) <= real(x))
                        break;
                    roots[i + 1] = roots[i];
                }
                roots[i + 1] = x;
            }
        }
    }

};

};
#endif