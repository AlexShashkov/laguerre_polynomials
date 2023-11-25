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
using std::complex;
using std::vector;

template<typename T>
class Laguerre13 : public BaseSolver<T>{
private:
    //
    static constexpr T frac[9] = {0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0};
    static constexpr T eps   = std::numeric_limits<T>::epsilon();
    static const int MT = 10;
    static const int lambda = 6;
    

    /**
     * \brief Laguerre's method to find a root of a polynomial.
        * \param a   Polynomial object.
        * \param x   Initial guess for the root.
     */
    inline void laguer(const std::vector<std::complex<T>>& a, std::complex<T>& x, int& converged, int itmax){        
        converged = 1;
        std::complex<T> dx, x1, b, d, f, g, h, sq, gp, gm, g2;
        T err, abx, abp, abm;
        int m = a.size() - 1;
        std::cout<<"iteration:\n";
        for (int iter = 1; iter <= itmax; iter++) {
            b = std::complex<T>(1.0, 0.0);
            err = std::abs(b);
            d = f = std::complex<T>(0.0, 0.0);
            abx = std::abs(x);

            for (int j = m - 1; j >= 0; j--) {
                f = fma(x, f, d);
                d = fma(x, d, b);
                b = fma(x, b, a[j]);
                err = fma(err, abx, std::abs(b));
            }
            std::cout<<"b\t" << std::abs(b)<<"\tepsilon\t"<<err*eps<<"\n";
            if (std::abs(b) <= err*eps)
                return;  // We are on the root.

            g = d / b;
            g2 = g * g;
            h = fma(f / b, -static_cast<T>(2), g2);
            sq = std::sqrt(static_cast<T>(lambda-1) * (fma(h, static_cast<T>(lambda), -g2)));
            gp = g + sq;
            gm = g - sq;
            abp = std::abs(gp);
            abm = std::abs(gm);

            gp = abp < abm ? gm : gp;
            dx = (std::max(abp, abm) > static_cast<T>(0.0)) ? (static_cast<T>(lambda) / gp) : std::polar(static_cast<T>(1) + abx, static_cast<T>(iter));
            x1 = x - dx;
            std::cout<<x1<<"\t";
            std::cout<<dx<<"\n";
            if (x == x1)
                return;  // Converged.
            x = iter % MT != 0 ? x1: x1;//fma(dx, -frac[iter / MT], x);
        
        }
        converged = -1;
        // throw "Too many iterations!";
    }


public:
    Laguerre13(){
        // 
    }

    /**
     * \brief Find the roots of a polynomial using Laguerre's method.
         * \param poly Polynomial object.
         * \param roots Vector to store the roots.
         * \param conv Vector to store convergence status of each root.
         * \param itmax Maximum number of iterations.
     */
    void operator()(std::vector<T>& poly, std::vector<std::complex<T>>& roots, std::vector<int>& conv, int itmax=80) override{
        std::complex<T> x, _b, _c;
        int m = poly.size() - 1;
        std::vector<std::complex<T>> ad(m + 1);

        // Copy of coefficients for successive deflation.
        for (int i = 0; i <= m; i++)
            ad[i] = poly[i];

        std::vector<std::complex<T>> ad_v;
        for (int j = m - 1; j > 0; --j) {
            x = std::complex<T>(0.0, 0.0);

            // Start at zero to favor convergence to the smallest remaining root.
            ad_v = std::vector<std::complex<T>>(ad.cbegin(), ad.cbegin() + j + 2);
            laguer(ad_v, x, conv[j], itmax);

            if (std::abs(imag(x)) <= std::abs(real(x)) * eps)
                x.imag(static_cast<T>(0));

            roots[j] = x;
            _b = ad[j + 1];
            for (int jj = j; jj >= 0; jj--) {
                _c = ad[jj];
                ad[jj] = _b;
                _b = fma(x, _b, _c);
            }
        }
        roots[0] = -ad[0];
        conv[0] = conv[1];
        std::cout<< "test:\n";
        for (auto &el : roots){
            std::cout<< el - 1.0 << "\n";
        }

        
       
    }

};

};
#endif
