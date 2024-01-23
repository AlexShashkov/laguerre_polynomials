// Александр, Дмитрий

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
using std::complex;
using std::vector;

template<typename T>
class Original : public BaseSolver<T>{
private:
    //
    static constexpr T eps   = std::numeric_limits<T>::epsilon();

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

            if (std::abs(b) <= err * eps)
                return;  // We are on the root.

            g = d / b;
            g2 = g * g;
            h = fma(f / b, -static_cast<T>(2), g2);
            sq = std::sqrt(static_cast<T>(m - 1) * (fma(h, static_cast<T>(m), -g2)));
            gp = g + sq;
            gm = g - sq;
            abp = std::abs(gp);
            abm = std::abs(gm);

            // den (in this case gp) = max{|G + sq|, |G − sq|}
            gp = abp < abm ? gm : gp; // |G + sq| < |G - sq| ? G − sq : G + sq
            // std::cout << "!!! " << gp << "\n"; 
            // dx = std::complex<T>(m, 0)/gp;
            dx = (std::max(abp, abm) > static_cast<T>(0.0)) ? (static_cast<T>(m) / gp) : std::polar(static_cast<T>(1) + abx, static_cast<T>(iter));
            // std::cout << "??? " << dx << "\n";
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
    void operator()(std::vector<T>& poly, std::vector<std::complex<T>>& roots, std::vector<int>& conv, int itmax=80) override{
        std::complex<T> x, _b, _c;
        int m = roots.size();
        // std::cout << "M SIZE " << m << "\n";
        std::vector<std::complex<T>> ad(m + 1);

        // Copy of coefficients for successive deflation.
        for (int i = 0; i <= m; i++)
            ad[i] = poly[i];

        std::vector<std::complex<T>> ad_v;
        for (int j = m - 1; j > 0; --j) {
            // std::cout << "j:" << j << "\n";
            x = std::complex<T>(0.0, 0.0);

            // Start at zero to favor convergence to the smallest remaining root.
            ad_v = std::vector<std::complex<T>>(ad.cbegin(), ad.cbegin() + j + 2);
            laguer(ad_v, x, conv[j], itmax);

            x.imag(std::fabs(imag(x)) <= std::fabs(real(x)) * eps ? static_cast<T>(0) : x.imag());
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
        
        // Polishing
        if (anycomplex(roots)) {
            std::cout << "has complex!\n";
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
        // std::cout << "SOLVED EVERYTHING\n";
    }

};

};
#endif