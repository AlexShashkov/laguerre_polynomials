#ifndef LAGUERRE
#define LAGUERRE

#include <vector>
#include <limits>
#include <numbers>
#include <complex>
#include <iostream>
#include <algorithm>

#include "Polynomial.h"
#include "ExtendedFunctions.h"

namespace Laguerre{
using std::complex;
using std::vector;

template<typename T>
class Original{
private:
    //
    static constexpr T eps   = std::numeric_limits<T>::epsilon();
    static const int MT = 10;

    /**
     * \brief Laguerre's method to find a root of a polynomial.
        * \param a   Polynomial object.
        * \param x   Initial guess for the root.
     */
    inline void laguer(const std::vector<std::complex<T>>& a, std::complex<T>& x, bool& converged, int m, int itmax){
        static const T frac[9] = {0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0};
        converged = 1;
        std::complex<T> dx, x1, b, d, f, g, h, sq, gp, gm, g2;
        T err, abx, abp, abm;
        std::complex<T> a_m = a[m];

        for (int iter = 1; iter <= itmax; iter++) {
            b = a_m;
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

            if (abp < abm)
                gp = gm;

            dx = (std::max(abp, abm) > static_cast<T>(0.0)) ? (static_cast<T>(m) / gp) : std::polar(static_cast<T>(1) + abx, static_cast<T>(iter));
            x1 = x - dx;

            if (x == x1)
                return;  // Converged.

            if (iter % MT != 0)
                x = x1;
            else
                x = fma(x, -frac[iter / MT], dx);
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
         * \param polish Whether to polish the roots.
         * \param itmax Maximum number of iterations.
     */
    void operator()(Polynomial<T>& poly, std::vector<std::complex<T>>& roots, std::vector<int>& conv, bool polish, int itmax=80){
       std::complex<T> x, _b, _c;
        int m = poly.degree();
        bool conv_status = true;
        std::vector<std::complex<T>> ad(m + 1);

        // Copy of coefficients for successive deflation.
        for (int i = 0; i <= m; i++)
            ad[i] = poly[i];

        for (int j = m - 1; j >= 0; j--) {
            x = std::complex<T>(0.0, 0.0);

            // Start at zero to favor convergence to the smallest remaining root.
            std::vector<std::complex<T>> ad_v(ad.cbegin(), ad.cbegin() + j + 2);
            laguer(ad_v, x, conv_status, m, itmax);
            conv[j] = conv_status ? 1 : -1;

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

        if (polish) {
            for (int j = 1; j < m; j++) {
                x = roots[j];
                int i;

                for (i = j - 1; j < m; j++) {
                    if (real(roots[i]) <= real(x))
                        break;
                    roots[i + j] = roots[i];
                }
                roots[i + j] = x;
            }
        }
    }

};

};
#endif