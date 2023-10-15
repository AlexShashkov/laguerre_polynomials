#include <complex>
#include <vector>
#include <limits>
#include "ExtendedFMA.h"

template <typename T>
using VecT = std::vector<T>;

// Global EPS (epsilon) for convergence
const double EPS = std::numeric_limits<double>::epsilon();

/**
 * Laguerre's method to find a root of a polynomial.
 *
 * @param a   Coefficients of the polynomial.
 * @param x   Initial guess for the root.
 */
template <typename T>
void laguer(const VecT<std::complex<T>>& a, std::complex<T>& x) {
    const int MT = 10;
    static const double frac[9] = {0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0};
    std::complex<T> dx, x1, b, d, f, g, h, sq, gp, gm, g2;
    T err, abx, abp, abm;
    int m = a.size() - 1;
    std::complex<T> a_m = a[m];

    for (int iter = 1; iter <= 80; iter++) {
        b = a_m;
        err = std::abs(b);
        d = f = std::complex<T>(0.0, 0.0);
        abx = std::abs(x);

        for (int j = m - 1; j >= 0; j--) {
            f = implementations::fma(x, f, d);
            d = implementations::fma(x, d, b);
            b = implementations::fma(x, b, a[j]);
            err = fma(err, abx, std::abs(b));
        }

        if (std::abs(b) <= err * EPS)
            return;  // We are on the root.

        g = d / b;
        g2 = g * g;
        h = implementations::fma(f / b, -static_cast<T>(2), g2);
        sq = std::sqrt(static_cast<T>(m - 1) * (implementations::fma(h, static_cast<T>(m), -g2)));
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
            x = implementations::fma(x, -frac[iter / MT], dx);
    }
    throw "Too many iterations in laguer";
}

/**
 * Find all roots of a polynomial using the ZROOTS method.
 *
 * @param a      Coefficients of the polynomial.
 * @param roots  Output vector to store the roots.
 * @param polish Whether to polish the roots.
 * @param _frac  Array of fractional values for iteration (provided to avoid repeated initialization).
 */
template <typename T>
void zroots(const VecT<std::complex<T>>& a, VecT<std::complex<T>>& roots, const bool polish, const double* _frac) {
    std::complex<T> x, _b, _c;
    int m = a.size() - 1;
    VecT<std::complex<T>> ad(m + 1);

    // Copy of coefficients for successive deflation.
    ad = a;

    for (int j = m - 1; j >= 0; j--) {
        x = std::complex<T>(0.0, 0.0);

        // Start at zero to favor convergence to the smallest remaining root.
        std::vector<std::complex<T>> ad_v(ad.cbegin(), ad.cbegin() + j + 2);
        laguer(ad_v, x);

        if (std::abs(imag(x)) <= std::abs(real(x)) * EPS)
            x.imag(static_cast<T>(0));

        roots[j] = x;
        _b = ad[j + 1];

        for (int jj = j; jj >= 0; jj--) {
            _c = ad[jj];
            ad[jj] = _b;
            _b = implementations::fma(x, _b, _c);
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
