#ifndef LAGUERRE18
#define LAGUERRE18

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
using std::fma;

template<typename T>
class ModifiedLaguerre18 : public BaseSolver<T>{
private:
    //
    static constexpr T eps   = std::numeric_limits<T>::epsilon();
    static constexpr T big   = std::numeric_limits<T>::max();
    static constexpr T small = std::numeric_limits<T>::min();
    //
    static constexpr T PI    = std::numbers::pi_v<T>;
    static constexpr T pi2   = PI*static_cast<T>(2);


    /**
     * \brief Calculate the determinant between two points.
         * \param h Vector of integers representing points.
         * \param a Vector of Ts representing points' coordinates.
         * \param c Index for the current point.
         * \param i Index for the next point.
         * \return T Determinant between two points.
     */
    inline T cross(std::vector<int>& h, std::vector<T>& a, int c, int i) {
        // determinant
        int h_c = h[c];
        int h_c1 = h[c-1];
        T a_hc = a[h_c1];
        return fms(a[i]-a_hc, static_cast<T>(h_c-h_c1), a[h_c]-a_hc, static_cast<T>(i-h_c1));
    }

    /**
     * \brief Calculate the convex hull of a set of points. Dont forget to calculate Andrews alghoritm to vector a!
         * \param n Number of points.
         * \param a Vector of templated type T representing points coordinates.
         * \param h Vector to store the convex hull.
         * \param c Number of points in the convex hull.
     */
    inline void conv_hull(int n, std::vector<T>& a, std::vector<int>& h, int& c) {
        c = 0;

        for (int i = n - 1; i >= 0; --i) {
            while (c >= 2 && cross(h, a, c, i) < eps)
                --c;
            ++c;
            h[c] = i;
        }
    }

    /**
     * \brief Estimate roots of a polynomial.
         * \param alpha Polynomial object.
         * \param roots Vector to store the estimated roots.
         * \param conv Vector to store convergence status of each root.
         * \param nz Number of zeros.
     */
    inline void estimates(vector<T>& alpha, int deg, std::vector<std::complex<T>>& roots, std::vector<int>& conv, int& nz) {
        int c, i, j, k, nzeros;
        T a1, a2, ang, r, th;
        std::vector<int> h(deg + 1);
        std::vector<T> a(deg + 1);

        // Andrews starting position
        for (i = 0; i < deg + 1; ++i)
            a[i] = alpha[i] > 0 ? log(alpha[i]) : small;
        conv_hull(deg + 1, a, h, c);

        k = 0;
        th = pi2/static_cast<T>(deg);

        // Initial Estimates
        for (i = c - 1; i >= 0; --i) {
            nzeros = h[i] - h[i + 1];
            T one_div_nzeros = static_cast<T>(1.0 / nzeros);
            a1 = pow(alpha[h[i + 1]], one_div_nzeros);
            a2 = pow(alpha[h[i]], one_div_nzeros);

            if (a1 <= a2 * small){
                // r is too small
                r = 0.0;
                nz += nzeros;
                for (j = k; j < k + nzeros; ++j)
                    conv[j] = -1;
                    roots[j] = {0.0, 0.0};
            }
            else{
                T sigma_thh = fma(th, h[i], static_cast<T>(0.7));
                if (a1 >= a2 * big){
                    // r is too big
                    r = big;
                    nz += nzeros;
                    for (j = k; j < k + nzeros; ++j)
                        conv[j] = -1;
                } 
                else{
                    // r is ok
                    r = a1 / a2;
                }
                ang = pi2 / static_cast<T>(nzeros);
                for (j = 0; j < nzeros; ++j){
                    T arg = fma(ang, j, sigma_thh);
                    roots[k + j] = r * std::complex<T>(cos(arg), sin(arg));
                }
            }
            k += nzeros;
        }
    }


    /**
     * \brief Perform Laguerre correction, compute backward error, and condition.
         * \param p Polynomial object.
         * \param alpha Vector of coefficients.
         * \param b Output complex number representing root.
         * \param c Output complex number representing correction.
         * \param z Complex number used for correction.
         * \param r Value used for correction.
         * \param conv Output value indicating convergence status.
         * \param berr Output backward error.
         * \param cond Output condition number.
     */
    inline void rcheck_lag(std::vector<T>& p, vector<T>& alpha, int deg, complex<T>& b, complex<T>& c, complex<T> z, T r, int& conv, T& berr, T& cond) {
        int k;
        T rr;
        complex<T> a, zz, zz2;

        // Evaluate polynomial and derivatives
        zz = complex<T>(1, 0) / z;
        zz2 = zz*zz;
        rr = static_cast<T>(1) / r;
        a = p[0];
        b = 0;
        c = 0;
        berr = alpha[0];

        for (k = 1; k < deg + 1; ++k) {
            c = fma(zz, c, b);
            b = fma(zz, b, a);
            a = fma(zz, a, p[k]);
            berr = fma(rr, berr, alpha[k]);
        }

        // Laguerre correction/ backward error and condition
        bool condition = abs(a) > eps * berr;
        b = condition ? b/a : b;
        c = condition ? fms(zz2, fma(-static_cast<T>(2)*zz, b, static_cast<T>(deg)), zz2*zz2, fms(complex<T>(2.0, 0), c / a, b, b)) : c;
        b = condition ? fms(zz, complex<T>(deg, 0) , zz2, b) : b;

        cond = (!condition)*(berr / fabs(fms(complex<T>(deg, 0),  a, zz, b))) + (condition)*cond;
        berr = (!condition)*(fabs(a) / berr) + (condition)*berr;

        conv = condition ? 
        (complexnotfinite(b, big) || complexnotfinite(c, big) ? -1 : conv) 
        : 1;

        /*
        if (fabs(a) > eps * berr) {
            b = b / a;
            c = fms(zz2, fma(-static_cast<T>(2)*zz, b, static_cast<T>(deg)), zz2*zz2, fms(complex<T>(2.0, 0), c / a, b, b));
            b = fms(zz, complex<T>(deg, 0) , zz2, b);

            if (complexnotfinite(b, big) || complexnotfinite(c, big))
                conv = -1;
        } else {
            cond = berr / fabs(fms(complex<T>(deg, 0),  a, zz, b));
            berr = fabs(a) / berr;
            conv = 1;
        }
        */
    }

    /**
     * \brief Perform Laguerre correction, compute backward error, and condition.
         * \param p Vector of polynomial coefficients.
         * \param alpha Vector of coefficients.
         * \param deg Degree of the polynomial.
         * \param b Output complex number representing root.
         * \param c Output complex number representing correction.
         * \param z Complex number used for correction.
         * \param r Value used for correction.
         * \param conv Output value indicating convergence status.
         * \param berr Output backward error.
         * \param cond Output condition number.
     */
    inline void check_lag(std::vector<T>& p, std::vector<T>& alpha, int deg, std::complex<T>& b, std::complex<T>& c, std::complex<T> z, T r, int& conv, T& berr, T& cond) {
        int k;
        std::complex<T> a;

        // Evaluate polynomial and derivatives
        a = p[deg];
        b = 0;
        c = 0;
        berr = alpha[deg];

        for (k = deg - 1; k >= 0; --k) {
            c = fma(z, c, b);
            b = fma(z, b, a);
            a = fma(z, a, p[k]);
            berr = fma(r, berr, alpha[k]);
        }

        // Laguerre correction/ backward error and condition
        bool condition = abs(a) > eps * berr;
        b = condition ? b/a : b;
        c = condition ? fms(b, b, std::complex<T>(2, 0), c/a) : c;

        cond = (!condition)*(berr / (r * fabs(b))) + (condition)*cond;
        berr = (!condition)*(fabs(a) / berr) + (condition)*berr;

        conv = condition ? 
        (complexnotfinite(b, big) || complexnotfinite(c, big) ? -1 : conv) 
        : 1;

        /*
        if (abs(a) > eps * berr) {
            b = b / a;
            c = fms(b, b, std::complex<T>(2, 0), c/a);

            if (complexnotfinite(b, big) || complexnotfinite(c, big))
                conv = -1;
        } else {
            cond = berr / (r * fabs(b));
            berr = fabs(a) / berr;
            conv = 1;
        }
        */
    }


    /**
     * \brief Modify Laguerre correction based on Aberth's method.
         * \param deg Degree of the polynomial.
         * \param b Complex number representing root.
         * \param c Complex number representing correction.
         * \param z Complex number used for correction.
         * \param j Index of the root to modify.
         * \param roots Vector of estimated roots.
     */
    void modify_lag(int deg, std::complex<T>& b, std::complex<T>& c, std::complex<T> z, int j, std::vector<std::complex<T>>& roots) {
        std::complex<T> t;

        // Aberth correction terms
        for (int k = 0; k < j - 1; ++k) {
            t = 1.0 / (z - roots[k]);
            b -= t;
            c = fma(-t, t, c);
        }
        // k != j
        for (int k = j + 1; k <= deg; ++k) {
            t = 1.0 / (z - roots[k]);
            b -= t;
            c = fma(-t, t, c);
        }

        // Laguerre correction
        t = sqrt(fms(complex<T>(static_cast<T>(deg - 1), 0), (static_cast<T>(deg) * c),  b, b));
        c = b + t;
        b -= t;
        std::complex<T> degc(deg, 0);
        c = abs(b) > abs(c) ? degc / b : degc / c;
    }

public:
    ModifiedLaguerre18(){
        // 
    }

    /**
     * \brief Find the roots of a polynomial using Laguerre's method.
         * \param poly Polynomial object.
         * \param roots Vector to store the roots.
         * \param conv Vector to store convergence status of each root.
         * \param itmax Maximum number of iterations.
     */
    void operator()(std::vector<T>& poly, std::vector<std::complex<T>>& roots, std::vector<int>& conv, int itmax) override{
        int deg = poly.size() - 1;
        std::vector<T> berr(deg); // Vector to store the backward errors.
        std::vector<T> cond(deg); // Vector to store the condition numbers.
        int i, j, nz;
        T r;
        std::vector<T> alpha(deg + 1);
        std::complex<T> b, c, z;

        // Precheck
        for (i = 0; i <= deg; i++)
            alpha[i] = abs(poly[i]);

        if (alpha[deg] < small) {
            // TODO: exception?
            std::cout << "Warning: leading coefficient is too small." << std::endl;
            return;
        } 

        conv.assign(deg, 0);
        nz = 0;
        estimates(alpha, deg, roots, conv, nz);

        for (i = 0; i <= deg; ++i)
            alpha[i] *= fma(static_cast<T>(3.8), static_cast<T>(i), static_cast<T>(1));
        // Main loop
        for (i = 1; i <= itmax; ++i)
        for (j = 0; j < deg; ++j) {
            if (conv[j] == 0) {
                z = roots[j];
                r = fabs(z);
                // 
                if (r > 1.0)
                    rcheck_lag(poly, alpha, deg, b, c, z, r, conv[j], berr[j], cond[j]);
                else
                    check_lag(poly, alpha, deg, b, c, z, r, conv[j], berr[j], cond[j]);

                if (conv[j] == 0) {
                    modify_lag(deg, b, c, z, j, roots);
                    roots[j] -= c;
                } 
                else {
                    nz++;
                    if (nz == deg){
                        if (*std::min_element(conv.begin(), conv.end()) == 1)
                            return;
                        // Display a warning
                        std::cout << "Some root approximations did not converge or experienced overflow/underflow.\n";

                        // Compute backward error and condition number for roots that did not converge.
                        for (j = 0; j < deg; ++j){
                            if (conv[j] != 1) {
                                z = roots[j];
                                r = abs(z);

                                if (r > 1.0) {
                                    z = 1.0 / z;
                                    r = 1.0 / r;
                                    c = 0;
                                    b = poly[0];
                                    berr[j] = alpha[0];

                                    for (i = 1; i < deg + 1; ++i) {
                                        c = fma(z, c, b);
                                        b = fma(z, b, poly[i]);
                                        berr[j] = fma(r, berr[j], alpha[i]);
                                    }

                                    cond[j] = berr[j] / fabs(fms(complex<T>(deg, 0), b, z, c));
                                    berr[j] = fabs(b) / berr[j];
                                }
                                else {
                                    c = 0.0;
                                    b = poly[deg];
                                    berr[j] = alpha[deg];

                                    for (i = deg - 1; i >= 0; --i) {
                                        c = fma(z, c, b);
                                        b = fma(z, b, poly[i]);
                                        berr[j] = fma(r, berr[j], alpha[i]);
                                    }

                                    cond[j] = berr[j] / (r * fabs(c));
                                    berr[j] = fabs(b) / berr[j];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

};

};
#endif
