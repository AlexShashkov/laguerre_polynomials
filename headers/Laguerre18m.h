// Александр, Дмитрий
// An effective implementation of a modified Laguerre method for the roots of a polynomial

#ifndef LAGUERRE18
#define LAGUERRE18

#include <vector>
#include <limits>
#include <numbers> // std::numbers::pi_v<fp_t>, requires -std=c++20
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

    // for backward error
    static constexpr T e_i   = fma(static_cast<T>(2), static_cast<T>(sqrtl(2L)), 1);

    /**
     * \brief Calculate the cross product between two points.
         * \param h Vector of integers representing points.
         * \param a Vector of Ts representing points' coordinates.
         * \param c Index for the current point.
         * \param i Index for the next point.
         * \return T cross product between two points.
     */
    inline T cross(std::vector<int>& h, std::vector<T>& a, int c, int i) {
        // determinant
        int h_c = h[c];
        int h_c1 = h[c-1];
        T a_hc = a[h_c1];
        return fms(a[i]-a_hc, static_cast<T>(h_c-h_c1), a[h_c]-a_hc, static_cast<T>(i-h_c1));
    }

    /**
     * \brief Calculate the convex hull of a set of points. 
         * \param n Number of points.
         * \param a Vector of templated type T representing points coordinates.
         * \param h Vector to store the convex hull.
         * \param c Number of points in the convex hull.
     */
    inline void conv_hull(int n, std::vector<T>& a, std::vector<int>& h, int& c) {
        c = 0;

        for (int i = n; i >= 0; --i) {
            while (c >= 2 && cross(h, a, c, i) < eps)
                --c;
            // std::cout << "i=" << i << ", c=" << c << "\n"; 
            h[c] = i;
            ++c;
            // ++c;
        }
    }

    /**
     * \brief Estimate roots of a polynomial.
         * \param alpha Polynomial object.
         * \param deg Degree of the polynomial.
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
        for (i = 0; i <= deg; ++i){
            a[i] = log(alpha[i]);
            a[i] = std::isfinite(a[i]) ? a[i] : -big;
        }
        conv_hull(deg, a, h, c);
        k = 0;
        th = pi2/static_cast<T>(deg);

        // Initial Estimates
        for (i = c - 2; i >= 0; --i) {
            nzeros = h[i] - h[i + 1];
            if(nzeros == static_cast<T>(0)) throw std::invalid_argument("Convex hull values cant be zero.");
            if(alpha[h[i]] == static_cast<T>(0)) throw std::invalid_argument("Alpha value cant be zero");
            T one_div_nzeros = static_cast<T>(1.0 / nzeros);
            a1 = pow(alpha[h[i + 1]], one_div_nzeros);
            a2 = pow(alpha[h[i]], one_div_nzeros);

            if (a1 <= a2 * small){
                // r is too small
                r = 0.0;
                nz += nzeros;
                for (j = k; j < k + nzeros; ++j){
                    conv[j] = -1;
                    roots[j] = {0.0, 0.0};
                }
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
                    // r is ok :big-finger-emoji:
                    if (a2 == static_cast<T>(0)) throw std::invalid_argument("Both correction values equal to zero.");
                    r = a1 / a2;
                }
                ang = pi2*one_div_nzeros;
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
        if (z == complex<T>(0, 0)) throw std::invalid_argument("correction absolute value cannot be zero!");
        if (r == static_cast<T>(0)) throw std::invalid_argument("correction value cannot be zero!");
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

        for (k = 1; k <= deg; ++k) {
            c = fma(zz, c, b);
            b = fma(zz, b, a);
            a = fma(zz, a, p[k]);
            berr = fma(rr, berr, alpha[k]);
        }
        if(a == static_cast<T>(0)) throw std::invalid_argument("Alpha value cant be zero");

        // Laguerre correction/ backward error and condition
        bool condition = abs(a) > (eps * berr);
        // std::cout << "berr is " << berr << "\n";
        b = condition ? b/a : b;
        // https://www.wolframalpha.com/input?i=z%5E2*%28-2*z*b%2Bd%29+-+%28z%5E2*z%5E2%29*%282*c%2Fa-b%5E2%29
        c = condition ? fms(zz2, fma(-static_cast<T>(2)*zz, b, static_cast<T>(deg)), zz2*zz2, fms(complex<T>(2.0, 0), c / a, b, b)) : c;
        b = condition ? fms(zz, complex<T>(deg, 0) , zz2, b) : b;

        // cond = (!condition)*(berr / abs(fms(complex<T>(deg, 0),  a, zz, b))) + (condition)*cond;
        // TODO: CHECK IF ITS CORRECT
        cond = fma(!condition, berr / std::abs(fms(std::complex<T>(deg, 0), a, zz, b)), condition * cond);
        berr = (!condition)*(abs(a) / berr) + (condition)*berr;

        conv = condition ? 
            (complexnotfinite(b, big) || complexnotfinite(c, big) ? -1 : conv)
            : 1;
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
        // b^2 + (-{2,0}*c/a)+{2,0}*c/a+(-{2,0}*c/a)
        c = condition ? fms(b, b, std::complex<T>(2, 0), c/a) : c;

        cond = (!condition)*(berr / (r * abs(b))) + (condition)*cond;
        berr = (!condition)*(abs(a) / berr) + (condition)*berr;

        conv = condition ? 
            (complexnotfinite(b, big) || complexnotfinite(c, big) ? -1 : conv) 
            : 1;
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

        // std::cout << "Entered modify lag c=" << c << "\n";
        // Aberth correction terms
        for (int k = 0; k < j; ++k) {
            t = static_cast<T>(1.0) / (z - roots[k]);
            b -= t;
            c = fma(-t, t, c);
            // std::cout << "Aberth k-1 c=" << c << "\n";
        }
        // k != j
        for (int k = j + 1; k < deg; ++k) {
            t = static_cast<T>(1.0) / (z - roots[k]);
            b -= t;
            c = fma(-t, t, c);
            // std::cout << "Aberth k+1 c=" << c << "\n";
        }

        // Laguerre correction
        complex<T> t_arg = fms(complex<T>(static_cast<T>(deg - 1), 0), (static_cast<T>(deg) * c),  complex<T>(static_cast<T>(deg - 1), 0), b*b);
        if(complexnotfinite(t_arg, big)) throw std::invalid_argument("Sqrt of not finite number is undefined");
        t = sqrt(t_arg);
        c = b + t;
        // std::cout << "c = b+t c=" << c << "\n";
        b -= t;
        // std::complex<T> degc(deg, 0);
        c = abs(b) > abs(c) ? static_cast<T>(deg) / b : static_cast<T>(deg) / c;
        // std::cout << "Exiting c=" << c << "\n";
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
        int deg = roots.size();
        std::vector<T> berr(deg); // Vector to store the backward errors.
        std::vector<T> cond(deg); // Vector to store the condition numbers.
        int i, j, nz;
        T r;
        std::vector<T> alpha(deg + 1);
        std::complex<T> b, c, z;

        // Precheck
        for (i = 0; i <= deg; ++i)
            alpha[i] = fabs(poly[i]);

        if (alpha[deg] < small) {
            std::cout << "Warning: leading coefficient is too small." << std::endl;
            return;
        } 

        // Initial estimates
        conv.assign(deg, 0);
        nz = 0;
        estimates(alpha, deg, roots, conv, nz);

        // Motivated by the backward error analysis of the Ruffini-Horner rule in (9), we let
        // alpha = (2√2 + 1)i + 1. Page 9 in paper
        for (i = 0; i <= deg; ++i)
            alpha[i] *= fma(static_cast<T>(e_i), static_cast<T>(i), static_cast<T>(1));
        // Main loop
        
        for (i = 1; i <= itmax; ++i){
        if (nz == deg) break;
        for (j = 0; j < deg; ++j) {
            if (conv[j] == 0) {
                z = roots[j];
                r = abs(z);
                // 
                // std::cout << "\nr is " << r << "\n";
                // std::cout << "z is " << z << "\n";
                if (r > 1.0)
                    rcheck_lag(poly, alpha, deg, b, c, z, r, conv[j], berr[j], cond[j]);
                else
                    check_lag(poly, alpha, deg, b, c, z, r, conv[j], berr[j], cond[j]);

                if (conv[j] == 0) {
                    modify_lag(deg, b, c, z, j, roots);
                    roots[j] -= c;
                    // std::cout << "ROOT #" << j << "=" << roots[j] << "\n";
                } 
                else {
                    ++nz;
                    if (nz == deg){
                        if (*std::min_element(conv.begin(), conv.end()) == 1)
                            break;
                        // Display a warning
                        std::cout << "Some root approximations did not converge or experienced overflow/underflow.\n";

                        // Compute backward error and condition number for roots that did not converge.
                        for (j = 0; j < deg; ++j){
                            if (conv[j] != 1) {
                                z = roots[j];
                                r = abs(z);

                                if (r > 1.0) {
                                    z = static_cast<T>(1.0) / z;
                                    r = static_cast<T>(1.0) / r;
                                    c = 0;
                                    b = poly[0];
                                    berr[j] = alpha[0];

                                    for (i = 1; i < deg + 1; ++i) {
                                        c = fma(z, c, b);
                                        b = fma(z, b, poly[i]);
                                        berr[j] = fma(r, berr[j], alpha[i]);
                                    }

                                    cond[j] = berr[j] / abs(fms(complex<T>(deg, 0), b, z, c));
                                    berr[j] = abs(b) / berr[j];
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

                                    cond[j] = berr[j] / (r * abs(c));
                                    berr[j] = abs(b) / berr[j];
                                }
                            }
                        }
                        break;
                    }
                }
            }
        }
        }
        int nanpos = -1;
        for(i=0; i<deg; ++i){
            if(complexnotfinite(roots[i], big)){
                std::cout << "One of the roots had overflow at position " << i << "\n";
            }
        }
        return;
    }

};

};
#endif