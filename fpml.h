#include <complex>
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>

const double eps = std::numeric_limits<double>::epsilon();
const double big = std::numeric_limits<double>::max();
const double small = std::numeric_limits<double>::min();

/**
 * @brief Calculate the determinant between two points.
 *
 * @param h Vector of integers representing points.
 * @param a Vector of doubles representing points' coordinates.
 * @param c Index for the current point.
 * @param i Index for the next point.
 * @return double Determinant between two points.
 */
double cross(std::vector<int>& h, std::vector<double>& a, int c, int i) {
    return (a[i] - a[h[c - 1]]) * (h[c] - h[c - 1]) - (a[h[c]] - a[h[c - 1]]) * (i - h[c - 1]);
}

/**
 * @brief Calculate the convex hull of a set of points.
 *
 * @param n Number of points.
 * @param a Vector of doubles representing points' coordinates.
 * @param h Vector to store the convex hull.
 * @param c Number of points in the convex hull.
 */
void conv_hull(int n, std::vector<double>& a, std::vector<int>& h, int& c) {
    c = 0;

    for (int i = n - 1; i >= 0; i--) {
        while (c >= 2 && cross(h, a, c, i) < eps)
            c--;
        c++;
        h[c] = i;
    }
}

/**
 * @brief Estimate roots of a polynomial.
 *
 * @param alpha Vector of coefficients.
 * @param deg Degree of the polynomial.
 * @param roots Vector to store the estimated roots.
 * @param conv Vector to store convergence status of each root.
 * @param nz Number of zeros.
 */
void estimates(std::vector<double>& alpha, int deg, std::vector<std::complex<double>>& roots, std::vector<int>& conv, int& nz) {
    int c, i, j, k, nzeros;
    double a1, a2, ang, r, th;
    std::vector<int> h(deg + 1);
    std::vector<double> a(deg + 1);
    const double pi2 = 6.2831853071795865;
    const double sigma = 0.7;

    // Log of absolute value of coefficients
    for (i = 0; i < deg + 1; i++) {
        if (alpha[i] > 0) {
            a[i] = log(alpha[i]);
        } else {
            a[i] = -1E+30;
        }
    }
    conv_hull(deg + 1, a, h, c);
    k = 0;
    th = pi2 / deg;

    // Initial Estimates
    for (i = c - 1; i >= 0; i--) {
        nzeros = h[i] - h[i + 1];
        a1 = pow(alpha[h[i + 1]], 1.0 / nzeros);
        a2 = pow(alpha[h[i]], 1.0 / nzeros);

        if (a1 <= a2 * small) {
            // r is too small
            r = 0.0;
            nz += nzeros;
            for (j = k; j < k + nzeros; j++)
                conv[j] = -1;
            for (j = k; j < k + nzeros; j++)
                roots[j] = {0.0, 0.0};
        } else if (a1 >= a2 * big) {
            // r is too big
            r = big;
            nz += nzeros;
            for (j = k; j < k + nzeros; j++)
                conv[j] = -1;
            ang = pi2 / nzeros;
            for (j = 0; j < nzeros; j++)
                roots[k + j] = r * std::complex<double>(cos(ang * j + th * h[i] + sigma), sin(ang * j + th * h[i] + sigma));
        } else {
            // r is just right
            r = a1 / a2;
            ang = pi2 / nzeros;
            for (j = 0; j < nzeros; j++)
                roots[k + j] = r * std::complex<double>(cos(ang * j + th * h[i] + sigma), sin(ang * j + th * h[i] + sigma));
        }
        k += nzeros;
    }
}

/**
 * @brief Check if a complex number contains NaN or Inf values.
 *
 * @param a Complex number to check.
 * @return bool True if NaN or Inf is found, false otherwise.
 */
bool check_nan_inf(std::complex<double> a) {
    double re_a = real(a);
    double im_a = imag(a);
    return std::isnan(re_a) || std::isnan(im_a) || (abs(re_a) > big) || (abs(im_a) > big);
}

/**
 * @brief Perform Laguerre correction, compute backward error, and condition.
 *
 * @param p Vector of polynomial coefficients.
 * @param alpha Vector of coefficients.
 * @param deg Degree of the polynomial.
 * @param b Output complex number representing root.
 * @param c Output complex number representing correction.
 * @param z Complex number used for correction.
 * @param r Value used for correction.
 * @param conv Output value indicating convergence status.
 * @param berr Output backward error.
 * @param cond Output condition number.
 */
void rcheck_lag(std::vector<std::complex<double>>& p, std::vector<double>& alpha, int deg, std::complex<double>& b, std::complex<double>& c, std::complex<double> z, double r, int& conv, double& berr, double& cond) {
    int k;
    double rr;
    std::complex<double> a, zz;

    // Evaluate polynomial and derivatives
    zz = 1.0 / z;
    rr = 1.0 / r;
    a = p[0];
    b = 0.0;
    c = 0.0;
    berr = alpha[0];

    for (k = 1; k < deg + 1; k++) {
        c = zz * c + b;
        b = zz * b + a;
        a = zz * a + p[k];
        berr = rr * berr + alpha[k];
    }

    // Laguerre correction/ backward error and condition
    if (abs(a) > eps * berr) {
        b = b / a;
        c = 2.0 * (c / a);
        c = pow(zz, 2) * ((double)deg - 2.0 * zz * b + pow(zz, 2) * (pow(b, 2) - c));
        b = zz * ((double)deg - zz * b);

        if (check_nan_inf(b) || check_nan_inf(c))
            conv = -1;
    } else {
        cond = berr / abs((double)deg * a - zz * b);
        berr = abs(a) / berr;
        conv = 1;
    }
}

/**
 * @brief Perform Laguerre correction, compute backward error, and condition.
 *
 * @param p Vector of polynomial coefficients.
 * @param alpha Vector of coefficients.
 * @param deg Degree of the polynomial.
 * @param b Output complex number representing root.
 * @param c Output complex number representing correction.
 * @param z Complex number used for correction.
 * @param r Value used for correction.
 * @param conv Output value indicating convergence status.
 * @param berr Output backward error.
 * @param cond Output condition number.
 */
void check_lag(std::vector<std::complex<double>>& p, std::vector<double>& alpha, int deg, std::complex<double>& b, std::complex<double>& c, std::complex<double> z, double r, int& conv, double& berr, double& cond) {
    int k;
    std::complex<double> a;

    // Evaluate polynomial and derivatives
    a = p[deg];
    b = 0.0;
    c = 0.0;
    berr = alpha[deg];

    for (k = deg - 1; k >= 0; k--) {
        c = z * c + b;
        b = z * b + a;
        a = z * a + p[k];
        berr = r * berr + alpha[k];
    }

    // Laguerre correction/ backward error and condition
    if (abs(a) > eps * berr) {
        b = b / a;
        c = pow(b, 2) - 2.0 * (c / a);

        if (check_nan_inf(b) || check_nan_inf(c))
            conv = -1;
    } else {
        cond = berr / (r * abs(b));
        berr = abs(a) / berr;
        conv = 1;
    }
}

/**
 * @brief Modify Laguerre correction based on Aberth's method.
 *
 * @param deg Degree of the polynomial.
 * @param b Complex number representing root.
 * @param c Complex number representing correction.
 * @param z Complex number used for correction.
 * @param j Index of the root to modify.
 * @param roots Vector of estimated roots.
 */
void modify_lag(int deg, std::complex<double>& b, std::complex<double>& c, std::complex<double> z, int j, std::vector<std::complex<double>>& roots) {
    std::complex<double> t;

    // Aberth correction terms
    for (int k = 0; k < j - 1; k++) {
        t = 1.0 / (z - roots[k]);
        b = b - t;
        c = c - pow(t, 2);
    }

    for (int k = j + 1; k <= deg; k++) {
        t = 1.0 / (z - roots[k]);
        b = b - t;
        c = c - pow(t, 2);
    }

    // Laguerre correction
    t = sqrt((double)(deg - 1) * ((double)deg * c - pow(b, 2)));
    c = b + t;
    b = b - t;

    if (abs(b) > abs(c)) {
        c = (double)deg / b;
    } else {
        c = (double)deg / c;
    }
}

/**
 * @brief Find the roots of a polynomial using Laguerre's method.
 *
 * @param poly Vector of polynomial coefficients.
 * @param deg Degree of the polynomial.
 * @param roots Vector to store the roots.
 * @param berr Vector to store the backward errors.
 * @param cond Vector to store the condition numbers.
 * @param conv Vector to store convergence status of each root.
 * @param itmax Maximum number of iterations.
 */
void laguerre(std::vector<std::complex<double>>& poly, int deg, std::vector<std::complex<double>>& roots, std::vector<double>& berr, std::vector<double>& cond, std::vector<int>& conv, int itmax) {
    int i, j, nz;
    double r;
    std::vector<double> alpha(deg + 1);
    std::complex<double> b = 1.0, c = 1.0, z = 1.0;

    // Precheck
    for (i = 0; i <= deg; i++)
        alpha[i] = abs(poly[i]);

    if (alpha[deg] < small) {
        std::cout << "Warning: leading coefficient too small." << std::endl;
        return;
    } else if (deg == 1) {
        roots[0] = -poly[0] / poly[1];
        conv[0] = 1;
        berr[0] = 0.0;
        cond[0] = (alpha[0] + alpha[1] * abs(roots[0])) / (abs(roots[0]) * alpha[1]);
        return;
    } else if (deg == 2) {
        b = -poly[1] / (2.0 * poly[2]);
        c = sqrt(pow(poly[1], 2) - 4.0 * poly[2] * poly[0]) / (2.0 * poly[2]);
        roots[0] = b - c;
        roots[1] = b + c;
        conv[0] = 1;
        conv[1] = 1;
        berr[0] = 0.0;
        berr[1] = 0.0;
        cond[0] = (alpha[0] + alpha[1] * abs(roots[0]) + alpha[2] * pow(abs(roots[0]), 2)) / (abs(roots[0]) * abs(poly[1] + 2.0 * poly[2] * roots[0]));
        cond[1] = (alpha[0] + alpha[1] * abs(roots[1]) + alpha[2] * pow(abs(roots[1]), 2)) / (abs(roots[1]) * abs(poly[1] + 2.0 * poly[2] * roots[1]));
        return;
    }

    // Initial estimates
    conv.assign(deg, 0);
    nz = 0;
    estimates(alpha, deg, roots, conv, nz);

    // Main loop
    for (i = 1; i <= deg + 1; i++)
        alpha[i - 1] = alpha[i - 1] * (3.8 * (i - 1) + 1);

    for (i = 1; i <= itmax; i++) {
        for (j = 0; j < deg; j++) {
            if (conv[j] == 0) {
                z = roots[j];
                r = abs(z);
                if (r > 1.0) {
                    rcheck_lag(poly, alpha, deg, b, c, z, r, conv[j], berr[j], cond[j]);
                } else {
                    check_lag(poly, alpha, deg, b, c, z, r, conv[j], berr[j], cond[j]);
                }

                if (conv[j] == 0) {
                    modify_lag(deg, b, c, z, j, roots);
                    roots[j] = roots[j] - c;
                } else {
                    nz++;
                    if (nz == deg)
                        goto label10;
                }
            }
        }
    }

    // Final check
    label10:
    if (*std::min_element(conv.begin(), conv.end()) == 1) {
        return;
    } else {
        // Display a warning
        std::cout << "Some root approximations did not converge or experienced overflow/underflow." << std::endl;

        // Compute backward error and condition number for roots that did not converge.
        // Note that this may produce overflow/underflow.
        for (j = 0; j < deg; j++) {
            if (conv[j] != 1) {
                z = roots[j];
                r = abs(z);

                if (r > 1.0) {
                    z = 1.0 / z;
                    r = 1.0 / r;
                    c = 0.0;
                    b = poly[0];
                    berr[j] = alpha[0];

                    for (i = 1; i < deg + 1; i++) {
                        c = z * c + b;
                        b = z * b + poly[i];
                        berr[j] = r * berr[j] + alpha[i];
                    }

                    cond[j] = berr[j] / abs((double)deg * b - z * c);
                    berr[j] = abs(b) / berr[j];
                } else {
                    c = 0.0;
                    b = poly[deg];
                    berr[j] = alpha[deg];

                    for (i = deg - 1; i >= 0; i--) {
                        c = z * c + b;
                        b = z * b + poly[i];
                        berr[j] = r * berr[j] + alpha[i];
                    }

                    cond[j] = berr[j] / (r * abs(c));
                    berr[j] = abs(b) / berr[j];
                }
            }
        }
    }
}