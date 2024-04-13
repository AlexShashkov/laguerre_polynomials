// Александр, Дмитрий

#ifndef POLYGEN
#define POLYGEN

#define PR_NUMBERS_OF_ROOTS_EQUAL     0
#define PR_AT_LEAST_ONE_ROOT_LOST    -1
#define PR_AT_LEAST_ONE_ROOT_IS_FAKE -2
#define PR_2_INFINITE_ROOTS          -3

#include "ttmath/ttmath.h" // Bignum C++ library by Tomasz Sowa
#include "ExtendedFunctions.h"
#include <cmath>
#include <numbers> // std::numbers::pi_v<fp_t>, requires -std=c++20
#include <chrono> 
#include <random>
#include <cstdlib> 
#include <cassert>
#include <iostream>
#include <iomanip>
#include <complex>

/* computes (a*b - c*d) with precision not worse than 1.5*(unit of least precision) suggested in Claude-Pierre Jeannerod,
Nicolas Louvet, and Jean-Michel Muller, "Further Analysis of Kahan's Algorithm for the Accurate Computation of 2x2 Determinants".
Mathematics of Computation, Vol. 82, No. 284, Oct. 2013, pp. 2245-2264 */
template <typename fp_t> inline fp_t pr_product_difference(fp_t a, fp_t b, fp_t c, fp_t d)
{
  auto tmp = d * c;
  return fma(a, b, -tmp) + fma(-d, c, tmp);
}

template<typename T>
inline T min_value(T a, T b) {
    return (a < b) ? a : b;
}

template<typename T>
inline T max_value(T a, T b) {
    return (a > b) ? a : b;
}

/**
 * \brief Generates a test polynomial of any given degree P.
 * 
 * The polynomial is represented both in the form of roots, e.g. (x-roots[0])*(x-roots[1])*(quadratic polynomial with no real roots),
 * as well as by its coefficients, e.g. (coefficients[n]=1)*x^n + ... + coefficients[4]*x^4 + coefficients[3]*x^3 + coefficients[2]*x^2 + coefficients[1]*x + coefficients[0].
 * The highest-degree coefficient always equals 1.
 * 
 * \tparam fp_t The floating point type.
 * \tparam exponent The exponent of the floating point type for bignum.
 * \tparam mantissa The mantissa of the floating point type for bignum.
 * 
 * \param P Polynomial degree.
 * \param N_pairs_of_complex_roots How many pairs of complex conjugate roots to introduce.
 * \param N_clustered_roots How many clustered roots to introduce; all the clustered roots are real.
 * \param N_multiple_roots How many multiple roots to introduce; all multiple roots are real.
 * \param max_distance_between_clustered_roots Maximal distance between the closest of the clustered roots.
 * \param root_sweep_low Low boundary of real roots; imaginary parts of complex conjugate roots are in the same range.
 * \param root_sweep_high High boundary of real roots; imaginary parts of complex conjugate roots are in the same range.
 * \param roots Storage where to put the roots; size should exceed P-1.
 * \param coefficients Storage where to put the coefficients; size should exceed P.
 * 
 * \return The actual number of different real roots placed into the vector (roots) (complex roots are not placed there).
 * Negative return values may mean internal implementation error.
 */
template<typename fp_t, int exponent, int mantissa> 
int generate_polynomial(
unsigned P,
unsigned N_pairs_of_complex_roots,
unsigned N_clustered_roots,
unsigned N_multiple_roots,
fp_t max_distance_between_clustered_roots,
fp_t root_sweep_low, fp_t root_sweep_high,
std::vector<fp_t> &roots,
std::vector<fp_t> &coefficients)
{
int N_simple_roots=P-2*N_pairs_of_complex_roots-N_clustered_roots-N_multiple_roots;
assert(N_clustered_roots!=1); assert(N_multiple_roots!=1); assert(N_simple_roots>=0); assert(P>0);
assert(max_distance_between_clustered_roots>static_cast<fp_t>(0.0L));
assert(root_sweep_high-root_sweep_low>2*P*max_distance_between_clustered_roots);

coefficients[P]=static_cast<fp_t>(1.0L); // invariant
unsigned long long seed=std::chrono::system_clock::now().time_since_epoch().count()+std::rand(); // counts milliseconds
std::mt19937_64 rng(seed); // randomize seed from the clock
std::uniform_real_distribution<fp_t> rnr(root_sweep_low, root_sweep_high); 
std::uniform_real_distribution<fp_t> rnc(static_cast<fp_t>(0.0L), max_distance_between_clustered_roots); // uniform random data generator for root clusters

fp_t re, im, root_mid_sweep=root_sweep_low+0.5*(root_sweep_high-root_sweep_low);

switch (P)
  {
  case 0:
    coefficients[0]=rnr(rng); return 0;
  case 1:
    coefficients[0]=-(roots[0]=rnr(rng)); return 1;
  default:
        unsigned crnt_idx=0;
        unsigned i = 0;

        // current root id
        re = rnr(rng); // Generate the first root

        //  Generate clustered roots
        // std::cout << "\nGenerating clustered roots: " <<  N_clustered_roots-first4_clustered_roots;
        for (i = 0; i < N_clustered_roots; ++i) {
            im = rnc(rng); // Generate a small random number
            roots[i] = (re > root_mid_sweep) ? re - im : re + im; // Generate the next root close to the previous one
            re = roots[i];
        }
        crnt_idx = i;

        // Generate multiple roots
        // std::cout << "\nGenerating multiple roots: " <<  N_multiple_roots;
        re = N_clustered_roots ? rnr(rng): re;
        for(; i < crnt_idx + N_multiple_roots; ++i){
          roots[i] = re;
        }
        crnt_idx = i;

        // Ggenerate simple roots
        // std::cout << "\nGenerating simple roots: " <<  N_simple_roots;
        // std::cout << "\n FROM " << crnt_idx << " TO " << crnt_idx + N_simple_roots;
        for(; i < crnt_idx + N_simple_roots; ++i){
          roots[i] = rnr(rng);
        }
        
        // Calculate resulting coefficients
        std::vector<ttmath::Big<exponent,mantissa>> big_coeffs, big_coeffs_new, big_roots;
        for (const fp_t num : coefficients) {
            big_coeffs.push_back(ttmath::Big<exponent,mantissa>(std::to_string(num)));
        }
        for (const auto num : roots) {
            big_roots.push_back(ttmath::Big<exponent,mantissa>(std::to_string(num)));
        }
        big_coeffs_new = big_coeffs;

        for (i = 0; i < P; ++i){
          // std::cout << "\nROOT #" << i << "\n";
          for (int j = P-2; j >= 0; --j){
              // big_coeffs_new[j] = ttmath::Big<exponent,mantissa>(std::fma(-big_coeffs[j+1].ToDouble(), big_roots[i].ToDouble(), big_coeffs[j].ToDouble()));
              // On some extreme cases even double type will overflow (on degree >= 200), as there is no direct casting from bignum to long double (the only way is to use string convertation)
              // We will allow bignum library to calculate a*b+c with high precision
              big_coeffs_new[j] = -big_coeffs[j+1]*big_roots[i] + big_coeffs[j];
          }
          big_coeffs_new[P-1] -= big_roots[i];
          big_coeffs = big_coeffs_new;
        }
        
        // complex roots generation
        for (i = 0; i < N_pairs_of_complex_roots; ++i) {
          re = rnr(rng); while ((im = rnr(rng)) == static_cast<fp_t>(0.0L)) {}
          std::cout << "(x-(" << re << " + " << im << "i))" << "\n";

          // TODO: Change to convertion from string
          big_coeffs_new = {
              ttmath::Big<exponent, mantissa>(static_cast<double>(re * re + im * im)),
              ttmath::Big<exponent, mantissa>(static_cast<double>(-2 * re)),
              ttmath::Big<exponent, mantissa>(1)
          };

          bool isNotFullyComplex = N_pairs_of_complex_roots*2 != P;
          big_coeffs = isNotFullyComplex || i != 0 ? 
            Laguerre::multiply(big_coeffs, big_coeffs_new, i == 0 && isNotFullyComplex ? N_pairs_of_complex_roots*2 : 0) 
            : big_coeffs_new;

          roots[P - i * 2 - 1] = re; 
          roots[P - i * 2 - 2] = re;
      }

        for (int i=0; i < P+1; ++i) {
            // Print out coeffs to show that it fails on float or double in extreme cases
            // std::cout << big_coeffs[i] << " ";
            // bignum to string, and then casting back to fp_t from long double value (on high degrees coeffs exist that are greater than 1e+308)
            coefficients[i] = static_cast<fp_t>(std::strtold(big_coeffs[i].ToString().c_str(), nullptr));

            // coefficients[i] = static_cast<fp_t>(big_coeffs[i].ToDouble());
        }
        return P - (N_pairs_of_complex_roots*2);
  }
return -1; // unreachable, means a flaw in control here
}

/**
 * \brief Compares two vectors of roots.
 * 
 * Root orderings play no role. For each entry in (roots_ground_truth), the closest entry in (roots_to_check) is found and corresponding distance found.
 * Among such distances, the largest will be stored to (max_deviation).
 * 
 * \tparam fp_t The floating point type.
 * 
 * \param N_roots_to_check Number of roots in (roots_to_check).
 * \param N_roots_ground_truth Number of roots in (roots_ground_truth).
 * \param roots_to_check One should take into account only first (N_roots_to_check) roots here.
 * \param roots_ground_truth One should take into account only first (N_roots_ground_truth) roots here.
 * \param max_absolute_error Here the greatest among the smallest deviations of the roots in (roots_to_check) and (roots_ground_truth) will be placed.
 * \param max_relative_error The maximum relative error.
 * 
 * \return Negative return values may mean internal implementation error.
 */
template<typename fp_t>
int compare_roots(unsigned N_roots_to_check,
                  unsigned N_roots_ground_truth,
                  std::vector<fp_t> &roots_to_check,
                  std::vector<fp_t> &roots_ground_truth,
                  fp_t &max_absolute_error,
                  fp_t &max_relative_error) {
    int rv = (N_roots_to_check < N_roots_ground_truth) ? PR_AT_LEAST_ONE_ROOT_LOST :
             ((N_roots_to_check > N_roots_ground_truth) ? PR_AT_LEAST_ONE_ROOT_IS_FAKE : PR_NUMBERS_OF_ROOTS_EQUAL);

    long double global_ae = std::numeric_limits<long double>::min();
    long double global_re = std::numeric_limits<long double>::min();

    long double abs;
    long double re;
    long double absLoc;
    long double reLoc;

    for (size_t i = 0; i < N_roots_ground_truth; ++i) {
        abs = std::numeric_limits<long double>::max();
        for (size_t j = 0; j < N_roots_to_check; ++j) {
          // std::cout << "Comparing GT " << roots_ground_truth[i] << " against " << roots_to_check[j] << "\n";
            absLoc = std::abs(static_cast<long double>(roots_ground_truth[i]) -
                                    static_cast<long double>(roots_to_check[j]));
            reLoc = std::abs((absLoc + std::numeric_limits<long double>::epsilon()) /
                                               (static_cast<long double>(std::max(roots_ground_truth[i], roots_to_check[j])) + std::numeric_limits<long double>::epsilon()));
            // std::cout << "Absloc: " << absLoc << ", abs: " << abs << " => abs=min(abs, absLoc)";
            abs = std::min(absLoc, abs);
            re = std::min(reLoc, re);
            // std::cout << "=" << abs << "\n";
        }

        // std::cout << "Global abs=" << global_ae << ", abs=" << abs << " => ";
        global_ae = std::max(abs, global_ae);
        // std::cout << "global = " << global_ae << "\n";
        global_re = std::max(re, global_re);
    }

    max_absolute_error = static_cast<fp_t>(global_ae);
    max_relative_error = static_cast<fp_t>(global_re);

    return rv;
}


#endif