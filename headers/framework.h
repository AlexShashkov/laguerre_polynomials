// Александр, Дмитрий

#ifndef POLYGEN
#define POLYGEN

#define PR_NUMBERS_OF_ROOTS_EQUAL     0
#define PR_AT_LEAST_ONE_ROOT_LOST    -1
#define PR_AT_LEAST_ONE_ROOT_IS_FAKE -2
#define PR_2_INFINITE_ROOTS          -3

#include "ttmath/ttmath.h" // Bignum C++ library by Tomasz Sowa
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
        int crnt_idx=0;
        int i = 0;

        // current root id
        re = rnr(rng); // Generate the first root

        //  Generate clustered roots
        // std::cout << "\nGenerating clustered roots: " <<  N_clustered_roots-first4_clustered_roots;
        for (i = 0; i < N_clustered_roots; ++i) {
            im = rnc(rng); // Generate a small random number
            re = roots[i] = (re > root_mid_sweep) ? re - im : re + im; // Generate the next root close to the previous one
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

        /*
        std::cout << "\n GENERATED ROOTS ";
        for(auto &el : roots){
            std::cout<< el << ',';
        }
        */
        
        // Calculate resulting coefficients
        std::vector<ttmath::Big<exponent,mantissa>> big_coeffs, big_coeffs_new, big_roots;
        for (const fp_t num : coefficients) {
            ttmath::Big<exponent,mantissa> bignum(num);
            big_coeffs.push_back(bignum);
        }
        // roots = {-0.594628, -0.706948, -0.743337, -0.702158, -0.608897};
        for (const fp_t num : roots) {
            ttmath::Big<exponent,mantissa> bignum(num);
            big_roots.push_back(bignum);
        }
        big_coeffs_new = big_coeffs;

        for (i = 0; i < P; ++i){
          // std::cout << "\nROOT #" << i << "\n";
          for (int j = P-2; j >= 0; --j){
              // std::cout << "j:" << j << "\n";
              // std::cout << "-coeff[j+1]:" << -big_coeffs[j+1] << "; roots[i]:" << big_roots[i] << "; " << "coeff[j]: " << big_coeffs[j];
              // std::cout << "; fma: " << std::fma(-big_coeffs[j+1].ToDouble(), big_roots[i].ToDouble(), big_coeffs[j].ToDouble());
              big_coeffs_new[j] = ttmath::Big<exponent,mantissa>(std::fma(-big_coeffs[j+1].ToDouble(), big_roots[i].ToDouble(), big_coeffs[j].ToDouble()));
              // std::cout << "\n coeff_new[j]:" << big_coeffs_new[j] << "\n";
          }
          big_coeffs_new[P-1] -= big_roots[i];
          // std::cout << "\ncoeff_new[P-1]: " << big_coeffs_new[P-1];
          // coefficients = coefficients_new;
          big_coeffs = big_coeffs_new;
        }

        // Generate last complex roots
        // std::cout << "\nGenerating complex roots: " <<  N_pairs_of_complex_roots << " pairs";
        for (i = 0; i < N_pairs_of_complex_roots; ++i) {
            re=rnr(rng); while ((im=rnr(rng))==static_cast<fp_t>(0.0L)) {}
            auto c1=static_cast<fp_t>(-2.0L*re); // -2*re
            auto c2=static_cast<fp_t>(pr_product_difference(re, re, -im, im)); // re*re+im*im

            for (int j = P-2; j >= 0; --j){
              big_coeffs_new[j] = ttmath::Big<exponent,mantissa>(std::fma(big_coeffs[j].ToDouble(), c1,
                                                   std::fma(big_coeffs[j+1].ToDouble(), c2, big_coeffs[j+2].ToDouble())));
            }
            big_coeffs_new[P] *= ttmath::Big<exponent,mantissa>(c2) ; // first not null element
            big_coeffs = big_coeffs_new;
            // std::cout << "\nIM: " << im << "\n";
            // roots[P-i*2-1] = re; // In future we can add support of complex numbers 
            // roots[P-i*2-2] = re;
        }
        std::cout << "Bignum coeffs: \n";
        for (int i=0; i < P+1; ++i) {
            // Print out coeffs to show that it fails on float
            std::cout << big_coeffs[i] << " ";
            coefficients[i] = static_cast<fp_t>(big_coeffs[i].ToDouble());
        }
        std::cout << "\n";
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
int compare_roots(
unsigned N_roots_to_check,
unsigned N_roots_ground_truth,
std::vector<fp_t> &roots_to_check,
std::vector<fp_t> &roots_ground_truth,
fp_t &max_absolute_error,
fp_t &max_relative_error){
    int rv = (N_roots_to_check<N_roots_ground_truth) ? PR_AT_LEAST_ONE_ROOT_LOST :
      ( (N_roots_to_check>N_roots_ground_truth) ? PR_AT_LEAST_ONE_ROOT_IS_FAKE : PR_NUMBERS_OF_ROOTS_EQUAL );
    long double abs = std::numeric_limits<long double >::max();
    long double  rel = std::numeric_limits<long double >::max();
    auto size = roots_to_check.size();
    for(int j = 0;j<size; j++)
    for(int i = 0;i < size; i++){
        long double  absLoc = std::abs((long double)(roots_ground_truth[i])-(long double)(roots_to_check[(i + j) % size]));
        abs = std::min(absLoc,abs);
        rel = std::min(std::abs(
                (long double)(absLoc + std::numeric_limits<fp_t>::epsilon())/
                        (long double)(std::max(roots_to_check[(i + j) % size],roots_ground_truth[i]) + std::numeric_limits<fp_t>::epsilon())),rel);
    }
    max_absolute_error = abs;
    max_relative_error = rel;
    return rv;
}

#endif