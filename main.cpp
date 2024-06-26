// Александр Шашков cracecrunch@gmail.com
// Дмитрий Балашов dimabalash0v@yandex.ru

#include "headers/framework.h"
#include "headers/Polynomial.h"
#include "headers/Laguerre.h"
#include "headers/Laguerre13m.h"
#include "headers/Laguerre18m.h"

#include <fstream>

#define PR_NUMBERS_OF_ROOTS_EQUAL     0
#define PR_AT_LEAST_ONE_ROOT_LOST    -1
#define PR_AT_LEAST_ONE_ROOT_IS_FAKE -2
#define PR_2_INFINITE_ROOTS          -3


#ifndef SOLVER
    #define SOLVER 2
    // 0 - Original
    // 1 - 2013 mod
    // 2 - 2018 mod
#endif
#ifndef NUMBER
    #define NUMBER double
#endif


// Exponent and mantissa for root generation
// You can check recommended precision for ttmath on https://www.ttmath.org/online_calculator/
// lets use big precision
// FOR 64-bit 
// mantissa = 2048 bits / 64 bits per word = 32 words
// exponent = 256 bits / 64 bits per word = 4 words
// FOR 32-bit
// mantissa = 2048 bits / 32 bits per word = 64 words
// exponent = 256 bits / 32 bits per word = 8 words
#ifndef EXPONENT
    #define EXPONENT 4 // Exponent for bignum notation
#endif
#ifndef MANTISSA
    #define MANTISSA 32 // Mantissa for bignum notation
#endif


#ifndef DEGREE
    #define DEGREE 100 // polynomial degree
#endif
#ifndef N_TESTS
    #define N_TESTS 10000 // count of tests 
#endif
// Generator params. If all N_*_ROOTS are zero, then generator will create only simple roots
#ifndef N_PAIRS_OF_COMPLEX_ROOTS
    #define N_PAIRS_OF_COMPLEX_ROOTS 0
#endif 
#ifndef N_CLUSTERED_ROOTS
    #define N_CLUSTERED_ROOTS 0
#endif
#ifndef N_MULTIPLE_ROOTS
    #define N_MULTIPLE_ROOTS 0
#endif
// Used only for clustered roots
#ifndef MAX_DISTANCE_BETWEEN_CLUSTERED
    #define MAX_DISTANCE_BETWEEN_CLUSTERED 1e-5
#endif
#ifndef ROOT_SWEEP_LOW
    #define ROOT_SWEEP_LOW -1
#endif
#ifndef ROOT_SWEEP_HIGH
    #define ROOT_SWEEP_HIGH 1
#endif
#ifndef TRIES
    #define TRIES 80 // For how many iterations Laguerre should try to solve polynomial
#endif

using std::cout;
using std::vector;
using std::complex;
using std::for_each;
using std::fill;
using std::exception;
using std::invalid_argument;
using std::endl;
using std::cerr;

int main(){
    using Laguerre::Polynomial;
    
    // Laguerre solvers
    #if SOLVER == 0
        cout << "Using original version of Laguerre\n";
        Laguerre::Original<NUMBER>* solver = new Laguerre::Original<NUMBER>();
    #elif SOLVER == 1
        cout << "Using modified version of Laguerre, 2013\n";
        Laguerre::ModifiedLaguerre13<NUMBER>* solver = new Laguerre::ModifiedLaguerre13<NUMBER>();
    #elif SOLVER == 2
        cout << "Using modified version of Laguerre, 2018\n";
        Laguerre::ModifiedLaguerre18<NUMBER>* solver = new Laguerre::ModifiedLaguerre18<NUMBER>();
    #endif

    // Generator stuff
    int rv,                                                                  // status of comparing roots
        N_roots_found_this_test, N_roots_gt_this_test,                       // amount of found roots in each test & gt roots
        N_true_roots_lost=0, N_fake_roots_added=0,                           // total counters of {lost, fake} roots aover all tests
        N_roots_gt_all_tests_ae_worst=0, N_roots_found_all_tests_ae_worst=0, // {ground truth, found} NUMBER of roots of the worst polynomial for absolute error
        N_roots_gt_all_tests_re_worst=0, N_roots_found_all_tests_re_worst=0; // {ground truth, found} NUMBER of roots of the worst polynomial for relative error

    vector<NUMBER> roots(DEGREE, 0.0), roots_found_this_test(DEGREE), a(DEGREE+1, 0.0),                                       // current polynomial
        roots_gt_all_tests_ae_worst(DEGREE), roots_found_all_tests_ae_worst(DEGREE), coefficients_all_tests_ae_worst(DEGREE+1),    // {gt roots, found roots, coefficients} of the worst polynomial for absolute error
        roots_gt_all_tests_re_worst(DEGREE), roots_found_all_tests_re_worst(DEGREE), coefficients_all_tests_re_worst(DEGREE+1);    // {gt roots, found roots, coefficients} of the worst polynomial for relative error

    NUMBER ae_this_test_worst, ae_all_tests_worst=static_cast<NUMBER>(-1.0L),   // absolute error: an impossible value that will be updated for sure
           re_this_test_worst, re_all_tests_worst=static_cast<NUMBER>(-1.0L);   // relative error: an impossible value that will be updated for sure
    
    Polynomial<NUMBER> pol;
    try{
        pol.setSolver(solver);
        for(int i=0; i < N_TESTS; ++i){
            roots_found_this_test.clear();
            N_roots_gt_this_test = generate_polynomial<NUMBER, EXPONENT, MANTISSA>(DEGREE, N_PAIRS_OF_COMPLEX_ROOTS, N_CLUSTERED_ROOTS,
                N_MULTIPLE_ROOTS, MAX_DISTANCE_BETWEEN_CLUSTERED, ROOT_SWEEP_LOW, ROOT_SWEEP_HIGH, roots, a);

            cout << "GENERATED POLY & ROOTS:\n";
            Laguerre::printVec(a);
            Laguerre::printVec(roots);
            pol.setCoeffs(a);

            vector<complex<NUMBER>> solved_roots(DEGREE);
            vector<int> conv(DEGREE);

            pol.solve(solved_roots, conv, TRIES);

            cout << "found roots: ";
            for_each(solved_roots.begin(), solved_roots.end(), [&roots_found_this_test](complex<NUMBER> x){ // solved_roots = b_roots
                cout << x.real() << "; ";
                roots_found_this_test.push_back(x.real());
            });

            N_roots_found_this_test = roots_found_this_test.size();

            rv=compare_roots<NUMBER>(N_roots_found_this_test, N_roots_gt_this_test, roots_found_this_test, roots,
                                                                                    ae_this_test_worst, re_this_test_worst);


            cout << endl << "P=" << DEGREE << ", test No " << i << " out of " << N_TESTS << endl
                << "comparison return value=" << (rv==PR_NUMBERS_OF_ROOTS_EQUAL ? "PR_NUMBERS_OF_ROOTS_EQUAL" :
                (rv==PR_AT_LEAST_ONE_ROOT_LOST ? "PR_AT_LEAST_ONE_ROOT_LOST" : 
                (rv==PR_AT_LEAST_ONE_ROOT_IS_FAKE ? "PR_AT_LEAST_ONE_ROOT_IS_FAKE" : "UNKNOWN VALUE")))
                << ", N_roots_gt_this_test=" << N_roots_gt_this_test << ", N_roots_found_this_test=" << N_roots_found_this_test
                << ", N_roots_verified_this_test=" << DEGREE
                << "ae_this_test_worst=" << ae_this_test_worst
                << ", re_this_test_worst=" << re_this_test_worst << " ------------------------------" << endl;
            
            N_true_roots_lost+=(N_roots_found_this_test<N_roots_gt_this_test)*(N_roots_gt_this_test-N_roots_found_this_test);
            N_fake_roots_added+=(N_roots_found_this_test>N_roots_gt_this_test)*(N_roots_found_this_test-N_roots_gt_this_test);
            if (ae_this_test_worst>ae_all_tests_worst && ( N_roots_found_this_test>=N_roots_gt_this_test))
            {
                roots_gt_all_tests_ae_worst=roots; N_roots_gt_all_tests_ae_worst=N_roots_gt_this_test;
                roots_found_all_tests_ae_worst=roots_found_this_test; N_roots_found_all_tests_ae_worst=N_roots_found_this_test;
                coefficients_all_tests_ae_worst=a; ae_all_tests_worst=ae_this_test_worst;
            }
            if (re_this_test_worst>re_all_tests_worst && (N_roots_found_this_test>=N_roots_gt_this_test))
            {
                roots_gt_all_tests_re_worst=roots; N_roots_gt_all_tests_re_worst=N_roots_gt_this_test;
                roots_found_all_tests_re_worst=roots_found_this_test; N_roots_found_all_tests_re_worst=N_roots_found_this_test;
                coefficients_all_tests_re_worst=a; re_all_tests_worst=re_this_test_worst;
            }
            if(Laguerre::anynotfinite<NUMBER>(roots_found_this_test)){
                cout << "\nNOT FINITE ROOTS!!!\n";
                delete solver;
                return 1;
            }
            if(Laguerre::anynotfinite<NUMBER>(a)){
                cout << "NOT FINITE COEFFS!!!\n";
                delete solver;
                return 1;
            }
            fill(roots.begin(), roots.end(), 0.0);
            fill(a.begin(), a.end(), 0.0);

            cout << "\n=======================\n";
        }

    }
    catch (const invalid_argument &exc)
    {
        cerr << "Bad argument: " << exc.what();
    }
    catch (const exception &exc)
    {
        cerr << exc.what();
    }
    cout << "\nDONE!\n";
    cout << endl << "--------------- settings ----------------" << endl
        << "P=" << DEGREE << ", N_tests=" << N_TESTS; 
        
    cout << endl << endl << "--------------- statistics ----------------" << endl
        << ", N_true_roots_lost=" << N_true_roots_lost << ", N_fake_roots_added=" << N_fake_roots_added
        << endl;

    // output the worst cases
    cout << endl << "--------------- worst case: absolute error: ----------------" << endl
        << "ae_all_tests_worst=" << ae_all_tests_worst
        << ", N_roots_gt_all_tests_ae_worst=" << N_roots_gt_all_tests_ae_worst
        << ", N_roots_found_all_tests_ae_worst=" << N_roots_found_all_tests_ae_worst << endl;

    cout << endl << "--------------- worst case: relative error: ----------------" << endl
        << "re_all_tests_worst=" << re_all_tests_worst
        << ", N_roots_gt_all_tests_re_worst=" << N_roots_gt_all_tests_re_worst
        << ", N_roots_found_all_tests_re_worst=" << N_roots_found_all_tests_re_worst << endl;

    cout << "\nWorst case method results ae: ";
    for (auto &root : roots_found_all_tests_ae_worst){
        cout << root << " ";
    }
    cout << "\nGround truth worst ae: {";
    for (auto &root : roots_gt_all_tests_ae_worst){
        cout << root << ", ";
    }
    cout << "}\nCoeffs worst ae: {";
    for (auto &root : coefficients_all_tests_ae_worst){
        cout << root << ", ";
    }

    cout << "}\nWorst case method results re: ";
    for (auto &root : roots_found_all_tests_re_worst){
        cout << root << " ";
    }
    cout << "\nGround truth worst re: {";
    for (auto &root : roots_gt_all_tests_re_worst){
        cout << root << ", ";
    }
    cout << "}\nCoeffs worst re: {";
    for (auto &root : coefficients_all_tests_re_worst){
        cout << root << ", ";
    }   
    delete solver;
    return 0;
}