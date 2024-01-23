#include "headers/framework.h"

#include "headers/Polynomial.h"
#include "headers/PolynomialGenerator.h"

#include "headers/Laguerre.h"
#include "headers/Laguerre13m.h"
#include "headers/Laguerre18m.h"
// #include "headers/L18_or.h"

#define PR_NUMBERS_OF_ROOTS_EQUAL     0
#define PR_AT_LEAST_ONE_ROOT_LOST    -1
#define PR_AT_LEAST_ONE_ROOT_IS_FAKE -2
#define PR_2_INFINITE_ROOTS          -3

#ifndef SOLVER
    #define SOLVER 1
    // 0 - Original
    // 1 - 2013 mod
    // 2 - 2018 mod
#endif
#ifndef NUMBER
    #define NUMBER float
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
    #define EXPONENT 4 // Exponent for big NUMBER notation
#endif
#ifndef MANTISSA
    #define MANTISSA 32 // Mantissa for big NUMBER notation
#endif


#ifndef DEGREE
    #define DEGREE 30 // polynomial degree
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
    #define N_MULTIPLE_ROOTS DEGREE
#endif
// Used only for clustered roots
#ifndef MAX_DISTANCE_BETWEEN_CLUSTERED
    #define MAX_DISTANCE_BETWEEN_CLUSTERED 1e-5
#endif
#ifndef ROOT_SWEEP_LOW
    #define ROOT_SWEEP_LOW 5.0
#endif
#ifndef ROOT_SWEEP_HIGH
    #define ROOT_SWEEP_HIGH 10.0
#endif
#ifndef TRIES
    #define TRIES 80 // For how many iterations Laguerre should trying to solve polynomial
#endif

int main(){
    using Laguerre::Polynomial;
    std::streamsize fp_precision_original=std::cout.precision(); // save default precision to provide maximal reasonable output
    
    // Laguerre solvers
    #if SOLVER == 0
        std::cout << "Using original version of Laguerre\n";
        Laguerre::Original<NUMBER>* solver = new Laguerre::Original<NUMBER>();
    #elif SOLVER == 1
        std::cout << "Using modified version of Laguerre, 2013\n";
        Laguerre::ModifiedLaguerre13<NUMBER>* solver = new Laguerre::ModifiedLaguerre13<NUMBER>();
    #elif SOLVER == 2
        std::cout << "Using modified version of Laguerre, 2018\n";
        Laguerre::ModifiedLaguerre18<NUMBER>* solver = new Laguerre::ModifiedLaguerre18<NUMBER>();
    #endif

    // Generator stuff
    int rv,                                                                  // status of comparing roots
        N_roots_found_this_test, N_roots_gt_this_test,                       // amount of found roots in each test & gt roots
        N_true_roots_lost=0, N_fake_roots_added=0,                           // total counters of {lost, fake} roots aover all tests
        N_roots_gt_all_tests_ae_worst=0, N_roots_found_all_tests_ae_worst=0, // {ground truth, found} NUMBER of roots of the worst polynomial for absolute error
        N_roots_gt_all_tests_re_worst=0, N_roots_found_all_tests_re_worst=0; // {ground truth, found} NUMBER of roots of the worst polynomial for relative error

    std::vector<NUMBER> roots(DEGREE, 0.0), roots_found_this_test(DEGREE), a(DEGREE+1, 0.0),                                       // current polynomial
        roots_gt_all_tests_ae_worst(DEGREE), roots_found_all_tests_ae_worst(DEGREE), coefficients_all_tests_ae_worst(DEGREE+1),    // {gt roots, found roots, coefficients} of the worst polynomial for absolute error
        roots_gt_all_tests_re_worst(DEGREE), roots_found_all_tests_re_worst(DEGREE), coefficients_all_tests_re_worst(DEGREE+1);    // {gt roots, found roots, coefficients} of the worst polynomial for relative error

    NUMBER ae_this_test_worst, ae_all_tests_worst=static_cast<NUMBER>(-1.0L),   // absolute error: an impossible value that will be updated for sure
           re_this_test_worst, re_all_tests_worst=static_cast<NUMBER>(-1.0L);   // relative error: an impossible value that will be updated for sure
    try{
        for(int i=0; i < N_TESTS; ++i){
            roots_found_this_test.clear();
            N_roots_gt_this_test = generate_polynomial<NUMBER, EXPONENT, MANTISSA>(DEGREE, N_PAIRS_OF_COMPLEX_ROOTS, N_CLUSTERED_ROOTS,
                N_MULTIPLE_ROOTS, MAX_DISTANCE_BETWEEN_CLUSTERED, ROOT_SWEEP_LOW, ROOT_SWEEP_HIGH, roots, a);
            std::cout << "GENERATED POLY & ROOTS:\n";
            Laguerre::printVec(a);
            Laguerre::printVec(roots);
            Polynomial<NUMBER> pol(a);

            std::vector<std::complex<NUMBER>> solved_roots(DEGREE);
            std::vector<int> conv(DEGREE);

            pol.setSolver(solver);
            pol.solve(solved_roots, conv, TRIES);

            std::cout << "found roots: ";
            std::for_each(solved_roots.begin(), solved_roots.end(), [&roots_found_this_test](std::complex<NUMBER> x){ // solved_roots = b_roots
                std::cout << x.real() << "; ";
                roots_found_this_test.push_back(x.real());
            });

            N_roots_found_this_test = roots_found_this_test.size();

            rv=compare_roots<NUMBER>(N_roots_found_this_test, N_roots_gt_this_test, roots_found_this_test, roots,
                                                                                    ae_this_test_worst, re_this_test_worst);


            std::cout << std::endl << "P=" << DEGREE << ", test No " << i << " out of " << N_TESTS << std::endl
                << "comparison return value=" << (rv==PR_NUMBERS_OF_ROOTS_EQUAL ? "PR_NUMBERS_OF_ROOTS_EQUAL" :
                (rv==PR_AT_LEAST_ONE_ROOT_LOST ? "PR_AT_LEAST_ONE_ROOT_LOST" : 
                (rv==PR_AT_LEAST_ONE_ROOT_IS_FAKE ? "PR_AT_LEAST_ONE_ROOT_IS_FAKE" : "UNKNOWN VALUE")))
                << ", N_roots_gt_this_test=" << N_roots_gt_this_test << ", N_roots_found_this_test=" << N_roots_found_this_test
                << ", N_roots_verified_this_test=" << DEGREE
                << "ae_this_test_worst=" << ae_this_test_worst
                << ", re_this_test_worst=" << re_this_test_worst << " ------------------------------" << std::endl;
            
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
            std::fill(roots.begin(), roots.end(), 0.0);
            std::fill(a.begin(), a.end(), 0.0);
            std::cout << "\n=======================\n";
        }

    }
    catch (const std::invalid_argument &exc)
    {
        std::cerr << "Bad argument: " << exc.what();
    }
    catch (const std::exception &exc)
    {
        std::cerr << exc.what();
    }
    std::cout << "DONE!\n";
    std::cout << std::endl << "--------------- settings ----------------" << std::endl
        << "P=" << DEGREE << ", N_tests=" << N_TESTS; 
        
    std::cout << std::endl << std::endl << "--------------- statistics ----------------" << std::endl
        << "N_true_roots_lost=" << N_true_roots_lost << ", N_fake_roots_added=" << N_fake_roots_added
        << std::endl;

    // output the worst cases
    std::cout << std::endl << "--------------- worst case: absolute error: ----------------" << std::endl
        << "ae_all_tests_worst=" << ae_all_tests_worst
        << ", N_roots_gt_all_tests_ae_worst=" << N_roots_gt_all_tests_ae_worst
        << ", N_roots_found_all_tests_ae_worst=" << N_roots_found_all_tests_ae_worst << std::endl;

    std::cout << std::endl << "--------------- worst case: relative error: ----------------" << std::endl
        << "re_all_tests_worst=" << re_all_tests_worst
        << ", N_roots_gt_all_tests_re_worst=" << N_roots_gt_all_tests_re_worst
        << ", N_roots_found_all_tests_re_worst=" << N_roots_found_all_tests_re_worst << std::endl;

    std::cout << "\nWorst case method results ae: ";
    for (auto &root : roots_found_all_tests_ae_worst){
        std::cout << root << " ";
    }
    std::cout << "\nGround truth worst ae: {";
    for (auto &root : roots_gt_all_tests_ae_worst){
        std::cout << root << ", ";
    }
    std::cout << "}\nCoeffs worst ae: {";
    for (auto &root : coefficients_all_tests_ae_worst){
        std::cout << root << ", ";
    }


    std::cout << "}\nWorst case method results re: ";
    for (auto &root : roots_found_all_tests_re_worst){
        std::cout << root << " ";
    }
    std::cout << "\nGround truth worst re: {";
    for (auto &root : roots_gt_all_tests_re_worst){
        std::cout << root << ", ";
    }
    std::cout << "}\nCoeffs worst re: {";
    for (auto &root : coefficients_all_tests_re_worst){
        std::cout << root << ", ";
    }
    std::cout << "}" << std::endl;
    std::cout << std::setprecision(fp_precision_original); // restore default precision
    
    delete solver;
    return 0;
}