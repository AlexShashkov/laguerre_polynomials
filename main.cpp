#include "framework/framework.h"

#include "headers/Polynomial.h"
#include "headers/PolynomialGenerator.h"

#include "headers/Laguerre.h"
#include "headers/Laguerre18m.h"

#define number double

#define PR_NUMBERS_OF_ROOTS_EQUAL     0
#define PR_AT_LEAST_ONE_ROOT_LOST    -1
#define PR_AT_LEAST_ONE_ROOT_IS_FAKE -2
#define PR_2_INFINITE_ROOTS          -3

int main(){
    using Laguerre::Polynomial;
    // Laguerre solvers
    Laguerre::Original<double>* solver = new Laguerre::Original<double>();
    Laguerre::ModifiedLaguerre18<double>* solver18 = new Laguerre::ModifiedLaguerre18<double>();

    // Generator stuff
    int l = 3, N_TESTS = 1000000,                                               // polynomial degree, count of tests 
        rv,                                                                  // status of comparing roots
        N_roots_found_this_test, N_roots_gt_this_test,                       // amount of found roots in each test & gt roots
        N_true_roots_lost=0, N_fake_roots_added=0,                           // total counters of {lost, fake} roots aover all tests
        N_roots_gt_all_tests_ae_worst=0, N_roots_found_all_tests_ae_worst=0, // {ground truth, found} number of roots of the worst polynomial for absolute error
        N_roots_gt_all_tests_re_worst=0, N_roots_found_all_tests_re_worst=0; // {ground truth, found} number of roots of the worst polynomial for relative error

    std::vector<number> roots(l, 0.0), roots_found_this_test(l), a(l+1, 0.0),                                       // current polynomial
        roots_gt_all_tests_ae_worst(l), roots_found_all_tests_ae_worst(l), coefficients_all_tests_ae_worst(l+1),    // {gt roots, found roots, coefficients} of the worst polynomial for absolute error
        roots_gt_all_tests_re_worst(l), roots_found_all_tests_re_worst(l), coefficients_all_tests_re_worst(l+1);    // {gt roots, found roots, coefficients} of the worst polynomial for relative error

    number ae_this_test_worst, ae_all_tests_worst=static_cast<number>(-1.0L),   // absolute error: an impossible value that will be updated for sure
           re_this_test_worst, re_all_tests_worst=static_cast<number>(-1.0L);   // relative error: an impossible value that will be updated for sure
    try{
        for(int i=0; i < N_TESTS; ++i){
            roots_found_this_test.clear();
            N_roots_gt_this_test = generate_polynomial(l, 0, 3, 0, 1.0, 0.0, 20.0, roots, a);
            std::cout << "GENERATED POLY & ROOTS:\n";
            Laguerre::printVec(a);
            Laguerre::printVec(roots);
            Polynomial<number> pol(a);

            // Очень странная штука - в цикле внезапно может оказаться, что вектор
            // для корней размером n крашится, хотя проверено, что нигде за массив не записывает;
            // повторить вне цикла не получается - те же коэффициенты спокойно решаются
            // с тем же правильным размером вектора n, причем сам последний элемент n+1 НЕ ТРОГАЕТ;
            // костыль, но откуда копать причину ошибки пока не знаю, просто "free(): invalid pointer"
            // по завершению функции случайно. Конечно же необходимо исправить, лишнее место ни к чему
            //  - Саша
            std::vector<std::complex<double>> solved_roots(l+1);
            std::vector<int> conv(l);

            std::cout << "ORIGINAL LAGUERRE:\n";
            pol.setSolver(solver);
            pol.solve(solved_roots, conv, 80);
            Laguerre::printVec(solved_roots);

            // std::cout << "2018 LAGUERRE MODIFICATION:\n";
            // pol.setSolver(solver18);
            // pol.solve(solved_roots, conv, 80);
            // Laguerre::printVec(solved_roots);

            std::cout << "=======================\n";

            //                                                    ! убрать -1 когда будет найдено решение ошибки с roots
            std::for_each(solved_roots.begin(), solved_roots.end()-1, [&roots_found_this_test](std::complex<number> x){ // solved_roots = b_roots
                roots_found_this_test.push_back(x.real());
            });

            N_roots_found_this_test = roots_found_this_test.size();

            rv=compare_roots<number>(N_roots_found_this_test, N_roots_gt_this_test, roots_found_this_test, roots,
                                                                                    ae_this_test_worst, re_this_test_worst);


            std::cout << std::endl << "P=" << l << ", test No " << i << " out of " << N_TESTS << std::endl
                << "comparison return value=" << (rv==PR_NUMBERS_OF_ROOTS_EQUAL ? "PR_NUMBERS_OF_ROOTS_EQUAL" :
                (rv==PR_AT_LEAST_ONE_ROOT_LOST ? "PR_AT_LEAST_ONE_ROOT_LOST" : 
                (rv==PR_AT_LEAST_ONE_ROOT_IS_FAKE ? "PR_AT_LEAST_ONE_ROOT_IS_FAKE" : "UNKNOWN VALUE")))
                << ", N_roots_gt_this_test=" << N_roots_gt_this_test << ", N_roots_found_this_test=" << N_roots_found_this_test
                << ", N_roots_verified_this_test=" << l
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
        << "P=" << l << ", N_tests=" << N_TESTS; // << ", N_pairs_of_complex_roots=" << N_pairs_of_complex_roots
        // << ", N_clustered_roots=" << N_clustered_roots << ", N_multiple_roots=" << N_multiple_roots << std::endl
        // << "max_distance_between_clustered_roots=" << max_distance_between_clustered_roots
        // << ", root_sweep=[" << root_sweep_low << "..." << root_sweep_high << "]\ncomplex_value_tolerance=" << complex_value_tolerance
        // << ", value_at_root_tolerance=" << value_at_root_tolerance << ", exclude_lost_root_cases_from_statistics=" << exclude_lost_root_cases_from_statistics << std::endl
        //<< "pr_lost_roots_saving_1=" << pr_lost_roots_saving_1 << ", pr_lost_roots_saving_2=" << pr_lost_roots_saving_2;
    // output total statistics
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
    delete solver;
    delete solver18;
    return 0;
}