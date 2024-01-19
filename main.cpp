#include "framework/framework.h"

#include "headers/Polynomial.h"
#include "headers/PolynomialGenerator.h"

#include "headers/Laguerre.h"
#include "headers/Laguerre18m.h"

#define number double

int main(){
    using Laguerre::Polynomial;
    // Laguerre solvers
    Laguerre::Original<double>* solver = new Laguerre::Original<double>();
    Laguerre::ModifiedLaguerre18<double>* solver18 = new Laguerre::ModifiedLaguerre18<double>();
    try{
        // Generator stuff
        int l = 4, cnt = 100;
        std::vector<double> roots(l, 0.0);
        std::vector<double> a(l+1, 0.0);
        for(int i=0; i < cnt; ++i){
            std::cout << "ITER #" << i << "\n";
            generate_polynomial(l, 0, 3, 0, 1.0, 0.0, 20.0, roots, a);
            std::cout << "GENERATED POLY & ROOTS:\n";
            Laguerre::printVec(a);
            Laguerre::printVec(roots);
            Polynomial<number> pol(a);
            pol.print();
            Laguerre::printVec(roots);

            std::cout << "ORIGINAL LAGUERRE:\n";

            pol.setSolver(solver);
            auto deg = pol.degree();
            std::vector<std::complex<double>> roots2(deg);
            std::vector<int> conv(deg);

            std::cout << "CREATED original - " << deg << " roots\n";

            pol.solve(roots2, conv, 80);
            std::cout << "Original solved!\n";
            Laguerre::printVec(roots2);

            std::cout << "2018 LAGUERRE MODIFICATION:\n";

            pol.setSolver(solver18);
            pol.solve(roots2, conv, 80);
            std::cout << "Modification solved!\n";
            Laguerre::printVec(roots2);
            std::cout << "=======================\n";

            std::fill(roots.begin(), roots.end(), 0);
            std::fill(a.begin(), a.end(), 0);
        }
        std::cout << "DONE!\n";
    }
    catch (const std::invalid_argument &exc)
    {
        std::cerr << "Bad argument: " << exc.what();
    }
    catch (const std::exception &exc)
    {
        std::cerr << exc.what();
    }
    delete solver;
    delete solver18;
    return 0;
}