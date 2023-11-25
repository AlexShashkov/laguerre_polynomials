#include "headers/Polynomial.h"
#include "headers/PolynomialGenerator.h"

//#include "headers/Laguerre18m.h"
#include "headers/Laguerre.h"
#include "headers/Laguerre13.h"

#define number double

int main(){
    using Laguerre::Polynomial;
    Laguerre::Original<double>* solver = new Laguerre::Original<double>();
    Laguerre::Laguerre13<double>* solver13 = new Laguerre::Laguerre13<double>();
    //FLaguerre::ModifiedLaguerre18<double>* solver18 = new Laguerre::ModifiedLaguerre18<double>();
    try{
        std::cout << "ORIGINAL LAGUERRE:\n";
        Polynomial<number> genRoots = Laguerre::Generator<number>::createFromRoots(3.f, 3.f, 3.f, 6.f, 7.f, 7.f, 7.f);
        genRoots.print();
        int deg = genRoots.degree();
        std::vector<std::complex<double>> roots(deg);
        std::vector<int> conv(deg);

        // Set Laguerre solver for Polynomial and solve it
        
        genRoots.setSolver(solver);
        genRoots.solve(roots, conv, 80);
        Laguerre::printVec(roots);
        Laguerre::printVec(conv);

        std::cout << "ORIGINAL LAGUERRE:\n";
        genRoots.setSolver(solver13);
        genRoots.solve(roots, conv, 1000);
        Laguerre::printVec(roots);
        Laguerre::printVec(conv);
        /*
        Polynomial<number> genRoots2 = Laguerre::Generator<number>::gen(4);
        genRoots2.print();
        genRoots2.setSolver(solver);

        deg = genRoots2.degree();
        std::vector<std::complex<double>> roots2(deg);
        std::vector<int> conv2(deg);

        genRoots2.solve(roots2, conv2, 80);
        Laguerre::printVec(roots2);
        Laguerre::printVec(conv2);
        std::cout << "2018 LAGUERRE MODIFICATION:\n";

        deg = genRoots.degree(); 
        roots = std::vector<std::complex<double>>(deg);
        conv = std::vector<int>(deg);

        deg = genRoots2.degree(); 
        roots2 = std::vector<std::complex<double>>(deg);
        conv2 = std::vector<int>(deg);

        genRoots.print();
        genRoots.setSolver(solver18);
        genRoots.solve(roots, conv, 80);
        Laguerre::printVec(roots);
        Laguerre::printVec(conv);

        genRoots2.print();
        genRoots2.setSolver(solver18);
        genRoots2.solve(roots2, conv2, 80);
        Laguerre::printVec(roots2);
        Laguerre::printVec(conv2);*/
    }
    catch (const std::invalid_argument &exc)
    {
        std::cerr << "Bad argument: " << exc.what();
    }
    catch (const std::exception &exc)
    {
        std::cerr << exc.what();
    }
    delete solver13;
    delete solver;
    //delete solver18;
    return 0;
}