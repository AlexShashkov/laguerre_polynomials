#include "headers/Polynomial.h"
#include "headers/PolynomialGenerator.h"
#include "headers/Laguerre18m.h"

#define number double

int main(){
    using Laguerre::Polynomial;
    try{
        Polynomial<number> genRoots = Laguerre::Generator<number>::createFromRoots(1.f, 2.f, 3.f);
        genRoots.print();
        Polynomial<number> genRoots2 = Laguerre::Generator<number>::gen(4);
        genRoots2.print();

        int deg = genRoots.degree();
        std::vector<std::complex<double>> roots(deg);
        std::vector<double> berr(deg);
        std::vector<double> cond(deg);
        std::vector<int> conv(deg);
        Laguerre::ModifiedLaguerre18<double> solver18;
        solver18(genRoots, roots, berr, cond, conv, 80);
        Laguerre::printVec(roots);
        Laguerre::printVec(berr);
        Laguerre::printVec(cond);
        Laguerre::printVec(conv);
    }
    catch (const std::invalid_argument &exc)
    {
        std::cerr << "Bad argument: " << exc.what();
    }
    catch (const std::exception &exc)
    {
        std::cerr << exc.what();
    }
    return 0;
}
