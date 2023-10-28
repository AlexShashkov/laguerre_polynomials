#include "headers/Polynomial.h"
#include "headers/ExtendedFunctions.h"
#include "headers/PolynomialGenerator.h"

#define number double

int main(){
    using Laguerre::Polynomial;
    try{
        Polynomial<number> genRoots = Laguerre::Generator<number>::createFromRoots(1.f, 2.f, 3.f);
        genRoots.print();
        Polynomial<number> genRoots2 = Laguerre::Generator<number>::gen(4);
        genRoots2.print();
        std::cout << genRoots2[2] << "\n";
        genRoots2[2] = 0.0001;
        std::cout << genRoots2[2] << "\n";
    }
    catch (const std::invalid_argument &exc)
    {
        std::cerr << "Bad argument: " << exc.what();
    }
    catch (const std::exception &exc)
    {
        std::cerr << exc.what();
    }
    catch ( ... )
    {
        std::cerr << "UNKNOWN EXCEPTION!\n";
    }
    return 0;
}
