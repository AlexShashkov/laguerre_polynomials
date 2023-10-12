#include "headers/Polynomial.h"
#include "headers/PolynomialGenerator.h"

int main(){
    Laguerre::Polynomial<float> test(-1., 2., 3., 4., 5.);
    test.print();
    Laguerre::Polynomial<float> test2(-1., 2., 3., 4.);
    test *= test2;
    test.print();
    Laguerre::Polynomial<float> genRoots = Laguerre::Generators<float>::createFromRoots(1., 2., 3.);
    genRoots.print();
    std::cout << "Diff res: ";
    genRoots.diff().print(); 
    try{
        test.setRoots(2, 2, 8, 9);
        test.setCoeffs();
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

    test.print();
    return 0;
}
