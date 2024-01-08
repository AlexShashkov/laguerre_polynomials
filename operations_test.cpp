#include "headers/Polynomial.h"
#include "headers/PolynomialGenerator.h"

#define number double

int main(){
    using Laguerre::Polynomial;
    try{
        // Polynomial<number> divident(-42.f, 0.f, -12.f, 1.f);
        // Polynomial<number> divisor(-3.f, 1.f);
        Polynomial<number> divident(1.f, 4.f, 2.f, -3.f, 1.f);
        Polynomial<number> divisor(-2.f, 0.f, 1.f);
        divident.print();
        divisor.print();

        Polynomial<number> quotient, remainder;

        Polynomial<number> res = divident / divisor;
        res.print();

        std::cout << "Qoutient and remainder:\n";
        divident.divide(divisor, quotient, remainder);
        quotient.print();
        remainder.print();
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