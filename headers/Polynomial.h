// Александр, Дмитрий

#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <complex>
#include <iostream>

#include "Polynomial.h"
#include "BaseSolver.h"
#include "ExtendedFunctions.h"

using std::cout;

namespace Laguerre{
template<typename T>
class Polynomial{
protected:
    vector<complex<T>> roots;
    vector<T> coeffs;

    // Pointer to a base solver class
    BaseSolver<T>* Solver;
public:

    Polynomial(){}

    /**
     * \brief Create polynomial object from given vectors of coefficients and roots.
         * \param _coeffs Vector of coefficients.
         * \param _roots Vector of roots.
     */
    Polynomial(vector<T> _coeffs, vector<complex<T>> _roots){
        setCoeffs(_coeffs);
        setRoots(_roots);
    }

    /**
     * \brief Create polynomial object from given vector of coefficients.
         * \param args Vector of coefficients.
    */
    Polynomial(vector<T> args){
        setCoeffs(args);
    }

    /**
     * \brief Create polynomial object from given coefficients.
         * \param args Coefficients.
    */
    template<typename... Args>
    Polynomial(Args...args){
        setCoeffs(args...);
    }

    /**
     * \brief Get polynomials degree.
         * \returns Degree of the polynomial.
    */
    int degree() const {
        int size = coeffs.size();
        return size > 0 ? size - 1 : 0;
    }

    /**
     * \brief Set roots for this polynomial.
         * \param args Vector of roots.
     */
    void setRoots(vector<complex<T>> args){
        size_t count = args.size();
        if(count != degree())
            throw invalid_argument("Count of given roots is different from given polynomial degree. (" 
                    + std::to_string(count) + " vs " + std::to_string(degree()) + ")\n");
        roots = args;
    }

    /**
     * \brief Set roots for this polynomial.
         * \param args Given roots.
     */
    template<typename... Args>
    void setRoots(Args...args){
        setRoots(vector<complex<T>>{args...});
    }

    /**
     * \brief Set coefficients for this polynomial.
         * \param args Vector of coefficients.
     */
    void setCoeffs(vector<T> args){
        coeffs = args;
        roots.clear();
    }

    /**
     * \brief Set coefficients for this polynomial.
         * \param args Given coefficients.
     */
    template<typename... Args>
    void setCoeffs(Args...args){
        setCoeffs(vector<T>{args...});
    }

    /**
     * \brief Set solver for this polynomial.
         * \param solver Pointer to a Solver object.
     */
    void setSolver(BaseSolver<T>* solver) {
        Solver = solver;
    }

    /**
     * \brief Find the roots of this polynomial. NOTE: Set solver first, use function setSolver
         * \param roots Vector to store the roots.
         * \param conv Vector to store convergence status of each root.
         * \param itmax Maximum number of iterations.
     */
    void solve(vector<complex<T>>& roots, vector<int>& conv, int maxiter=80){
        if(Solver) {
            (*Solver)(coeffs, roots, conv, maxiter);
            // cout << "solved, exiting";
        }
        else{
            throw invalid_argument("Solver wasnt set!\n");
        }
    }

    Polynomial operator+(Polynomial& rhs){
        vector<T> result_coeffs(std::max(this->degree(), rhs.degree()) + 1, 0);
        for (int i = 0; i <= this->degree(); ++i) {
            result_coeffs[i] += this->coeffs[i];
        }

        for (int i = 0; i <= rhs.degree(); ++i) {
            result_coeffs[i] += rhs.coeffs[i];
        }

        return Polynomial<T>(result_coeffs);
    }

    Polynomial operator-(const Polynomial& rhs) const {
        vector<T> result_coeffs(std::max(this->degree(), rhs.degree()) + 1, 0);
        for (int i = 0; i <= this->degree(); ++i) {
            result_coeffs[i] += this->coeffs[i];
        }
        for (int i = 0; i <= rhs.degree(); ++i) {
            result_coeffs[i] -= rhs.coeffs[i];
        }
        return Polynomial<T>(result_coeffs);
    }

    Polynomial operator*(const Polynomial& other){
        int n = coeffs.size();
        int m = other.coeffs.size();
        vector<T> result(n + m - 1, 0);

        for (int i = 0; i < n; ++i) {
            auto this_coeff = coeffs[i];
            for (int j = 0; j < m; ++j) {
                // result[i + j] += this_coeff * other.coeffs[j];
                result[i + j] = fma(this_coeff, other.coeffs[j], result[i + j]);
            }
        }
        return Polynomial(result);
    }

    Polynomial operator/(Polynomial& divisor) {
        if (divisor.degree() > this->degree()) {
            throw invalid_argument("The degree of the divisor is greater than the dividend");
        }

        Polynomial<T> dividend = *this;
        vector<T> quotient_coeffs(this->degree() - divisor.degree() + 1, 0);

        T coeff;
        int degree_diff;

        while (dividend.degree() >= divisor.degree()) {
            degree_diff = dividend.degree() - divisor.degree();
            coeff = dividend.coeffs.back() / divisor.coeffs.back();
            if(anynotfinite(coeff)){
                throw invalid_argument("Highest degree of coefficient for divisor is probably zero: " << divisor.coeffs.back());
            }

            // Subtract and update the dividend
            for (int i = 0; i <= divisor.degree(); ++i) {
                // dividend.coeffs[i + degree_diff] -= coeff * divisor.coeffs[i];
                dividend.coeffs[i + degree_diff] = fma(-coeff, divisor.coeffs[i], dividend.coeffs[i + degree_diff]);
            }
            dividend.coeffs.pop_back();

            // Update the quotient
            quotient_coeffs[degree_diff] = coeff;
        }

        return Polynomial<T>(quotient_coeffs);
    }


    /**
     * \brief Divides one polynomial by another.
     * 
     * This function performs polynomial division. Given a divisor and a dividend (the instance calling the function),
     * it calculates the quotient and the remainder of the division.
     * 
     * \param divisor The polynomial to divide by.
     * \param quotient The result of the division.
     * \param remainder The remainder after division.
     */
    void divide(Polynomial& divisor, Polynomial& quotient, Polynomial& remainder) {
        if (divisor.degree() > this->degree()) {
            throw invalid_argument("The degree of the divisor is greater than the dividend");
        }

        Polynomial<T> dividend = *this;
        vector<T> quotient_coeffs(this->degree() - divisor.degree() + 1, 0);
        Polynomial<T> tmp;
        int degree_diff;
        T coeff;

        while (dividend.degree() >= divisor.degree()) {
            degree_diff = dividend.degree() - divisor.degree();
            coeff = dividend.coeffs.back() / divisor.coeffs.back();

            if(anynotfinite(coeff)){
                throw invalid_argument("Highest degree of coefficient for divisor is probably zero: " << divisor.coeffs.back());
            }

            // Update the temporary polynomial for subtraction
            tmp.coeffs = vector<T>(degree_diff + 1, 0);
            tmp.coeffs.back() = coeff;

            // Subtract and update the dividend
            for (int i = 0; i <= divisor.degree(); ++i) {
                // dividend.coeffs[i + degree_diff] -= coeff * divisor.coeffs[i];
                dividend.coeffs[i + degree_diff] = fma(-coeff, divisor.coeffs[i], dividend.coeffs[i + degree_diff]);
            }
            dividend.coeffs.pop_back();

            // Update the quotient
            quotient_coeffs[degree_diff] = coeff;
        }

        quotient = Polynomial<T>(quotient_coeffs);
        remainder = dividend;
    }

    Polynomial& operator*=(const Polynomial& other) {
        Polynomial temp(*this); // Temp copy so multiplying on self wont break polynomial
        *this = temp * other;
        this->roots.clear();
        return *this;
    }

    Polynomial& operator=(const Polynomial& other){
        if (this != &other) {
            this->setCoeffs(other.coeffs);
        }
        return *this;
    }

    T operator[](int i) const{
        return coeffs[i];
    }

    T& operator [](int i){
        return coeffs[i];
    }

    /**
     * \brief Calculates the derivative of a Polynomial object.
     * 
     * This function calculates the derivative of a given degree of the Polynomial object. 
     * The default degree of the derivative is 1, which means it calculates the first derivative if no degree is specified.
     * 
     * \tparam T The data type of the Polynomial coefficients.
     * 
     * \param deg The degree of the derivative. Default is 1.
     * 
     * \return A vector of coefficients representing the derivative of the Polynomial.
     */
    vector<T> diff(int deg=1){
        vector ret = coeffs;
        for(int j = 0; j < deg; ++j){
            for (int i = 1; i < ret.size(); ++i) {
                ret[i] *= static_cast<T>(i);
            }
            if (!ret.size()) {
                ret.erase(ret.begin());
            }
            else break;
        }
        return ret;
    }

    void print(){
        int _degree = degree();
        if(!_degree){
            cout << "This polynomial object doesnt contain any coefficients.\n";
        }
        else{
            for(int i = _degree; i >= 0; --i){
                cout << "(" << coeffs[i] << ")";
                cout << "x^" << i;
                if(i != 0) {
                    cout << " + ";
                }
            }
            cout << ". Degree: " << _degree << '\n';
            cout << "Roots: " << roots.size() <<  "; ";
            printVec(roots);
        }
    }       
};
};

#endif