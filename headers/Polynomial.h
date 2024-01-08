#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <complex>
#include <iostream>

#include "Polynomial.h"
#include "BaseSolver.h"
#include "ExtendedFunctions.h"

namespace Laguerre{
template<typename T>
class Polynomial{
protected:
    std::vector<std::complex<T>> roots;
    std::vector<T> coeffs;

    // Pointer to a base solver class
    BaseSolver<T>* Solver;
public:

    /**
     * \brief Create polynomial object from given vectors of coefficients and roots.
         * \param _coeffs Vector of coefficients.
         * \param _roots Vector of roots.
     */
    Polynomial(std::vector<T> _coeffs, std::vector<std::complex<T>> _roots){
        setCoeffs(_coeffs);
        setRoots(_roots);
    }

    /**
     * \brief Create polynomial object from given vector of coefficients.
         * \param args Vector of coefficients.
    */
    Polynomial(std::vector<T> args){
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
    void setRoots(std::vector<std::complex<T>> args){
        size_t count = args.size();
        if(count != degree())
            throw std::invalid_argument("Count of given roots is different from given polynomial degree. (" 
                    + std::to_string(count) + " vs " + std::to_string(degree()) + ")\n");
        roots = args;
    }

    /**
     * \brief Set roots for this polynomial.
         * \param args Given roots.
     */
    template<typename... Args>
    void setRoots(Args...args){
        setRoots(std::vector<std::complex<T>>{args...});
    }

    /**
     * \brief Set coefficients for this polynomial.
         * \param args Vector of coefficients.
     */
    void setCoeffs(std::vector<T> args){
        coeffs = args;
        roots.clear();
    }

    /**
     * \brief Set coefficients for this polynomial.
         * \param args Given coefficients.
     */
    template<typename... Args>
    void setCoeffs(Args...args){
        setCoeffs(std::vector<T>{args...});
    }

    /**
     * \brief Set solver for this polynomial.
         * \param solver Pointer to a Solver object.
     */
    void setSolver(BaseSolver<T>* solver) {
        Solver = solver;
    }

    /**
     * \brief Find the roots of this polynomial.
         * \param roots Vector to store the roots.
         * \param conv Vector to store convergence status of each root.
         * \param itmax Maximum number of iterations.
     */
    void solve(std::vector<std::complex<T>>& roots, std::vector<int>& conv, int maxiter=80){
        if(Solver) {
            (*Solver)(coeffs, roots, conv, 80);
        }
        else{
            throw std::invalid_argument("Solver wasnt set!\n");
        }
    }

    Polynomial operator+(Polynomial& rhs){
        std::vector<T> result_coeffs(std::max(this->degree(), rhs.degree()) + 1, 0);
        for (int i = 0; i <= this->degree(); ++i) {
            result_coeffs[i] += this->coeffs[i];
        }

        for (int i = 0; i <= rhs.degree(); ++i) {
            result_coeffs[i] += rhs.coeffs[i];
        }

        return Polynomial<T>(result_coeffs);
    }

    Polynomial operator-(const Polynomial& rhs) const {
        std::vector<T> result_coeffs(std::max(this->degree(), rhs.degree()) + 1, 0);
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
        std::vector<T> result(n + m - 1, 0);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                result[i + j] += coeffs[i] * other.coeffs[j];
            }
        }
        return Polynomial(result);
    }

    Polynomial operator/(Polynomial& divisor) {
        if (divisor.degree() > this->degree()) {
            throw std::invalid_argument("The degree of the divisor is greater than the dividend");
        }

        Polynomial<T> dividend = *this;
        std::vector<T> quotient_coeffs(this->degree() - divisor.degree() + 1, 0);
        Polynomial<T> tmp(std::vector<T>(this->degree() + 1, 0));

        while (dividend.degree() >= divisor.degree()) {
            int degree_diff = dividend.degree() - divisor.degree();
            T coeff = dividend.coeffs.back() / divisor.coeffs.back();

            // Update the temporary polynomial for subtraction
            tmp.coeffs[degree_diff] = coeff;

            // Subtract and update the dividend
            for (int i = 0; i <= divisor.degree(); ++i) {
                dividend.coeffs[i + degree_diff] -= coeff * divisor.coeffs[i];
            }
            dividend.coeffs.pop_back();

            // Update the quotient
            quotient_coeffs[degree_diff] = coeff;
        }

        return Polynomial<T>(quotient_coeffs);
    }


    void divide(Polynomial& divisor, Polynomial& quotient, Polynomial& remainder) {
        if (divisor.degree() > this->degree()) {
            throw std::invalid_argument("The degree of the divisor is greater than the dividend");
        }

        Polynomial<T> dividend = *this;
        std::vector<T> quotient_coeffs(this->degree() - divisor.degree() + 1, 0);
        Polynomial<T> tmp;

        while (dividend.degree() >= divisor.degree()) {
            int degree_diff = dividend.degree() - divisor.degree();
            T coeff = dividend.coeffs.back() / divisor.coeffs.back();

            // Update the temporary polynomial for subtraction
            tmp.coeffs = std::vector<T>(degree_diff + 1, 0);
            tmp.coeffs.back() = coeff;

            // Subtract and update the dividend
            for (int i = 0; i <= divisor.degree(); ++i) {
                dividend.coeffs[i + degree_diff] -= coeff * divisor.coeffs[i];
            }
            dividend.coeffs.pop_back();

            // Update the quotient
            quotient_coeffs[degree_diff] = coeff;
        }

        quotient = Polynomial<T>(quotient_coeffs);
        remainder = dividend;
    }

    Polynomial& operator*=(const Polynomial& other){
        *this = *this * other;
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

    std::vector<T> diff(int deg=1){
        std::vector ret = coeffs;
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
            std::cout << "This polynomial object doesnt contain any coefficients.\n";
        }
        else{
            for(int i = _degree; i >= 0; --i){
                std::cout << "(" << coeffs[i] << ")";
                std::cout << "x^" << i;
                if(i != 0) {
                    std::cout << " + ";
                }
            }
            std::cout << ". Degree: " << _degree << '\n';
            std::cout << "Roots: " << roots.size() <<  "; ";
            printVec(roots);
        }
    }       
};
};

#endif