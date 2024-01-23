#ifndef EXTENDEDFMA_H
#define EXTENDEDFMA_H


#ifndef EXPONENT
    #define EXPONENT 4 // Exponent for big NUMBER notation
#endif
#ifndef MANTISSA
    #define MANTISSA 32 // Mantissa for big NUMBER notation
#endif

#include "ttmath/ttmath.h" // Bignum C++ library by Tomasz Sowa
#include <cmath>
#include <numbers>
#include <vector>
#include <complex>
#include <iostream>
#include <algorithm>
#include <iomanip>

using namespace std::complex_literals;

namespace Laguerre{
	using std::complex;
	using std::vector;
	using std::fma;

	/** \brief Check if at least one number is not finite
	*/
	template<typename ... T>
	inline bool anynotfinite(T && ... t)
	{
		return ((!std::isfinite(t)) || ...);
	}

	/** \brief Check if a complex number contains NaN or Inf values.
	*/
	template <typename T>
    bool complexnotfinite(complex<T> a, T big){
        T re_a = real(a);
        T im_a = imag(a);
        return !std::isfinite(re_a) || !std::isfinite(im_a) || (fabs(re_a) > big) || (fabs(im_a) > big);
    }

	/** \brief Check if at least one complex number has imaginary part
	*/
	template<typename ... T>
	inline bool anycomplex(T && ... t)
	{
		return ((t.imag() > 0) || ...);
	}

	/** \brief Check if at least one number in complex vector has imaginary part
	*/
	template<typename T>
	inline bool anycomplex(vector<complex<T>> vec)
	{
		return std::any_of(vec.begin(), vec.end(), [](complex<T> t){return t.imag() > 0;});
	}


	/** \brief Sign
	*/
	template <typename number>
	inline int sign(number val) {
		return (number(0) < val) - (val < number(0));
	}

	/** \brief Fused multiply-substract a*b+(-d*c)+d*c+(-d*c)
	*/
	template <typename number>
	inline number fms(number a, number b, number c, number d) {
		auto tmp = -d * c;
		return fma(a, b, tmp) + fma(d, c, tmp);
	}

	/** \brief Fused multiply-substract for complex numbers
	*/
	template <typename number>
	inline complex<number> fms(std::complex<number> a, std::complex<number> b, std::complex<number> c,std::complex<number> d) {
		return {
			fms(a.real(),b.real(),c.real(),d.real())+fms(c.imag(),d.imag(),a.imag(),b.imag()),
			fms(a.real(),b.imag(),c.real(),d.imag())+fms(b.real(),a.imag(),d.real(),c.imag())
		};
	}

	/** \brief Fused multiply-add for complex numbers
	 * 	\return a*b+c
	*/
	template<typename number>
	inline complex<number> fma(complex<number> a, complex<number> b, complex<number> c) {
		return {
			fms(a.real(), b.real(), a.imag(), b.imag()) + c.real(),
			fma(a.real(), b.imag(), fma(a.imag(), b.real(), c.imag()))
		};
	}

	/** \brief Fused multiply-substract for complex&real numbers
	 * 	\return a*b+c
	*/
	template <typename number>
	inline complex<number> fma(std::complex<number> a, std::complex<number> b, number c){
		return {std::fma(a.real(),b.real(),std::fma(-b.imag(),a.imag(),c)),fms(a.real(),b.imag(),-b.real(),a.imag())};
	}

	/** \brief Fused multiply-substract for complex&real numbers
	 * 	\return a*b+c
	*/
	template <typename number>
	inline complex<number> fma(complex<number> a, number b, complex<number> c){
		return {std::fma(a.real(),b,c.real()),std::fma(a.imag(),b,c.imag())};
	}

	/** \brief Fused multiply-substract for complex&real numbers
	 * 	\return a*b+c
	*/
	template <typename number>
	inline complex<number> fma(number a,number b, complex<number> c){
		return {std::fma(a,b,c.real()),c.imag()};
	}

	/** \brief Fused multiply-substract for complex&real numbers
	 * 	\return a*b+c
	*/
	template <typename number>
	inline complex<number> fma(complex<number> a,number b, number c){
		return {std::fma(a.real(),b,c),a.imag()*b};
	}

	/** \brief Print out vector
	*/
	template <typename number>
	inline void printVec(vector<number> vec){
        for (number i: vec) std::cout << i << " ";
        std::cout << "\n";
    } 

	/** \brief Convert vector<T> to vector<number>
     * \return new vector of template type <number>
	*/
	template <typename number, typename T>
	inline vector<number> castVec(vector<T> vec){
        vector<number> res;
        for (T i: vec) res.push_back(static_cast<number>(i));
        return res;
    } 

	/** \brief Diff vector
     * \return new vector of template type <T>
	*/
	template <typename T>
	inline vector<T> diff(vector<T> coeffs, int deg=1){
        std::vector ret = coeffs;
        for(int j = 0; j < deg; ++j){
            for (int i = 1; i < ret.size(); ++i) {
                ret[i] *= static_cast<T>(i);
            }
            if (ret.size()) {
                ret.erase(ret.begin());
            }
            else break;
        }
        return ret;
    } 

	/** \brief Divide vectors (coefficients of polynomials)
	*/
	template <typename T>
	inline void divide(vector<T> dividend, vector<T> divisor,
					vector<T>& quotient, vector<T>& remainder) {
		if (divisor.size() > dividend.size()) {
            throw std::invalid_argument("The degree of the divisor is greater than the dividend");
        }

		std::vector<T> quotient_coeffs(dividend.size() - divisor.size() + 1, 0);
        std::vector<T> tmp;

        while (dividend.size() >= divisor.size()) {
            int degree_diff = dividend.size() - divisor.size();
            T coeff = dividend.back() / divisor.back();

            // Update the temporary polynomial for subtraction
            tmp = std::vector<T>(degree_diff + 1, 0);
            tmp.back() = coeff;

            // Subtract and update the dividend
            for (int i = 0; i <= divisor.size()-1; ++i) {
                dividend[i + degree_diff] -= coeff * divisor[i];
            }
            dividend.pop_back();

            // Update the quotient
            quotient_coeffs[degree_diff] = coeff;
        }

        quotient = quotient_coeffs;
        remainder = dividend;
	} 

	/** \brief Divide vectors (coefficients of polynomials)
	*/
	template <typename T>
	inline void getRemainder(vector<T> dividend, vector<T> divisor, vector<T>& remainder) {
		if (divisor.size() > dividend.size()) {
            throw std::invalid_argument("The degree of the divisor is greater than the dividend");
        }

        std::vector<T> tmp;
        while (dividend.size() >= divisor.size()) {
            int degree_diff = dividend.size() - divisor.size();
            T coeff = dividend.back() / divisor.back();

            // Update the temporary polynomial for subtraction
            tmp = std::vector<T>(degree_diff + 1, 0);
            tmp.back() = coeff;

            // Subtract and update the dividend
            for (int i = 0; i <= divisor.size()-1; ++i) {
                dividend[i + degree_diff] -= coeff * divisor[i];
            }
            dividend.pop_back();
        }
        remainder = dividend;
	} 

	/** \brief Divide vectors (coefficients of polynomials)
	*/
	template <typename T>
	inline void getRemainderFromBigNum(vector<T> dividend, vector<T> divisor, vector<T>& remainder) {
		if (divisor.size() > dividend.size()) {
            throw std::invalid_argument("The degree of the divisor is greater than the dividend");
        }

		std::vector<ttmath::Big<EXPONENT, MANTISSA>> big_dividend, big_divisor, tmp;
		for (const T num : dividend) {
            ttmath::Big<EXPONENT, MANTISSA> bignum(num);
            big_dividend.push_back(bignum);
        }
		for (const T num : divisor) {
            ttmath::Big<EXPONENT, MANTISSA> bignum(num);
            big_divisor.push_back(bignum);
        }

        while (big_dividend.size() >= big_divisor.size()) {
            int degree_diff = big_dividend.size() - big_divisor.size();
            ttmath::Big<EXPONENT, MANTISSA> coeff = big_dividend.back() / big_divisor.back();

            // Update the temporary polynomial for subtraction
            tmp = std::vector<ttmath::Big<EXPONENT, MANTISSA>>(degree_diff + 1, ttmath::Big<EXPONENT, MANTISSA>(0.0));
            tmp.back() = coeff;

            // Subtract and update the dividend
            for (int i = 0; i <= big_divisor.size()-1; ++i) {
                big_dividend[i + degree_diff] -= coeff * big_divisor[i];
            }
            big_dividend.pop_back();
        }
		remainder = std::vector<T>(big_dividend.size(), static_cast<T>(0));
		for (int i=0; i < big_dividend.size(); ++i) {
            remainder[i] = static_cast<T>(big_dividend[i].ToDouble());
        }
	} 
}
#endif
