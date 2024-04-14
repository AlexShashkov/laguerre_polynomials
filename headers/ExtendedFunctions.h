// Александр, Дмитрий

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
	using std::abs;
	using std::fma;
	using std::isfinite;
	using std::any_of;
	using std::invalid_argument;

	/** \brief Check if at least one number in vector is not finite
	*/
	template<typename T>
	inline bool anynotfinite(const vector<T>& vec)
	{
		return any_of(vec.begin(), vec.end(), [](T t){return !isfinite(t);});
	}

	/** \brief Check if at least one number is not finite
	*/
	template<typename ... T>
	inline bool anynotfinite(T && ... t)
	{
		return ((!isfinite(t)) || ...);
	}

	/** \brief Check if a complex number contains NaN or Inf values.
	*/
	template <typename T>
    bool complexnotfinite(const complex<T>& a, const T& big){
        T re_a = real(a);
        T im_a = imag(a);
        return !isfinite(re_a) || !isfinite(im_a) || (fabs(re_a) > big) || (fabs(im_a) > big);
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
	inline bool anycomplex(const vector<complex<T>>& vec)
	{
		return any_of(vec.begin(), vec.end(), [](complex<T> t){return t.imag() > 0;});
	}


	/** \brief Sign
	*/
	template <typename number>
	inline int sign(const number& val) {
		return (number(0) < val) - (val < number(0));
	}

	/** \brief Fused multiply-substract a*b+(-d*c)+d*c+(-d*c)
	*/
	template <typename number>
	inline number fms(const number& a, const number& b, const number& c, const number& d) {
		auto tmp = -d * c;
		return fma(a, b, tmp) + fma(d, c, tmp);
	}

	/** \brief Fused multiply-substract for complex numbers
	*/
	template <typename number>
	inline complex<number> fms(const complex<number>& a, const complex<number>& b,
		const complex<number>& c, const complex<number>& d) {
		return {
			fms(a.real(),b.real(),c.real(),d.real())+fms(c.imag(),d.imag(),a.imag(),b.imag()),
			fms(a.real(),b.imag(),c.real(),d.imag())+fms(b.real(),a.imag(),d.real(),c.imag())
		};
	}

	/** \brief Fused multiply-add for complex numbers
	 * 	\return a*b+c
	*/
	template<typename number>
	inline complex<number> fma(const complex<number>& a, const complex<number>& b, 
		const complex<number>& c) {
		return {
			fms(a.real(), b.real(), a.imag(), b.imag()) + c.real(),
			fma(a.real(), b.imag(), fma(a.imag(), b.real(), c.imag()))
		};
	}

	/** \brief Fused multiply-substract for complex&real numbers
	 * 	\return a*b+c
	*/
	template <typename number>
	inline complex<number> fma(const complex<number>& a, const complex<number>& b, const number& c){
		return {fma(a.real(),b.real(),fma(-b.imag(),a.imag(),c)),fms(a.real(),b.imag(),-b.real(),a.imag())};
	}

	/** \brief Fused multiply-substract for complex&real numbers
	 * 	\return a*b+c
	*/
	template <typename number>
	inline complex<number> fma(const complex<number>& a, const number& b, const complex<number>& c){
		return {fma(a.real(),b,c.real()),fma(a.imag(),b,c.imag())};
	}

	/** \brief Fused multiply-substract for complex&real numbers
	 * 	\return a*b+c
	*/
	template <typename number>
	inline complex<number> fma(const number& a, const number& b, const complex<number>& c){
		return {fma(a,b,c.real()),c.imag()};
	}

	/** \brief Fused multiply-substract for complex&real numbers
	 * 	\return a*b+c
	*/
	template <typename number>
	inline complex<number> fma(const complex<number>& a, const number& b, const number& c){
		return {fma(a.real(),b,c),a.imag()*b};
	}

	/** \brief Print out vector
	*/
	template <typename number>
	inline void printVec(const vector<number>& vec){
        for (number i: vec) std::cout << i << " ";
        std::cout << "\n";
    } 

	/** \brief Convert vector<T> to vector<number>
     * \return new vector of template type <number>
	*/
	template <typename number, typename T>
	inline vector<number> castVec(const vector<T>& vec){
        vector<number> res;
        for (T i: vec) res.push_back(static_cast<number>(i));
        return res;
    } 

	/** \brief Diff vector
     * \return new vector of template type <T>
	 * 
	 * 
	*/
	template <typename T>
	inline vector<T> diff(const vector<T>& coeffs, int deg=1){
        vector ret = coeffs;
        for (int j = 0; j < deg && ret.size(); ++j) {
			for (int i = 1; i < ret.size(); ++i) {
				ret[i] *= static_cast<T>(i);
			}
			ret.erase(ret.begin());
		}
        return ret;
    } 

	/** \brief Divide vectors (coefficients of polynomials)
	*/
	template <typename T>
	inline void divide(vector<T> dividend, const vector<T>& divisor,
					vector<T>& quotient, vector<T>& remainder) {
		if (divisor.size() > dividend.size()) {
            throw invalid_argument("The degree of the divisor is greater than the dividend");
        }

		vector<T> quotient_coeffs(dividend.size() - divisor.size() + 1, 0);

		int degree_diff;
		T coeff;

        while (dividend.size() >= divisor.size()) {
            degree_diff = dividend.size() - divisor.size();
            coeff = dividend.back() / divisor.back();

			if(anynotfinite(coeff)){
                throw invalid_argument("Highest degree of coefficient for divisor is probably zero: " << divisor.back());
            }

            // Subtract and update the dividend
            for (int i = 0; i <= divisor.size()-1; ++i) {
                // dividend[i + degree_diff] -= coeff * divisor[i];
                dividend[i + degree_diff] = fma(-coeff, divisor[i], dividend[i + degree_diff]);
            }
            dividend.pop_back();

            // Update the quotient
            quotient_coeffs[degree_diff] = coeff;
        }

        quotient = quotient_coeffs;
        remainder = dividend;
	} 

	/** \brief Get remainder from vectors division (coefficients of polynomials)
	*/
	template <typename T>
	inline void getRemainder(vector<T> dividend, const vector<T>& divisor, vector<T>& remainder) {
		if (divisor.size() > dividend.size()) {
            throw invalid_argument("The degree of the divisor is greater than the dividend");
        }

		int degree_diff;
		T coeff;

        while (dividend.size() >= divisor.size()) {
            degree_diff = dividend.size() - divisor.size();
            coeff = dividend.back() / divisor.back();

			if(anynotfinite(coeff)){
                throw invalid_argument("Highest degree of coefficient for divisor is probably zero: " << divisor.back());
            }

            // Subtract and update the dividend
            for (int i = 0; i <= divisor.size()-1; ++i) {
                // dividend[i + degree_diff] -= coeff * divisor[i];
				dividend[i + degree_diff] = fma(-coeff, divisor[i], dividend[i + degree_diff]);
            }
            dividend.pop_back();
        }
        remainder = dividend;
	} 

	/** \brief Get remainder from vectors division, but with bignum precision (coefficients of polynomials)
	*/
	template <typename T>
	inline void getRemainderFromBigNum(const vector<T>& dividend, const vector<T>& divisor, vector<T>& remainder) {
		if (divisor.size() > dividend.size()) {
            throw invalid_argument("The degree of the divisor is greater than the dividend");
        }

		vector<ttmath::Big<EXPONENT, MANTISSA>> big_dividend, big_divisor;
		for (const T num : dividend) {
            big_dividend.push_back(ttmath::Big<EXPONENT, MANTISSA>(num));
        }
		for (const T num : divisor) {
            big_divisor.push_back(ttmath::Big<EXPONENT, MANTISSA>(num));
        }

		ttmath::Big<EXPONENT, MANTISSA> coeff;
		if (big_divisor.back() == ttmath::Big<EXPONENT, MANTISSA>(0)) {
			throw invalid_argument("Highest degree of coefficient for divisor is probably zero.");
		}
        while (big_dividend.size() >= big_divisor.size()) {
            coeff = big_dividend.back() / big_divisor.back();

            // Subtract and update the dividend
            for (int i = 0; i <= big_divisor.size()-1; ++i) {
				// bignum doesnt support fma
                big_dividend[i + big_dividend.size() - big_divisor.size()] -= coeff * big_divisor[i];
            }
            big_dividend.pop_back();
        }
		remainder = vector<T>(big_dividend.size(), static_cast<T>(0));
		for (int i=0; i < big_dividend.size(); ++i) {
            remainder[i] = static_cast<T>(big_dividend[i].ToDouble());
        }
	} 

	/** \brief Multiply vectors (coefficients of polynomials)
	*/
	template <typename T>
	inline vector<T> multiply(const vector<T>&  polynom_A, const vector<T>&  polynom_B, int stop_on = 0, bool is_bn=false){
		vector<T> polynomial(polynom_A.size() + polynom_B.size() - 1 - stop_on, T(0));
		for(size_t i = stop_on; i < polynom_A.size(); ++i)
		{
			for(size_t j = 0; j < polynom_B.size(); ++j)
			{
				// polynomial[i + j - stop_on] += polynom_A[i] * polynom_B[j];
				polynomial[i + j - stop_on] = fma(polynom_A[i], polynom_B[j], polynomial[i + j - stop_on]);
				// std::cout << polynom_A[i] << " " << polynom_B[j] << " " << polynomial[i + j - stop_on].ToDouble() << "\n";
			}
		}
		return  polynomial;
	}

	/** \brief Multiply vectors of bignum (coefficients of polynomials), as it doesnt support fma
	*/
	template <typename T>
	inline vector<T> multiply_bn(const vector<T>&  polynom_A, const vector<T>&  polynom_B, int stop_on = 0){
		vector<T> polynomial(polynom_A.size() + polynom_B.size() - 1 - stop_on, T(0));
		for(size_t i = stop_on; i < polynom_A.size(); ++i)
		{
			for(size_t j = 0; j < polynom_B.size(); ++j)
			{
				polynomial[i + j - stop_on] += polynom_A[i] * polynom_B[j];
				// std::cout << polynom_A[i] << " " << polynom_B[j] << " " << polynomial[i + j - stop_on].ToDouble() << "\n";
			}
		}
		return  polynomial;
	}
}
#endif
