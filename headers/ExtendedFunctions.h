#ifndef EXTENDEDFMA_H
#define EXTENDEDFMA_H

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

	/** \brief Sign
	*/
	template <typename number>
	inline int sign(number val) {
		return (number(0) < val) - (val < number(0));
	}

	/** \brief Fused multiply-substract
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

	/** \brief Fused multiply-substract for complex numbers
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
}
#endif
