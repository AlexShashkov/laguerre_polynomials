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

	/** \brief Имплементация 'On the Cost of Floating-Point Computation Without Extra-Precise Arithmetic' \n
    Взято из https://github.com/KMBO19MCCO/qdrtcsKahanWepa \n
	\author Prof. W. Kahan

	Статья: https://people.eecs.berkeley.edu/~wkahan/Qdrtcs.pdf \n
	*/
	template <typename number>
    void KahanQuadratic(number a, number b, number c, vector<complex<number>> &roots){
        // x^2, x, c
        //Coefficients should be in CBA order
        try{
            b = b/static_cast<number>(-2);
            if(!std::isfinite(a)) throw std::invalid_argument("Коэффициент при x^2 равен нулю или б/м.");
            if(!std::isfinite(b)) throw std::invalid_argument("Коэффициент при x равен нулю или б/м.");
            number p = b*b;
            number q = a*c;
            //Use the hardware's FMA
            number dp = fma(b,b,-p);
            number dq = fma(a,c,-q);
            // дискриминант
            number d = (p-q) + (dp - dq);
            d = std::max(d,static_cast<number>(0));
            number S = b;
            S = std::sqrt(d)*(sign(S) + (S==0)) + S;
            number Z1 = S/a;
            number Z2 = c/S;
            if(anynotfinite(Z1, Z2)) throw std::invalid_argument("Полученные корни не определены.");

            roots = {Z1, Z2};
        }
        catch(const std::invalid_argument &err){
            std::cerr << "Error occured while working with " << a << " " << b << " " << c << " " << "\n";
            std::cerr << "Invalid argument was passed: " << err.what();
        }
        catch(const std::out_of_range &err){
            std::cerr << "Error occured while working with " << a << " " << b << " " << c << " " << "\n";
            std::cerr << "Out of range: " << err.what();
        }
        catch(...){
            std::cerr << "!!! Error occured while working with " << a << " " << b << " " << c << " " << "\n";
        }
    }

	/** \brief Случай для одного корня \n
	*/
	template <typename number>
    void simpleEquation(number a, number b, vector<complex<number>> &roots){
        // x, c
        number res = b/a;
        if(anynotfinite(a, b, res)){
            std::cerr << "Error occured while working with " << a << " " << b << " " << res << "\n";
            std::cerr << "Invalid argument was passed\n";
            return;
        }
        roots = {res};
    }
}
#endif

