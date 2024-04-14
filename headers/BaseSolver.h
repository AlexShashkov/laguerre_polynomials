#ifndef BASESOLVER_H
#define BASESOLVER_H

#include<vector>

template<typename T>
class BaseSolver {
protected:
    static constexpr T eps   = std::numeric_limits<T>::epsilon();
    static constexpr T big   = std::numeric_limits<T>::max();
    static constexpr T small = std::numeric_limits<T>::min();
public:
    /**
     * \brief Find the roots of a polynomial.
         * \param coeff Polynomials coefficients.
         * \param roots Vector to store the roots.
         * \param conv Vector to store convergence status of each root.
         * \param itmax Maximum number of iterations.
     */
    virtual void operator()(const std::vector<T>& coeff, std::vector<std::complex<T>>& roots, std::vector<int>& conv, int itmax) = 0;
    virtual ~BaseSolver() {}
};

#endif