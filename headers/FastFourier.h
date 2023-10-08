#ifndef FAST_FOURIER_H
#define FAST_FOURIER_H

#include <vector>
#include <iostream>

namespace PolyGen {

template<typename T>
class FourierMultiplication{
private:
    std::vector<T> roots;
    std::vector<T> coeffs;
    int degree = 0;
public:
    template<typename T>
    static void fft(vector<T> & a, bool invert) {
        int n = a.size();

        for (int i = 1, j = 0; i < n; i++) {
            int bit = n >> 1;
            for (; j & bit; bit >>= 1)
                j ^= bit;
            j ^= bit;

            if (i < j)
                swap(a[i], a[j]);
        }

        for (int len = 2; len <= n; len <<= 1) {
            double ang = 2 * PI / len * (invert ? -1 : 1);
            T wlen(cos(ang), sin(ang));
            for (int i = 0; i < n; i += len) {
                T w(1);
                for (int j = 0; j < len / 2; j++) {
                    T u = a[i+j], v = a[i+j+len/2] * w;
                    a[i+j] = u + v;
                    a[i+j+len/2] = u - v;
                    w *= wlen;
                }
            }
        }

        if (invert) {
            for (T & x : a)
                x /= n;
        }
    }

    template<typename T>
    int reverse(T num, int lg_n) {
        T res = 0;
        for (int i = 0; i < lg_n; i++) {
            if (num & (1 << i))
                res |= 1 << (lg_n - 1 - i);
        }
        return res;
    }

    static void operator(){

    }

    FourierMultiplication() = delete;
};
}

#endif
