#ifndef FAST_FOURIER_H
#define FAST_FOURIER_H

#include <vector>
#include <iostream>
#include <complex>

namespace PolyGen {

    template<typename T>
    class FourierMultiplication {
    private:
        std::vector<T> roots;
        std::vector<T> coeffs;
        int degree = 0;

    public:
        static void fft(std::vector<T> &a, bool invert) {
            int n = a.size();

            for (int i = 1, j = 0; i < n; i++) {
                int bit = n >> 1;
                for (; j & bit; bit >>= 1)
                    j ^= bit;
                j ^= bit;

                if (i < j)
                    std::swap(a[i], a[j]);
            }

            for (int len = 2; len <= n; len <<= 1) {
                double ang = 2 * M_PI / len * (invert ? -1 : 1);
                T wlen = std::polar(1.0, ang);
                for (int i = 0; i < n; i += len) {
                    T w(1);
                    for (int j = 0; j < len / 2; j++) {
                        T u = a[i + j], v = a[i + j + len / 2] * w;
                        a[i + j] = u + v;
                        a[i + j + len / 2] = u - v;
                        w *= wlen;
                    }
                }
            }

            if (invert) {
                for (T &x : a)
                    x /= n;
            }
        }


        FourierMultiplication() = delete;
    };

}

#endif
