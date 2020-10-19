#include "base.h"
#include <algorithm>

namespace altruct {
namespace math {

/**
 * Slow implementation of a convolution:
 *   r[k] = Sum[f[i] * g[j], k == k_func(i, j)]
 *
 * @param r - arry to store the result; must not be `f` or `g`
 * @param n - size of the arrays
 * @param k_func - grouping function, e.g. `k = i ^ j`
 */
template<typename T, typename F>
void slow_k_convolution(T* r, const T* f, const T* g, int n, F k_func) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int k = k_func(i, j);
            r[k] += f[i] * g[j];
        }
    }
}

template<typename T>
void slow_and_convolution(T* r, const T* f, const T* g, int log_n) {
    slow_k_convolution(r, f, g, 1 << log_n, [](int i, int j){ return i & j; });
}
template<typename T>
void slow_or_convolution(T* r, const T* f, const T* g, int log_n) {
    slow_k_convolution(r, f, g, 1 << log_n, [](int i, int j){ return i | j; });
}
template<typename T>
void slow_xor_convolution(T* r, const T* f, const T* g, int log_n) {
    slow_k_convolution(r, f, g, 1 << log_n, [](int i, int j){ return i ^ j; });
}
template<typename T>
void slow_max_convolution(T* r, const T* f, const T* g, int n) {
    slow_k_convolution(r, f, g, n, [](int i, int j){ return std::max(i, j); });
}
template<typename T>
void slow_cyclic_convolution(T* r, const T* f, const T* g, int n) {
    slow_k_convolution(r, f, g, n, [&](int i, int j){ return (i + j) % n; });
}

/**
 * Fast Radix-2 Decimation-in-Frequency Transform.
 */
template <typename T, typename F>
void fast_radix2_dif_transform(T *f, int log_n, F tr) {
    const int n = 1 << log_n;
    for (int log_m = log_n; log_m >= 1; --log_m) {
        const int m = 1 << log_m, mh = m >> 1;
        for (int i = 0; i < n; i += m) {
            int k1 = i, k2 = i + mh;
            for (int j = 0; j < mh; ++j, ++k1, ++k2) {
                tr(f[k1], f[k2]);
            }
        }
    }
}

/**
 * Fast Walsh�Hadamard Transform.
 */
template<typename T>
void fast_walsh_hadamard_transform(T* f, int log_n) {
    fast_radix2_dif_transform(f, log_n, [](T& u, T& v){
        // (u, v) <-- (u + v, u - v)
        T t = u - v; u += v; v = t;
    });
}

/**
 * Fast Arithmetic Transform (positive sign).
 */
template <typename T>
void fast_arith_transform_plus(T *f, int log_n) {
    fast_radix2_dif_transform(f, log_n, [](T& u, T& v){
        // (u, v) <-- (u, v + u)
        v += u;
    });
}

/**
 * Fast Arithmetic Transform (negative sign).
 */
template <typename T>
void fast_arith_transform_minus(T *f, int log_n) {
    fast_radix2_dif_transform(f, log_n, [](T& u, T& v){
        // (u, v) <-- (u, v - u)
        v -= u;
    });
}

/**
 * And-convolution:
 *   r[k] = Sum[f[i] * g[j], k == i and j]
 *
 * Note: all `f`, `g` and `r` are modified.
 * It is allowed for `f`, `g` and `r` to be the same array.
 *
 * @param log_n - base-2 logarithm of the array length
 */
template<typename T>
void and_convolution(T* r, T* f, T* g, int log_n) {
    const int n = 1 << log_n;
    std::reverse(f, f + n), fast_arith_transform_plus(f, log_n);
    if (g != f) std::reverse(g, g + n), fast_arith_transform_plus(g, log_n);
    for (int k = 0; k < n; ++k) r[k] = f[k] * g[k];
    fast_arith_transform_minus(r, log_n), std::reverse(r, r + n);
}

/**
 * Or-convolution:
 *   r[k] = Sum[f[i] * g[j], k == i or j]
 *
 * Note: all `f`, `g` and `r` are modified.
 * It is allowed for `f`, `g` and `r` to be the same array.
 *
 * @param log_n - base-2 logarithm of the array length
 */
template<typename T>
void or_convolution(T* r, T* f, T* g, int log_n) {
    const int n = 1 << log_n;
    fast_arith_transform_plus(f, log_n);
    if (g != f) fast_arith_transform_plus(g, log_n);
    for (int k = 0; k < n; ++k) r[k] = f[k] * g[k];
    fast_arith_transform_minus(r, log_n);
}

/**
 * Xor-convolution (Dyadic convolution):
 *   r[k] = Sum[f[i] * g[j], k == i xor j]
 *
 * Note: all `f`, `g` and `r` are modified.
 * It is allowed for `f`, `g` and `r` to be the same array.
 *
 * @param log_n - base-2 logarithm of the array length
 */
template<typename T>
void xor_convolution(T* r, T* f, T* g, int log_n) {
    const int n = 1 << log_n;
    fast_walsh_hadamard_transform(f, log_n);
    if (g != f) fast_walsh_hadamard_transform(g, log_n);
    for (int k = 0; k < n; ++k) r[k] = f[k] * g[k];
    fast_walsh_hadamard_transform(r, log_n);
    for (int k = 0; k < n; ++k) r[k] /= n;
}

/**
 * Max-convolution:
 *   r[k] = Sum[f[i] * g[j], k == max(i, j)]
 *
 * Note: all `f`, `g` and `r` are modified.
 * It is allowed for `f`, `g` and `r` to be the same array.
 *
 * @param n - the array length
 */
template<typename T>
void max_convolution(T* r, T* f, T* g, int n) {
    T e0 = zeroT<T>::of(*f);
    T sf = e0, sg = e0;
    for (int k = 0; k < n; ++k) {
        T r_new = f[k] * g[k] + sf * g[k] + sg * f[k];
        sf += f[k], sg += g[k], r[k] = r_new;
    }
}

} // math
} // altruct
