#include "base.h"
#include "bits.h"
#include <algorithm>
#include <vector>

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

/**
 * Subset-Sum:
 *   r[k] = Sum[f[i], i is bit-subset of k]
 *
 * Complexity: O(3^log_n) = O(n^1.585)
 * `f` and `r` must not be the same array.
 *
 * @param log_n - base-2 logarithm of the array length
 */
template<typename T>
void slow_subset_sum(T* r, T* f, int log_n) {
    const uint32_t n = uint32_t(1) << log_n;
    for (uint32_t w = 0; w < n; w++) {
        r[w] = zeroOf(f[0]);
        for (uint32_t sub = w; ; sub = (sub - 1) & w) {
            r[w] += f[sub];
            if (sub == 0) break;
        }
    }
}

/**
 * Subset-Sum (Moebius Transform on the Subset Lattice):
 *   r[k] = Sum[f[i], i is bit-subset of k]
 *
 * Complexity: O(log_n 2^log_n) = O(n log n)
 * It is allowed for `f` and `r` to be the same array.
 *
 * @param log_n - base-2 logarithm of the array length
 */
template<typename T>
void fast_subset_sum(T* r, T* f, int log_n) {
    const uint32_t n = uint32_t(1) << log_n;
    if (r != f) {
        for (uint32_t w = 0; w < n; w++) {
            r[w] = f[w];
        }
    }
    for (int i = 0; i < log_n; i++) {
        for (uint32_t w = 0; w < n; w++) {
            if (w & (uint32_t(1) << i)) {
                r[w] += r[w ^ (uint32_t(1) << i)];
            }
        }
    }
}

/**
 * Subset-Sum Inverse (Moebius Inversion on the Subset Lattice):
 *   f[k] = Sum[r[i], i is bit-subset of k]
 *
 * Complexity: O(log_n 2^log_n) = O(n log n)
 * It is allowed for `f` and `r` to be the same array.
 *
 * @param log_n - base-2 logarithm of the array length
 */
template<typename T>
void fast_subset_sum_inverse(T* r, T* f, int log_n) {
    const uint32_t n = uint32_t(1) << log_n;
    if (r != f) {
        for (uint32_t w = 0; w < n; w++) {
            r[w] = f[w];
        }
    }
    for (int i = 0; i < log_n; i++) {
        for (uint32_t w = 0; w < n; w++) {
            if (w & (uint32_t(1) << i)) {
                r[w] -= r[w ^ (uint32_t(1) << i)];
            }
        }
    }
}

/**
 * Subset-Sum Convolution (Moebius Transform on the Subset Lattice):
 *   r[k] = Sum[f[i] * g[k^i], i is bit-subset of k]
 *
 * Complexity: O(log_n^2 2^log_n) = O(n log^2 n)
 * `f` and `g` must not be the same array as `r`.
 *
 * @param log_n - base-2 logarithm of the array length
 */
template<typename T>
void slow_subset_convolution(T* r, T* f, T* g, int log_n) {
    const uint32_t n = uint32_t(1) << log_n;
    for (uint32_t w = 0; w < n; w++) {
        r[w] = zeroOf(f[0]);
        for (uint32_t sub = w; ; sub = (sub - 1) & w) {
            r[w] += f[sub] * g[w ^ sub];
            if (sub == 0) break;
        }
    }
}

/**
 * Subset-Sum Convolution (Uses Moebius Transform on the Subset Lattice):
 *   r[k] = Sum[f[i] * g[k^i], i is bit-subset of k]
 *
 * Complexity: O(log_n^2 2^log_n) = O(n log^2 n)
 * It is allowed for `f` and `r` to be the same array.
 *
 * @param log_n - base-2 logarithm of the array length
 */
template<typename T>
void fast_subset_convolution(T* r, T* f, T* g, int log_n) {
    const uint32_t n = uint32_t(1) << log_n;
    // n log_n memory allocation for ranked moebius transform and inversion
    std::vector<std::vector<T>> f1(log_n + 1, std::vector<T>(n, zeroOf(f[0]))), g1 = f1, h1 = f1;
    for (uint32_t w = 0; w < n; w++) {
        f1[bit_cnt1(w)][w] = f[w];
        g1[bit_cnt1(w)][w] = g[w];
    }
    for (int i = 0; i <= log_n; i++) {
        fast_subset_sum(f1[i].data(), f1[i].data(), log_n);
        fast_subset_sum(g1[i].data(), g1[i].data(), log_n);
    }
    for (int k = 0; k <= log_n; k++) {
        for (int j = 0; j <= k; j++) {
            for (uint32_t w = 0; w < n; w++) {
                h1[k][w] += f1[j][w] * g1[k - j][w];
            }
        }
    }
    for (int i = 0; i <= log_n; i++) {
        fast_subset_sum_inverse(h1[i].data(), h1[i].data(), log_n);
    }
    for (uint32_t w = 0; w < n; w++) {
        r[w] = h1[bit_cnt1(w)][w];
    }
}

} // math
} // altruct
