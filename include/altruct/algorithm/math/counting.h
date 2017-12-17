#pragma once

#include <vector>

#include "altruct/algorithm/math/ranges.h"
#include "altruct/structure/math/polynom.h"
#include "altruct/structure/math/series.h"

namespace altruct {
namespace math {

/**
 * Factorial of `n`. A signle term.
 *
 * Complexity: `O(n)`
 */
template<typename T>
T factorial(int n, T id = T(1)) {
    T f = id;
    for (int i = 2; i <= n; i++) {
        f *= i;
    }
    return f;
}

/**
 * Stirling numbers of the first kind; Table of the first `n` rows and `k` columns.
 *
 * Complexity: `O(n k)`
 */
template<typename T>
std::vector<std::vector<T>> stirling_s1_all(int n, int k, T id = T(1)) {
    auto e0 = zeroT<T>::of(id);
    std::vector<std::vector<T>> s1(n);
    for (int i = 0; i < n; i++) {
        s1[i].resize(std::min(i + 1, k), e0);
        for (int j = 1; j < std::min(i, k); j++) {
            s1[i][j] = s1[i - 1][j - 1] - s1[i - 1][j] * (i - 1);
        }
        if (i < k) s1[i][i] = id;
    }
    return s1;
}

/**
 * Stirling numbers of the first kind; `k`-th column up to `n`.
 *
 * Complexity: `O(n log n log k)`
 */
template<typename T>
std::vector<T> stirling_s1_all_n_for_k(int n, int k, T id = T(1)) {
    typedef seriesX<T> ser;
    auto egf_s1_k = powT(ser{ id, id }.resize(n).ln(), k) / factorial(k, id);
    auto s1_k = egf_s1_k.make_ordinary();
    s1_k.p.resize(n);
    return s1_k.p.c;
}

/**
 * Stirling numbers of the first kind; `n`-th row.
 *
 * Complexity: `O(n log^2 n)`
 */
template<typename T>
std::vector<T> stirling_s1_all_k_for_n(int n, T id = T(1)) {
    if (n == 0) return{ id };
    typedef polynom<T> poly;
    std::vector<poly> vp(n);
    for (int i = 0; i < n; i++) vp[i] = {-id * i, id};
    for (int l = 1; l < n; l *= 2) {
        for (int i = 0; i + l < n; i += l * 2) {
            vp[i] *= vp[i + l];
            vp[i + l] = {};
        }
    }
    return vp[0].c;
}

/**
 * Stirling number of the first kind; A single term.
 *
 * Complexity: `O(n log n log k)`; TODO: is there a better way?
 */
template<typename T>
T stirling_s1(int n, int k, T id = T(1)) {
    T e0 = zeroT<T>::of(id);
    if (k == n) return id;
    if (k <= 0 || k > n) return e0;
    return stirling_s1_all_k_for_n(n, id)[k];
}

/**
 * Stirling numbers of the second kind; Table of the first `n` rows and `k` columns.
 *
 * Complexity: `O(n k)`
 */
template<typename T>
std::vector<std::vector<T>> stirling_s2_all(int n, int k, T id = T(1)) {
    auto e0 = zeroT<T>::of(id);
    std::vector<std::vector<T>> s1(n);
    for (int i = 0; i < n; i++) {
        s1[i].resize(std::min(i + 1, k), e0);
        for (int j = 1; j < std::min(i, k); j++) {
            s1[i][j] = s1[i - 1][j - 1] + s1[i - 1][j] * j;
        }
        if (i < k) s1[i][i] = id;
    }
    return s1;
}

/**
 * Stirling numbers of the second kind; `k`-th column up to `n`.
 *
 * Complexity: `O(n log n log k)`
 */
template<typename T>
std::vector<T> stirling_s2_all_n_for_k(int n, int k, T id = T(1)) {
    typedef seriesX<T> ser;
    auto egf_s2_k = powT(ser::exp(id, n) - ser{ id }, k) / factorial(k, id);
    auto s2_k = egf_s2_k.make_ordinary();
    s2_k.p.resize(n);
    return s2_k.p.c;
}

/**
 * Stirling numbers of the second kind; `n`-th row.
 *
 * Complexity: `O(n log n)`
 */
template<typename T>
std::vector<T> stirling_s2_all_k_for_n(int n, T id = T(1)) {
    typedef seriesX<T> ser;
    auto v = range(n + 1, id); power(v.begin(), v.end(), n);
    auto s2_n = ser::exp(-id, n + 1) * ser(v, n + 1).make_exponential();
    s2_n.p.resize(n + 1);
    return s2_n.p.c;
}

/**
 * Stirling number of the second kind. A signle term.
 *
 * Complexity: `O(k log n)`
 */
template<typename T, typename I>
T stirling_s2(I n, int k, T id = T(1)) {
    T e0 = zeroT<T>::of(id);
    if (k == n) return id;
    if (k <= 0 || k > n) return e0;
    T b = id, f = id, r = e0;
    for (int j = 1; j <= k; j++) {
        b *= (k + 1 - j), b /= j, f *= j; // TODO: `j` inverse may add another O(log k) factor! Precompute?
        T t = b * powT(id * j, n);
        // b = -b;
        r += ((k - j) % 2 == 1) ? -t : t;
    }
    return r / f;
}

/**
 * Series expansion of the Euler function up to `n`.
 *
 *  Euler[x] = Product[1 - x^k, {k, 1, inf}]
 *
 * Complexity: `O(n)`
 */
template<typename T>
std::vector<T> euler_function(int n, T id = T(1)) {
    // See https://en.wikipedia.org/wiki/Pentagonal_number_theorem
    T e0 = zeroT<T>::of(id);
    std::vector<T> e(n, e0);
    e[0] = id;
    for (int k = 4, m = 1; m < n; m += k, k += 3) {
        e[m] = (k % 2 == 0) ? -id : id;
    }
    for (int k = 5, m = 2; m < n; m += k, k += 3) {
        e[m] = (k % 2 == 0) ? id : -id;
    }
    return e;
}

/**
 * Partition numbers up to `n`.
 *
 *   p[n] = number of partitions of n
 *
 * Complexity: `O(n sqrt n)`
 * Note: this can be calculated in `O(M(n))` by `partitions_p`
 */
template<typename T>
std::vector<T> partitions_p_slow(int n, T id = T(1)) {
    // See https://en.wikipedia.org/wiki/Pentagonal_number_theorem
    T e0 = zeroT<T>::of(id);
    std::vector<T> p(n, e0);
    p[0] = id;
    for (int i = 1; i < n; i++) {
        for (int k = 4, m = 1; m <= i; m += k, k += 3) {
            p[i] += (k % 2 == 0) ? p[i - m] : -p[i - m];
        }
        for (int k = 5, m = 2; m <= i; m += k, k += 3) {
            p[i] += (k % 2 == 0) ? -p[i - m] : p[i - m];
        }
    }
    return p;
}

/**
 * Distinct partition numbers up to `n`.
 *
 *   q[n] = number of partitions of n into distinct parts
 *
 * Complexity: `O(n sqrt n)`
 * Note: this can be calculated in `O(M(n))` by `partitions_p_distinct`
 */
template<typename T>
std::vector<T> partitions_p_distinct_slow(int n, T id = T(1)) {
    // See https://oeis.org/A000009 Jerome Malenfant, Feb 16 2011
    auto q = partitions_p_slow(n, id);
    for (int i = n - 1; i >= 0; i--) {
        for (int k = 8, m = 2; m <= i; m += k, k += 6) {
            q[i] += (k % 4 == 0) ? -q[i - m] : q[i - m];
        }
        for (int k = 10, m = 4; m <= i; m += k, k += 6) {
            q[i] += (k % 4 == 0) ? q[i - m] : -q[i - m];
        }
    }
    return q;
}

/**
 * Partition numbers up to `n`.
 *
 *   p[n] = number of partitions of n
 *
 * Complexity: `O(n log n)`
 */
template<typename T>
std::vector<T> partitions_p(int n, T id = T(1)) {
    auto e = seriesX<T>(euler_function(n, id));
    return e.inverse().p.c;
}

/**
 * Distinct partition numbers up to `n`.
 *
 *   q[n] = number of partitions of n into distinct parts
 *
 * Complexity: `O(n log n)`
 */
template<typename T>
std::vector<T> partitions_p_distinct(int n, T id = T(1)) {
    auto e = seriesX<T>(euler_function(n, id));
    return (e.sub_pow(2) / e).p.c;
}

/**
 * Distinct odd partition numbers up to `n`.
 *
 *   q1[n] = number of partitions of n into distinct odd parts
 *
 * Complexity: `O(n log n)`
 */
template<typename T>
std::vector<T> partitions_p_distinct_odd(int n, T id = T(1)) {
    auto e = seriesX<T>(euler_function(n, id));
    return (e.sub_mul(-id) / e.sub_pow(2)).p.c;
}

} // math
} // altruct
