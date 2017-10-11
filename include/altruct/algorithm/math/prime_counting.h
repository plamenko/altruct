#pragma once

#include "altruct/algorithm/math/base.h"
#include "altruct/algorithm/math/sums.h"
#include "altruct/structure/container/sqrt_map.h"

namespace altruct {
namespace math {

/**
 * Calculates `PrimePowerSum[z, n / k]` for each `k` in `[1, n]` in `O(n^(5/7))`.
 *
 * Where:
 *  `PrimePowerSum[z, n] := Sum[If[IsPrime[k], k^z, 0], {k, 1, n}]`
 *
 * Note, there is only `O(sqrt n)` different values and the result is given as `sqrt_map`.
 */
template<typename T, typename I>
container::sqrt_map<I, T> prime_power_sum_sqrt(int z, I n, T id) {
    // Initially, we start with the sum of all powers:
    // s[i] = Sum[k^z, {2 <= k <= i}]
    // After each round j, all multiples of a prime p[j] get eliminated:
    // s[i] = Sum[k^z, {2 <= k <= i, SmallestPrimeFactorOf[k] > p[j] || IsPrime[k]}]
    I q = sqrtT(n) + 1;
    container::sqrt_map<I, T> s(q - 1, n);
    for (I i = 1; i < q; i++) {
        s[i] = sum_pow(z, i, id) - id;
    }
    for (I k = n / q; k >= 1; k--) {
        I i = n / k;
        s[i] = sum_pow(z, i, id) - id;
    }
    for (I p = 2; p < q; p++) {
        if (s[p - 1] == s[p]) continue;
        T t = s[p - 1];
        I p2 = sqT(p);
        I k_max = std::min(n / q, n / p2);
        T pz = powT(castOf(id, p), z);
        for (I k = 1; k <= k_max; k++) {
            //I i = n / k;
            //s.el(i) -= (s.el(i / p) - t) * pz;
            I j = n / (k * p);
            s.hi(k) -= (s.el(j) - t) * pz;
        }
        for (I i = q - 1; i >= p2; i--) {
            s.lo(i) -= (s.lo(i / p) - t) * pz;
        }
    }
    return s;
}

/**
 * Calculates `PrimeSum[n / k]` for each `k` in `[1, n]` in `O(n^(5/7))`.
 *
 * Where:
 *  `PrimeSum[n] := Sum[If[IsPrime[k], k, 0], {k, 1, n}]`
 *
 * Note, there is only `O(sqrt n)` different values and the result is given as `sqrt_map`.
 *
 * @param castT - casts I to T
 */
template<typename T, typename I>
container::sqrt_map<I, T> prime_sum_sqrt(I n, T id) {
    return prime_power_sum_sqrt(1, n, id);
}

/**
 * Calculates `PrimeSum[n]` in `O(n^(5/7))`.
 */
template<typename T, typename I>
T prime_sum(I n, T id) {
    return (n < 1) ? zeroT<T>::of(id) : prime_sum_sqrt<T, I>(n, id)[n];
}

/**
 * Calculates `PrimePi[n / k]` for each `k` in `[1, n]` in `O(n^(5/7))`.
 *
 * Where:
 *  `PrimePi[n] := Sum[If[IsPrime[k], 1, 0], {k, 1, n}]`
 *
 * Note, there is only `O(sqrt n)` different values and the result is given as `sqrt_map`.
 */
template<typename I = int64_t>
container::sqrt_map<I, I> prime_pi_sqrt(I n) {
    return prime_power_sum_sqrt(0, n, I(1));
}

/**
 * Calculates `PrimePi[n]` in `O(n^(5/7))`.
 */
template<typename I = int64_t>
I prime_pi(I n) {
    return (n < 1) ? 0 : prime_pi_sqrt<I>(n)[n];
}

/**
 * Calculates `PrimePi1[n / k]` and `PrimePi3[n / k` for each `k` in `[1, n]` in `O(n^(5/7))`.
 *
 * Where:
 *  `PrimePi1[n] := Sum[If[IsPrime[k] && (k % 4) == 1, 1, 0], {k, 1, n}]`
 *  `PrimePi3[n] := Sum[If[IsPrime[k] && (k % 4) == 3, 1, 0], {k, 1, n}]`
 * Note:
 *  `PrimePi[n] = PrimePi1[n] + PrimePi3[n] + [n >= 2]`
 *
 * Note, there is only `O(sqrt n)` different values and the result is given as two `sqrt_map`.
 */
template<typename I = int64_t>
std::vector<container::sqrt_map<I, I>> prime_pi13_sqrt(I n) {
    I q = sqrtT(n) + 1;
    std::vector<container::sqrt_map<I, I>> res{
        container::sqrt_map<I, I>(q, n), // pi1
        container::sqrt_map<I, I>(q, n), // pi3
    };
    if (n <= 1) return res;
    std::vector<I> keys;
    for (I i = 1; i <= q; i++) keys.push_back(n / i);
    for (I k = n / q - 1; k > 0; k--) keys.push_back(k);
    for (I k : keys) {
        res[0][k] = (k - 1) / 4;
        res[1][k] = (k + 1) / 4;
    }
    for (I p = 2; p < q; p++) {
        auto s0 = res[0][p - 1];
        auto s1 = res[1][p - 1];
        if (s0 == res[0][p] && s1 == res[1][p]) continue;
        I p2 = sqT(p);
        for (I k : keys) {
            if (k < p2) break;
            res[p % 4 == 3][k] -= res[0][k / p] - s0;
            res[p % 4 == 1][k] -= res[1][k / p] - s1;
        }
    }
    return res;
}

/**
 * Calculates `PrimePi1[n]` in `O(n^(5/7))`.
 */
template<typename I = int64_t>
I prime_pi1(I n) {
    return (n < 1) ? 0 : prime_pi13_sqrt<I>(n)[0][n];
}

/**
 * Calculates `PrimePi3[n]` in `O(n^(5/7))`.
 */
template<typename I = int64_t>
I prime_pi3(I n) {
    return (n < 1) ? 0 : prime_pi13_sqrt<I>(n)[1][n];
}

} // math
} // altruct
