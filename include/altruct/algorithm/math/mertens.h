#pragma once

#include "altruct/algorithm/math/divisor_sums.h"

namespace altruct {
namespace math {

/**
 * Sieves Mertens up to `n` in `O(n log log n)`.
 * Prefix sum of MoebiusMu: `Sum[MoebiusMu[k], {k, 1, n-1}]`
 *
 * @param M - table to store the calculated values; accessed via [] operator
 * @param n - bound up to which to sieve (exclusive)
 * @param pa - table of all `m` prime numbers up to `n`
 */
template<typename T, typename TBL>
void sieve_mertens(TBL& M, int n, int* pa, int m, T id = T(1)) {
    auto one = [&](int k){ return id; };
    sieve_m_multiplicative(M, one, one, n, pa, m);
}

/**
 * Sieves MertensOdd up to `n` in `O(n log log n)`.
 * Sums only odd terms: `Sum[MoebiusMu[k], {k, 1, n-1, 2}]`
 *
 * @param M - table to store the calculated values; accessed via [] operator
 * @param n - bound up to which to sieve (exclusive)
 * @param pa - table of all `m` prime numbers up to `n`
 */
template<typename T, typename TBL>
void sieve_mertens_odd(TBL& M1, int n, int* pa, int m, T id = T(1)) {
    auto zero = zeroOf(id);
    auto t = [&](int k) { return id; };
    auto p = [&](int k) { return (k % 2 == 1) ? id : zero; };
    sieve_m_multiplicative(M1, t, p, n, pa, m);
}

/**
 * Sieves MertensEven up to `n` in `O(n log n)`.
 * Sums only even terms: `Sum[MoebiusMu[k], {k, 2, n-1, 2}]`
 *
 * @param M - table to store the calculated values; accessed via [] operator
 * @param n - bound up to which to sieve (exclusive)
 * @param pa - table of all `m` prime numbers up to `n`
 */
template<typename T, typename TBL>
void sieve_mertens_even(TBL& M0, int n, T id = T(1)) {
    auto zero = zeroOf(id);
    auto t = [&](int k) { return (k > 1) ? -id : zero; };
    auto p = [&](int k) { return (k % 2 == 1) ? id : zero; };
    sieve_m(M0, t, p, n); // not multiplicative!
}

/**
 * Sieves MertensEven and MertensOdd up to `n` in `O(n log log n)`.
 *
 * @param M - table to store the calculated values; accessed via [] operator
 * @param n - bound up to which to sieve (exclusive)
 * @param pa - table of all `m` prime numbers up to `n`
 */
template<typename T, typename TBL>
void sieve_mertens_even_odd(TBL& M0, TBL& M1, int n, int* pa, int m, T id = T(1)) {
    sieve_mertens_odd(M1, n, pa, m, id);
    sieve_mertens(M0, n, pa, m, id); // M0 = M - M1
    for (int k = 0; k < n; k++) M0[k] -= M1[k];
}

/**
 * Mertens function: `Sum[MoebiusMu(k), {k, 1, n}]` in `O(n^(3/4))` or `O(n^(2/3))`.
 *
 * Note, to achieve the better `O(n^(2/3))` complexity, `tbl` values
 * smaller than `O(n^(2/3))` have to be precomputed with sieve in advance.
 *
 * @deprecated - `mertens` implementation that uses `pi_tbl` is faster in practice
 *               because this implementation requires too much memory
 *
 * @param n - argument at which to evaluate `M`
 * @param tbl - table to store the calculated values
 */
template<typename T, typename I, typename TBL>
T mertens(I n, TBL& tbl, T id = T(1)) {
    // p = 1, f = mu, g = delta, t = 1
    auto one = [&](I k){ return id; };
    return sum_m<T>(one, n, tbl, id);
}
/**
 * Mertens function in `O(n^(2/3))`
 *
 * @param n - argument at which to evaluate `M`
 * @param pi_tbl - table of prime_pi(n/i) for all i
 * @param pa - table of all `psz` prime numbers up to `sqrt(n)`, or `sqrt(n log n)` (it's faster)`
 */
template<typename T, typename PI_TBL>
altruct::container::sqrt_map<int64_t, T> mertens(int64_t n, PI_TBL& pi_tbl, const int* pa, int psz, T id = T(1)) {
    auto mu = [&](T f_pe1, int p, int e) { return castOf(id, (e > 1) ? 0 : -1); };
    auto s1 = [&](int64_t n) { return -castOf(id, pi_tbl[n]); };
    return sum_multiplicative<T>(s1, mu, n, pa, psz, id);
}

/**
 * MertensOdd function: `Sum[MoebiusMu(k), {k, 1, n, 2}]` in `O(n^(3/4))` or `O(n^(2/3))`.
 *
 * Note, to achieve the better `O(n^(2/3))` complexity, `tbl` values
 * smaller than `O(n^(2/3))` have to be precomputed with sieve in advance.
 *
 * @deprecated - `mertens_odd` implementation that uses `pi_tbl` is faster in practice
 *               because this implementation requires too much memory
 *
 * @param n - argument at which to evaluate `M`
 * @param tbl - table to store the calculated values
 */
template<typename T, typename I, typename TBL>
T mertens_odd(I n, TBL& tbl, T id = T(1)) {
    auto t = [&](I k) { return id; };
    auto s = [&](I k) { return castOf<T>(id, (k + 1) / 2); };
    return sum_m<T>(t, s, n, tbl, id);
}
/**
 * MertensOdd function in `O(n^(2/3))`
 *
 * @param n - argument at which to evaluate `M`
 * @param pi_tbl - table of prime_pi(n/i) for all i
 * @param pa - table of all `psz` prime numbers up to `sqrt(n)`, or `sqrt(n log n)` (it's faster)`
 */
template<typename T, typename PI_TBL>
altruct::container::sqrt_map<int64_t, T> mertens_odd(int64_t n, PI_TBL& pi_tbl, const int* pa, int psz, T id = T(1)) {
    T zero = zeroOf(id);
    auto mu = [&](T f_pe1, int p, int e) { return castOf(id, (p == 2 || e > 1) ? 0 : -1); };
    auto s1 = [&](int64_t n) { return -castOf(id, pi_tbl[n]) + ((n >= 2) ? id : zero); };
    return sum_multiplicative<T>(s1, mu, n, pa, psz, id);
}

/**
 * MertensEven function: `Sum[MoebiusMu(k), {k, 2, n, 2}]` in `O(n^(3/4))` or `O(n^(2/3))`.
 *
 * Note, to achieve the better `O(n^(2/3))` complexity, `tbl` values
 * smaller than `O(n^(2/3))` have to be precomputed with sieve in advance.
 *
 * @deprecated - `mertens_even` implementation that uses `pi_tbl` is faster in practice
 *               because this implementation requires too much memory
 *
 * @param n - argument at which to evaluate `M`
 * @param tbl - table to store the calculated values
 */
template<typename T, typename I, typename TBL>
T mertens_even(I n, TBL& tbl, T id = T(1)) {
    auto zero = zeroOf(id);
    auto t = [&](int k) { return (k > 1) ? -id : zero; };
    auto s = [&](I k) { return castOf<T>(id, (k + 1) / 2); };
    return sum_m<T>(t, s, n, tbl, id);
}
/**
 * MertensEven function in `O(n^(2/3))`
 *
 * @param n - argument at which to evaluate `M`
 * @param pi_tbl - table of prime_pi(n/i) for all i
 * @param pa - table of all `psz` prime numbers up to `sqrt(n)`, or `sqrt(n log n)` (it's faster)`
 */
template<typename T, typename PI_TBL>
altruct::container::sqrt_map<int64_t, T> mertens_even(int64_t n, PI_TBL& pi_tbl, const int* pa, int psz, T id = T(1)) {
    auto M1 = mertens_odd(n, pi_tbl, pa, psz, id);
    auto M0 = mertens(n, pi_tbl, pa, psz, id);
    int64_t q = isqrt(n), n_q = n / (q + 1);
    for (int64_t i = 1; i <= q; i++) M0[i] -= M1[i];
    for (int64_t i = 1; i <= n_q; i++) M0[n / i] -= M1[n / i];
    return M0;
}

} // math
} // altruct
