#pragma once

#include "altruct/algorithm/math/base.h"
#include "altruct/algorithm/math/primes.h"
#include "altruct/algorithm/math/factorization.h"
#include "altruct/structure/math/modulo.h"

#include <vector>
#include <algorithm>

namespace altruct {
namespace math {

/**
 * Based on the following paper:
 *   Binomial coefficients modulo prime powers by Andrew Granville
 */


/**
 * Generates a table of factorials modulo prime power `p^e` up to `n`
 *   where multiples of `p` are skipped (i.e. not taken into the product).
 *
 * Complexity: O(sz * M(log p^e))
 *   where M(n) is complexity of n-bit multiplication
 *
 * Note, multiples of `p` are skipped as required by `factorial_mod_pp`.
 */
template<typename R>
std::vector<R> factorials_mod_pp_skipped(int n, int p, int e = 1) {
    moduloX<R> r(1, powT(castOf<R>(p), e));
    std::vector<R> tbl(n + 1);
    for (int i = 0; i <= n; i++) {
        if (i % p != 0) r *= castOf<R>(i);
        tbl[i] = r.v;
    }
    return tbl;
}

/**
 * Factorial of `n` modulo prime power `p^e`
 *   where multiples of `p` are skipped (i.e. not taken into the product).
 *
 * `f = n! / (u! * p^u)`, where `u = floor(n/p)`.
 *
 * Complexity: O(n * M(log n + log p^e)),
 *   where M(n) is complexity of n-bit multiplication
 *
 * @param n - number to take the factorial of
 * @return - the pair {f, u}
 */
template<typename R, typename I>
std::pair<R, I> factorial_mod_pp_skipped_slow(I n, int p, int e) {
    moduloX<R> r(1, powT(castOf<R>(p), e));
    for (I i = 1; i <= n; i++) {
        if (i % p != 0) r *= castOf<R>(i);
    }
    return { r.v, n / p };
}

/**
 * Factorial of `n` modulo prime power `p^e`
 *   where multiples of `p` are reduced (i.e. factored out).
 *
 * `f = (n! / p^a) % p^e`, where `p^a` is the largest power of `p` dividing `n!`.
 * Note, before modulo operation gets applied, all factors `p` are removed.
 *
 * Complexity: O(n * M(log n + log p^e)),
 *   where M(n) is complexity of n-bit multiplication
 *
 * @param I - type of the number `n`
 * @param n - number to take the factorial of
 * @return - the pair {f, a}
 */
template<typename R, typename I>
std::pair<R, I> factorial_mod_pp_reduced_slow(I n, int p, int e) {
    I a = 0;
    moduloX<R> r(1, powT(castOf<R>(p), e));
    for (I i = 1; i <= n; i++) {
        r *= castOf<R>(factor_out(i, p, a));
    }
    return { r.v, a };
}

/**
 * Binomial of (n, k) modulo prime power `p^e`
 *   where multiples of `p` are reduced (i.e. factored out).
 *
 * `b = (binomial(n, k) / p^a) % p^e`, where `p^a` is the largest power of `p`.
 * Note, before modulo operation gets applied, all factors `p` are removed.
 *
 * Complexity: O(k * M(log n + log p^e)),
 *   where M(n) is complexity of n-bit multiplication
 *
 * @param I - type of the numbers `n` and `k`
 * @param n, k - numbers to take the binomial of
 * @return - the pair {b, a}
 */
template<typename R, typename I>
std::pair<R, I> binomial_mod_pp_reduced_slow(I n, I k, int p, int e) {
    k = std::min<I>(k, n - k);
    moduloX<R> r(1, powT(castOf<R>(p), e)), s = r;
    I ar = 0, as = 0;
    for (I i = 1; i <= k; i++, n--) {
        r *= castOf<R>(factor_out(n, p, ar));
        s *= castOf<R>(factor_out(i, p, as));
    }
    r /= s;
    ar -= as;
    return{ r.v, ar };
}

/**
 * Factorial of `n` modulo prime power `p^e`
 *   where multiples of `p` are reduced (i.e. factored out).
 *
 * `f = (n! / p^a) % p^e`, where `p^a` is the largest power of `p` dividing `n!`.
 * Note, before modulo operation gets applied, all factors `p` are removed.
 *
 * This is a faster variant of `factorial_mod_pp_reduced_slow` and uses a
 * look-up table up to `p^e`. See also `factorial_mod_pp_reduced` which is
 * even faster and uses a smaller look-up table up to `p*e` only.
 *
 * Complexity: O(log n / log p * M(log n + log p^e)),
 *   where M(n) is complexity of n-bit multiplication
 *
 * @param I - type of the number `n`
 * @param n - number to take the factorial of
 * @param fact_table - look-up table of `n! % p^e` up to 'p^e',
 *   see `factorials_mod_pp_skipped`
 * @return - the pair {f, a}
 */
template<typename R, typename I>
std::pair<moduloX<R>, I> factorial_mod_pp_reduced_2_modx(I n, int p, int e, const R* fact_table) {
    int pe = powT(p, e);
    // fact_table[p^e] == (p == 2 && k != 2) +1 : -1;
    bool sign = !(p == 2 && e != 2);
    moduloX<R> f(fact_table[0], castOf<R>(pe));
    I a = 0;
    while (n > 1) {
        I q = n / pe, r = n % pe;
        if (sign && (q % 2 != 0)) f = -f;
        f *= fact_table[castOf<int>(r)]; n /= p;
        a += n;
    }
    return{ f, a };
}
template<typename R, typename I>
std::pair<R, I> factorial_mod_pp_reduced_2(I n, int p, int e, const R* fact_table) {
    auto r = factorial_mod_pp_reduced_2_modx(n, p, e, fact_table);
    return{ r.first.v, r.second };
}

/**
 * Binomial of (n, k) modulo prime power `p^e`
 *   where multiples of `p` are reduced (i.e. factored out).
 *
 * Calculated as `n! / k! / (n - k)!`. See `factorial_mod_pp_reduced_2`.
 * `b = (binomial(n, k) / p^a) % p^e`, where `p^a` is the largest power of `p`.
 * Note, before modulo operation gets applied, all factors `p` are removed.
 *
 * This is a faster variant of `binomial_mod_pp_reduced_slow` and uses a
 * look-up table up to `p^e`. See also `binomial_mod_pp_reduced` which is
 * even faster and uses a smaller look-up table up to `p*e` only.
 *
 * Complexity: O(log n / log p * M(log n + log p^e)),
 *   where M(n) is complexity of n-bit multiplication
 *
 * @param I - type of the numbers `n` and `k`
 * @param n, k - numbers to take the binomial of
 * @param fact_table - look-up table of `n! % p^e` up to 'p^e',
 *   see `factorials_mod_pp_skipped`
 * @return - the pair {b, a}
 */
template<typename R, typename I>
std::pair<R, I> binomial_mod_pp_reduced_2(I n, I k, int p, int e, const R* fact_table) {
    auto fn = factorial_mod_pp_reduced_2_modx(n, p, e, fact_table);
    auto fk = factorial_mod_pp_reduced_2_modx(k, p, e, fact_table);
    auto fl = factorial_mod_pp_reduced_2_modx(I(n - k), p, e, fact_table);
    auto b = fn.first / (fk.first * fl.first);
    auto a = fn.second - (fk.second + fl.second);
    return{ b.v, a };
}

/**
 * Factorial of `n` modulo prime power `p^e`
 *   where multiples of `p` are skipped (i.e. not taken into the product).
 *
 * `f = n! / (u! * p^u)`, where `u = floor(n/p)`.
 *
 * Complexity: O(e^3 * M(log n) + e^2 * log p * M(log p^e)), ?
 *   where M(n) is complexity of n-bit multiplication
 *
 * @param I - type of the number `n`
 * @param n - number to take the factorial of
 * @param fact_table - look-up table of `n! % p^e` up to 'p*e',
 *   see `factorials_mod_pp_skipped`
 * @return - the pair {f, u}
 */
template<typename R, typename I>
std::pair<R, I> factorial_mod_pp_skipped(I n, int p, int e, const R* fact_table) {
    R m = powT<R>(p, e);
    R me = m / p * (p - 1);
    I u = n / p; int v = castOf<int, I>(n % p);
    int r = (e + 1) / 2;
    if (n <= p * e) return{ fact_table[castOf<int, I>(n)] % m, u };

    // gcd look-up table
    std::vector<std::vector<int>> gcd_table(e + 1, std::vector<int>(e + 1));
    for (int i = 0; i <= e; i++) {
        for (int j = 0; j <= e; j++) {
            gcd_table[i][j] = gcd(i, j);
        }
    }
    auto gcd_f = [&](I n, int d) { return gcd_table[castOf<int, I>(n % d)][d]; };

    std::vector<I> numerators;
    std::vector<int> denominators;
    // alpha coefficients
    std::vector<R> alphas(e);
    for (int j = 1; j < e; j++) {
        moduloX<R> alpha(1, me);
        numerators.clear();
        denominators.clear();
        for (int i = 0; i < e; i++) {
            if (i == j) continue;
            if (j < i) alpha = -alpha;
            numerators.push_back(u - i);
            denominators.push_back(absT(j - i));
        }
        fraction_reduce(numerators, denominators, gcd_f);
        // denominators should all be 1 now
        for (const auto& num : numerators) {
            alpha *= castOf<R, I>(num % me);
        }
        alphas[j] = alpha.v;
    }
    // beta coefficients
    std::vector<R> betas(r + 1);
    for (int j = 1; j <= r; j++) {
        moduloX<R> beta(1, me);
        numerators.clear();
        denominators.clear();
        for (int i = 0; i <= r; i++) {
            if (i == j) continue;
            if (j < i) beta = -beta;
            numerators.push_back(u - i);
            denominators.push_back(absT(j - i));
            if (i == 0) continue;
            numerators.push_back(u + i);
            denominators.push_back(j + i);
        }
        fraction_reduce(numerators, denominators, gcd_f);
        // denominators should all be 1 now
        for (const auto& num : numerators) {
            beta *= castOf<R, I>(num % me);
        }
        betas[j] = beta.v;
    }

    moduloX<R> f(1, m);
    // fact(u*p, p)
    for (int j = 1; j <= r; j++) {
        f *= modulo_power<R, R>(fact_table[j * p], betas[j], m);
    }
    if (p == 2 && f.v % 4 != (u | 1) % 4) f = -f;
    // fact(v, p)
    f *= moduloX<R>(fact_table[v], m);
    // bin(u*p+v, v, p)
    for (int j = 1; j < e; j++) {
        moduloX<R> bin_num(fact_table[j * p + v], m);
        moduloX<R> bin_d1(fact_table[j * p], m);
        moduloX<R> bin_d2(fact_table[v], m);
        moduloX<R> bin = bin_num / (bin_d1 * bin_d2);
        f *= modulo_power<R, R>(bin.v, alphas[j], m);
    }
    // fact(u*p+v, p) = fact(u*p, p) * fact(v, p) * bin(u*p+v, v, p)
    return{ f.v, u };
}

/**
 * Factorial of `n` modulo prime power `p^e`
 *   where multiples of `p` are reduced (i.e. factored out).
 *
 * `f = (n! / p^a) % p^e`, where `p^a` is the largest power of `p` dividing `n!`.
 * Note, before modulo operation gets applied, all factors `p` are removed.
 *
 * Complexity: O(e^3 * log^2 n / log p + e^2 * log n * M(log p^e)), ?
 *   where M(n) is complexity of n-bit multiplication
 *
 * @param I - type of the number `n`
 * @param n - number to take the factorial of
 * @param fact_table - look-up table of `n! % p^e` up to 'p*e',
 *   see `factorials_mod_pp_skipped`
 * @return - the pair {f, a}
 */
template<typename R, typename I>
std::pair<R, I> factorial_mod_pp_reduced(I n, int p, int e, const R* fact_table) {
    moduloX<R> r(1, powT(castOf<R>(p), e));
    I a = 0;
    for (int i = 0; n > 1; i++) {
        r *= factorial_mod_pp_skipped(n, p, e, fact_table).first, n /= p, a += n;
    }
    return{ r.v, a };
}

/**
 * Binomial of (n, k) modulo prime power `p^e`
 *   where multiples of `p` are reduced (i.e. factored out).
 *
 * Calculated as `n! / k! / (n - k)!`. See `factorial_mod_pp_reduced`.
 * `b = (binomial(n, k) / p^a) % p^e`, where `p^a` is the largest power of `p`.
 * Note, before modulo operation gets applied, all factors `p` are removed.
 *
 * Complexity: O(e^5 * log n + e^3 * log n * M(log p^e)), ?
 *   where M(n) is complexity of n-bit multiplication
 *
 * @param I - type of the numbers `n` and `k`
 * @param n, k - numbers to take the binomial of
 * @param fact_table - look-up table of `n! % p^e` up to 'p*e',
 *   see `factorials_mod_pp_skipped`
 * @return - the pair {b, a}
 */
template<typename R, typename I>
std::pair<R, I> binomial_mod_pp_reduced(I n, I k, int p, int e, const R* fact_table) {
    I l = n - k;
    R m = powT<R>(castOf<R>(p), e);
    moduloX<R> r(1, m), s = r;
    I a = 0;
    bool sign = !(p == 2 && e != 2);
    for (int i = 0; n > 1; i++) {
        // Since we are multiplying factorial_mod_pp_skipped(n % m)
        // instead of simply doing factorial_mod_pp_reduced(n),
        // we need to adjust the sign. This is faster than the later.
        if (sign && i >= e - 1 && (n % p) < (k % p) + (l % p)) r = -r;
        if (n > 0) r *= factorial_mod_pp_skipped(castOf<R, I>(n % m), p, e, fact_table).first, n /= p, a += n;
        if (k > 0) s *= factorial_mod_pp_skipped(castOf<R, I>(k % m), p, e, fact_table).first, k /= p, a -= k;
        if (l > 0) s *= factorial_mod_pp_skipped(castOf<R, I>(l % m), p, e, fact_table).first, l /= p, a -= l;
    }
    r /= s;
    return{ r.v, a };
}

} // math
} // altruct
