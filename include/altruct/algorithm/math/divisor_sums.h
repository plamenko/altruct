#pragma once

#include <type_traits>
#include <vector>

#include "altruct/algorithm/math/base.h"
#include "altruct/algorithm/math/polynoms.h"
#include "altruct/structure/container/sqrt_map.h"
#include "altruct/structure/math/fenwick_tree.h"

namespace altruct {
namespace math {

/**
 * Useful identities involving Dirichlet convolution.
 *
 * Denote `*` as dirichlet convolution and `.` as pointwise multiplication.
 *
 * Definitions:
 *   `e`    : dirichlet multiplicative identity: `e(n) = [n == 1]`
 *   `1`    : constant function: `1(n) = 1`
 *   `Id`   : identity function: `Id(n) = n`
 *   `Id_k` : identity power function: `Id(n) = n^k`
 *   `d`    : divisor function (k=0)
 *   `s`    : divisor sigma function (k=1)
 *   `s_k`  : divisor sigma function (kt-h powers of divisors)
 *   `phi`  : euler totient function
 *   `mu`   : moebius mu function
 *   `L`    : liouville lambda function
 *   `f`,`g`: arbitrary functions
 *   `p`    : arbitrary completely multiplicative function
 *
 * Identities:
 *   f . 1 = 1 . f = f
 *   f . e = e . f = e  // f(1) = 1
 *   f * e = e * f = f
 *   f * f^-1 = f^-1 * f = e
 *   p . (f * g) = (p . f) * (p . g)
 *   mu * 1 = e
 *   mu = 1^-1
 *   L * |mu| = e
 *   L = |mu|^-1
 *   Id_k * (Id_k . mu) = e
 *   Id_k^-1 = (Id_k . mu)
 *   d * mu = 1
 *   d = 1 * 1
 *   s_k * mu = Id_k
 *   s_k = Id_k * 1
 *   phi * 1 = Id
 *   phi = Id * mu
 *   s = phi * d
 *
 * A short diagram:
 *        1       1
 *   phi ---> Id ---> sigma
 *   phi <--- Id <--- sigma
 *        mu      mu
 */

/**
 * Dirichlet convolution of `f` and `g` up to `n` in `O(n log n)`.
 *
 * Calculates `h` where `h[n] = Sum[f(n/d) * g(d), {d|n}]`
 *
 * Where:
 *   `f` and `g` are arbitrary arithmetic functions
 *
 * @param h - table to store the result; accessed via [] operator
 * @param f, g - functions as defined above; accessed via () operator
 * @param n - bound up to which to calculate `h`; exclusive
 */
template<typename F1, typename F2, typename TBL>
void dirichlet_convolution(TBL& h, F1 f, F2 g, int n) {
    auto e0 = zeroOf(castOf(h[0], f(1))); // TODO: pass via argument?
    for (int i = 0; i < n; i++) {
        h[i] = e0;
    }
    for (int d = 1; d < n; d++) {
        for (int e = 1, i = d; i < n; i += d, e++) {
            h[i] += castOf(e0, f(d)) * castOf(e0, g(e));
        }
    }
}

/**
 * Dirichlet division of `f` with `g` up to `n` in `O(n log n)`.
 *
 * Calculates `h` such that: `h = f * g^-1`
 *
 * Where:
 *   `f` is an arbitrary arithmetic function
 *   `g` is an arbitrary arithmetic function such that:
 *       g(1) != 0, and g(1) is invertible
 *
 * @param h - table to store the result; accessed via [] operator
 * @param f, g - functions as defined above; accessed via () operator
 * @param n - bound up to which to calculate `h`; exclusive
 */
template<typename F1, typename F2, typename TBL>
void dirichlet_division(TBL& h, F1 f, F2 g, int n) {
    auto e1 = identityOf(castOf(h[0], f(1))); // TODO: pass via argument?
    auto g1 = castOf(e1, g(1));
    auto ig1 = e1 / g1;
    for (int i = 1; i < n; i++) {
        h[i] = castOf(e1, f(i));
    }
    for (int d = 1; d < n; d++) {
        h[d] *= ig1;
        for (int j = 2, i = d * 2; i < n; i += d, j++) {
            h[i] -= castOf(e1, g(j)) * h[d];
        }
    }
}

/**
 * Dirichlet inverse of `f` up to `n` in `O(n log n)`.
 *
 * Calculates `f_inv` such that: `f * f_inv = e`
 *
 * Where:
 *   `f` is an arbitrary arithmetic function such that:
 *       f(1) != 0, and f(1) is invertible
 *   `e` is the dirichlet multiplicative identity: `e(n) = [n == 1]`
 *
 * @param f_inv - table to store the result; accessed via [] operator
 * @param f - function as defined above; accessed via () operator
 * @param n - bound up to which to calculate `f_inv`; exclusive
 */
template<typename F1, typename TBL>
void dirichlet_inverse(TBL& f_inv, F1 f, int n) {
    auto e1 = identityOf(castOf(f_inv[0], f(1))), e0 = zeroOf(e1); // TODO: pass via argument?
    auto e = [&](int n){ return (n == 1) ? e1 : e0; };
    dirichlet_division(f_inv, e, f, n);
}

/**
 * Calculates all the values of a multiplicative function `f` up to `n`,
 * from the values at prime powers, in `O(n log log n)`.
 *
 * @param f - table of values of `f`; accessed via [] operator
 *            must be set to the actual value for prime powers and 1 elsewhere
 * @param n - bound up to which to calculate `f`; exclusive
 * @param pa - table of all `m` prime numbers up to `n`
 */
template<typename TBL>
void calc_multiplicative(TBL& f, int n, int* pa, int m) {
    for (int i = 0; i < m && pa[i] < n; i++) {
        for (int64_t qq = pa[i]; qq < n; qq *= pa[i]) {
            for (int q = (int)qq, l = 2, m = 2 * q; m < n; m += q, l++) {
                if (l % pa[i] != 0) f[m] *= f[q];
            }
        }
    }
}

/**
 * Dirichlet convolution of `f` and `g` up to `n` in `O(n log log n)`.
 *
 * Calculates `h` where `h[n] = Sum[f(n/d) * g(d), {d|n}]`
 *
 * Where:
 *   `f` and `g` are arbitrary arithmetic functions such that
 *   `h = f * g` is a multiplicative function
 *
 * Note that only `h` needs to be multiplicative!
 *
 * @param h - table to store the result; accessed via [] operator
 * @param f, g - functions as defined above; accessed via () operator
 * @param n - bound up to which to calculate `h`; exclusive
 * @param pa - table of all `m` prime numbers up to `n`
 */
template<typename F1, typename F2, typename TBL>
void dirichlet_convolution_multiplicative(TBL& h, F1 f, F2 g, int n, int* pa, int m) {
    auto e1 = identityOf(castOf(h[0], f(1))), e0 = zeroOf(e1); // TODO: pass via argument?
    typedef decltype(e1) T;
    for (int i = 1; i < n; i++) {
        h[i] = e1;
    }
    int q[32]; T fq[32], gq[32];
    for (int i = 0; i < m && pa[i] < n; i++) {
        int m = 0;
        for (int64_t qq = 1; qq < n; qq *= pa[i]) {
            fq[m] = castOf(e1, f((int)qq));
            gq[m] = castOf(e1, g((int)qq));
            q[m++] = (int)qq;
        }
        for (int k = 0; k < m; k++) {
            auto hq_k = e0;
            for (int j = 0; j <= k; j++) {
                hq_k += fq[k - j] * gq[j];
            }
            h[q[k]] = hq_k;
        }
    }
    calc_multiplicative(h, n, pa, m);
}

/**
 * Dirichlet division of `f` with `g` up to `n` in `O(n log log n)`.
 *
 * Calculates `h` such that: `h = f * g^-1`
 *
 * Where:
 *   `h = f * g^-1` is a multiplicative function
 *   `f` is an arbitrary arithmetic function
 *   `g` is an arbitrary arithmetic function such that:
 *       g(1) != 0, and g(1) is invertible
 *
 * @param h - table to store the result; accessed via [] operator
 * @param f, g - functions as defined above; accessed via () operator
 * @param n - bound up to which to calculate `h`; exclusive
 * @param pa - table of all `m` prime numbers up to `n`
 */
template<typename F1, typename F2, typename TBL>
void dirichlet_division_multiplicative(TBL& h, F1 f, F2 g, int n, int* pa, int m) {
    auto e1 = identityOf(castOf(h[0], f(1))), e0 = zeroOf(e1); // TODO: pass via argument?
    typedef decltype(e0) T;
    for (int i = 1; i < n; i++) {
        h[i] = e1;
    }
    int q[32]; T gq[32], hq[32];
    for (int i = 0; i < m && pa[i] < n; i++) {
        int m = 0;
        for (int64_t qq = 1; qq < n; qq *= pa[i]) {
            gq[m] = castOf(e1, g((int)qq));
            q[m++] = (int)qq;
        }
        hq[0] = e1;
        for (int k = 1; k < m; k++) {
            hq[k] = castOf(e1, f(q[k]));
            for (int j = 0; j < k; j++) {
                hq[k] -= gq[k - j] * hq[j];
            }
            h[q[k]] = hq[k];
        }
    }
    calc_multiplicative(h, n, pa, m);
}

/**
 * Dirichlet inverse of `f` up to `n` in `O(n log log n)`.
 *
 * Calculates `f_inv` such that: `f * f_inv = e`
 *
 * Where:
 *   `f` is a multiplicative function, which means:
 *       its inverse `f_inv` is also multiplicative
 *   `e` is the dirichlet multiplicative identity: `e(n) = [n == 1]`
 *
 * @param f_inv - table to store the result; accessed via [] operator
 * @param f - function as defined above; accessed via () operator
 * @param n - bound up to which to calculate `f_inv`; exclusive
 * @param pa - table of all `m` prime numbers up to `n`
 */
template<typename F1, typename TBL>
void dirichlet_inverse_multiplicative(TBL& f_inv, F1 f, int n, int* pa, int m) {
    auto e1 = identityOf(castOf(f_inv[0], f(1))), e0 = zeroOf(e1); // TODO: pass via argument?
    auto e = [&](int n){ return (n == 1) ? e1 : e0; };
    dirichlet_division_multiplicative(f_inv, e, f, n, pa, m);
}

/**
 * Calculates all the values of a completely multiplicative function `f` up to `n`,
 * from the values at primes, in `O(n)`.
 *
 * param f - table of values of `f`; accessed via [] operator
 *           must be set to the actual value for primes and 1 elsewhere
 * @param n - bound up to which to calculate `f`; exclusive
 * @param pf - table of prime factors up to `n`; pf[k] = some_prime_factor_of(k)
 */
template<typename TBL>
void calc_completely_multiplicative(TBL& f, int n, int* pf) {
    for (int i = 2; i < n; i++) {
        int p = pf[i];
        if (pf[i] != i) f[i] = f[i / p] * f[p];
    }
}

/**
 * Dirichlet convolution of `f` and `g` up to `n` in `O(n)`.
 *
 * Calculates `h` where `h[n] = Sum[f(n/d) * g(d), {d|n}]`
 *
 * Where:
 *   `f` and `g` are arbitrary arithmetic functions such that
 *   `h = f * g` is a completely multiplicative function
 *
 * Note that only `h` needs to be completely multiplicative!
 * For example, `f(n) = mu(n)` and `g(n) = sigma1(n)` are not
 * completely multiplicative, but its convolution `h(n) = n` is!
 *
 * @param h - table to store the result; accessed via [] operator
 * @param f, g - functions as defined above; accessed via () operator
 * @param n - bound up to which to calculate `h`; exclusive
 * @param pf - table of prime factors up to `n`; pf[k] = some_prime_factor_of(k)
 */
template<typename F1, typename F2, typename TBL>
void dirichlet_convolution_completely_multiplicative(TBL& h, F1 f, F2 g, int n, int *pf) {
    auto e1 = identityOf(castOf(h[0], f(1))); // TODO: pass via argument?
    auto f1 = castOf(e1, f(1)), g1 = castOf(e1, g(1));
    for (int i = 1; i < n; i++) {
        h[i] = e1;
    }
    for (int p = 2; p < n; p++) {
        if (pf[p] == p) h[p] = castOf(e1, f(p)) * g1 + castOf(e1, g(p)) * f1;
    }
    calc_completely_multiplicative(h, n, pf);
}

/**
 * Dirichlet division of `f` with `g` up to `n` in `O(n)`.
 *
 * Calculates `h` such that: `h = f * g^-1`
 *
 * Where:
 *   `h = f * g^-1` is a completely multiplicative function
 *   `f` is an arbitrary arithmetic function
 *   `g` is an arbitrary arithmetic function such that:
 *       g(1) != 0, and g(1) is invertible
 *
 * Note that only `h` needs to be completely multiplicative!
 * For example, `phi(n)` and `mu(n)` are not completely multiplicative,
 * but their division `phi(n) / mu(n) = n` is!
 *
 * @param h - table to store the result; accessed via [] operator
 * @param f, g - functions as defined above; accessed via () operator
 * @param n - bound up to which to calculate `h`; exclusive
 * @param pf - table of prime factors up to `n`; pf[k] = some_prime_factor_of(k)
 */
template<typename F1, typename F2, typename TBL>
void dirichlet_division_completely_multiplicative(TBL& h, F1 f, F2 g, int n, int* pf) {
    auto e1 = identityOf(castOf(h[0], f(1))); // TODO: pass via argument?
    for (int i = 1; i < n; i++) {
        h[i] = e1;
    }
    for (int p = 2; p < n; p++) {
        if (pf[p] == p) h[p] = castOf(e1, f(p)) - castOf(e1, g(p));
    }
    calc_completely_multiplicative(h, n, pf);
}

/**
 * Dirichlet inverse of `f` up to `n` in `O(n)`.
 *
 * Calculates `f_inv` such that: `f * f_inv = e`
 *
 * Where:
 *   `f` is an arbitrary arithmetic function such that
 *   `f_inv = f^-1` is a completely multiplicative function
 *   `e` is the dirichlet multiplicative identity: `e(n) = [n == 1]`
 *
 * Note that only `f_inv` needs to be completely multiplicative!
 * For example, `f(n) = n mu(n)` is not completely multiplicative,
 * but its inverse `f_inv(n) = n` is!
 *
 * @param f_inv - table to store the result; accessed via [] operator
 * @param f - function as defined above; accessed via () operator
 * @param n - bound up to which to calculate `f_inv`; exclusive
 * @param pf - table of prime factors up to `n`; pf[k] = some_prime_factor_of(k)
 */
template<typename F1, typename TBL>
void dirichlet_inverse_completely_multiplicative(TBL& f_inv, F1 f, int n, int* pf) {
    auto e1 = identityOf(castOf(f_inv[0], f(1))), e0 = zeroOf(e1); // TODO: pass via argument?
    auto e = [&](int n){ return (n == 1) ? e1 : e0; };
    dirichlet_division_completely_multiplicative(f_inv, e, f, n, pf);
}

/**
 * Moebius transform of `f` up to `n` in `O(n log n)`.
 *
 * g[n] = Sum[mu(n/d) * f(d), {d|n}].
 *
 * @param g - table to store the result; accessed via [] operator
 * @param f - function to transform; accessed via () operator
 * @param n - bound up to which to calculate `g`; exclusive
 */
template<typename TBL, typename F>
void moebius_transform(TBL& g, F f, int n) {
    auto e1 = identityOf(castOf(g[0], f(1))); // TODO: pass via argument?
    dirichlet_division(g, f, [&](int n){ return e1; }, n);
}

/**
 * Moebius transform of `f` up to `n` in `O(n log log n)`.
 *
 * g[n] = Sum[mu(n/d) * f(d), {d|n}].
 * `g` is a multiplicative function.
 *
 * @param g - table to store the result; accessed via [] operator
 * @param f - a function to transform; accessed via () operator
 * @param n - bound up to which to calculate `g`; exclusive
 * @param pa - table of all `m` prime numbers up to `n`
 */
template<typename TBL, typename F>
void moebius_transform_multiplicative(TBL& g, F f, int n, int* pa, int m) {
    auto e1 = identityOf(castOf(g[0], f(1))); // TODO: pass via argument?
    dirichlet_division_multiplicative(g, f, [&](int n){ return e1; }, n, pa, m);
}

/**
 * Moebius transform of `f` up to `n` in `O(n)`.
 *
 * g[n] = Sum[mu(n/d) * f(d), {d|n}].
 * `g` is a completely multiplicative function.
 *
 * @param g - table to store the result; accessed via [] operator
 * @param f - a function to transform; accessed via () operator
 * @param n - bound up to which to calculate `g`; exclusive
 * @param pf - table of prime factors up to `n`; pf[k] = some_prime_factor_of(k)
 */
template<typename TBL, typename F>
void moebius_transform_completely_multiplicative(TBL& g, F f, int n, int* pf) {
    auto e1 = identityOf(castOf(g[0], f(1))); // TODO: pass via argument?
    dirichlet_division_completely_multiplicative(g, f, [&](int n){ return e1; }, n, pf);
}

/**
 * Sieves `M` up to `n` in `O(n log log n)`.
 *
 * Where:
 *   `M` is an arbitrary function such that:
 *       its backward difference `M'` is multiplicative
 *   `t` is an arbitrary function such that:
 *       t(n) = Sum[p(k) M(n/k), {k, 1, n}]
 *   `p` is a multiplicative function, which means:
 *       its inverse `p_inv` is also multiplicative
 *
 * Then the following holds and is used for computation:
 *   t'(n) = Sum[p(d) M'(n/d), {d|n}]
 *   M'(n) = Sum[p^-1(d) t'(n/d), {d|n}]
 *   t'(n) = t(n) - t(n-1)
 *   M'(n) = M(n) - M(n-1)
 *
 * @param M - table to store the calculated values; accessed via [] operator
 * @param t, p - functions as defined above; accessed via () operator
 * @param n - bound up to which to sieve (exclusive)
 * @param pa - table of all `m` prime numbers up to `n`
 */
template<typename F1, typename F2, typename TBL>
void sieve_m_multiplicative_inv(TBL& M, F1 t, F2 p_inv, int n, int* pa, int m) {
    auto dt = [&](int n){ return (n == 1) ? t(n) : t(n) - t(n - 1); };
    dirichlet_convolution_multiplicative(M, p_inv, dt, n, pa, m);
    for (int i = 1; i < n; i++) M[i] += M[i - 1];
}
template<typename F1, typename F2, typename TBL>
void sieve_m_multiplicative(TBL& M, F1 t, F2 p, int n, int* pa, int m) {
    auto p_inv = M;
    dirichlet_inverse_multiplicative(p_inv, p, n, pa, m);
    auto _p_inv = [&](int n){ return p_inv[n]; };
    sieve_m_multiplicative_inv(M, t, _p_inv, n, pa, m);
}

/**
 * Creates a map of values of `f(n / k)` for each `k` in `[1, n]` in `O(n^(1/2))` assuming `f` is `O(1)`.
 */
template<typename I, typename T, typename F>
altruct::container::sqrt_map<I, T> make_sqrt_map(F f, I n) {
    I q = sqrtT(n), n_q = n / q;
    altruct::container::sqrt_map<I, T> tbl(n_q, n);
    for (I k = 1; k <= n_q; k++) {
        tbl[k] = f(k);
    }
    for (I m = 1; m <= q; m++) {
        tbl[n / m] = f(n / m);
    }
    tbl[0] = zeroOf(tbl[1]);
    return tbl;
}

/**
 * Sieves `M` up to `n` in `O(n log n)`.
 *
 * Where:
 *   `t` is an arbitrary function such that:
 *       t(n) = Sum[p(k) M(n/k), {k, 1, n}]
 *   `p` is an arbitrary function such that:
 *       p(1) != 0, and p(1) is invertible
 *
 * Then the following holds and is used for computation:
 *   t'(n) = Sum[p(d) M'(n/d), {d|n}]
 *   t'(n) = t(n) - t(n-1)
 *   M'(n) = M(n) - M(n-1)
 *
 * @param M - table to store the calculated values; accessed via [] operator
 * @param t, p - functions as defined above; accessed via () operator
 * @param n - bound up to which to sieve (exclusive)
 */
template<typename F1, typename F2, typename TBL>
void sieve_m(TBL& M, F1 t, F2 p, int n) {
    auto e1 = identityOf(castOf(M[0], p(1))); // TODO: pass via argument?
    auto ip1 = e1 / castOf(e1, p(1));
    M[1] = castOf(e1, t(1));
    for (int i = 2; i < n; i++) {
        M[i] = castOf(e1, t(i)) - castOf(e1, t(i - 1));
    }
    for (int d = 1; d < n; d++) {
        M[d] *= ip1;
        for (int j = 2, i = d * 2; i < n; i += d, j++) {
            M[i] -= castOf(e1, p(j)) * M[d];
        }
        if (d > 1) M[d] += M[d - 1];
    }
}

/**
 * Sieves `M` up to `n` in `O(n log n)`.
 *
 * Where:
 *   `t` is an arbitrary function such that:
 *       t(n) = Sum[M(n/k), {k, 1, n}]
 *
 * Same as `sieve_m(n, t, p, M)` with `p(n) = 1`.
 *
 * @param M - table to store the calculated values; accessed via [] operator
 * @param t - function as defined above; accessed via () operator
 * @param n - bound up to which to sieve (exclusive)
 */
template<typename F1, typename TBL>
void sieve_m(TBL& M, F1 t, int n) {
    auto e1 = identityOf(castOf(M[0], t(1))); // TODO: pass via argument?
    M[1] = castOf(e1, t(1));
    for (int i = 2; i < n; i++) {
        M[i] = castOf(e1, t(i)) - castOf(e1, t(i - 1));
    }
    for (int d = 1; d < n; d++) {
        for (int i = d * 2; i < n; i += d) {
            M[i] -= M[d];
        }
        if (d > 1) M[d] += M[d - 1];
    }
}

/**
 * Calculates `M(n)` in `O(n^(3/4))` or `O(n^(2/3))`.
 *
 * Where:
 *   `t` is an arbitrary function such that:
 *       t(n) = Sum[p(k) M(n/k), {k, 1, n}]
 *   `p` is an arbitrary function such that:
 *       p(1) != 0, and p(1) is invertible
 *   `s` is a partial sum of `p`:
 *       s(n) = Sum[p(k), {k, 1, n}]
 *
 * This allows us to compute `M` in sublinear time given
 * that we can efficiently compute `t` and `s`.
 *
 * Note, to achieve the better `O(n^(2/3))` complexity, `tbl` values
 * smaller than `O(n^(2/3))` have to be precomputed with sieve in advance.
 * More precisely, given the precomputed values up to `U`, the complexity
 * is `O(n / sqrt(U))`. I.e. if `U` is `(n/v)^(2/3)` for some `v`, the
 * complexity is `O(n^(2/3) * v^(1/3))`. Sieving is then usually `O(U v(U))`.
 * The most common scenarios are:
 *   Sieving up to `U` | Otpimal value for `U`     | Calculating `M(n)`
 *    O(U)             |  O(n^(2/3))               |  O(n^(2/3))
 *    O(U log log U)   |  O((n / log log n)^(2/3)) |  O(n^(2/3) (log log n)^(1/3))
 *    O(U log U)       |  O((n / log n)^(2/3))     |  O(n^(2/3) (log n)^(1/3))
 * Sieving can always be done in `O(U log U)` or better; See `sieve_m`.
 *
 * @param t, s - functions as defined above; accessed via () operator
 * @param n - argument at which to evaluate `M`
 * @param tbl - table to store the calculated values; accessed via [] operator
 */
template<typename T, typename I, typename F1, typename F2, typename TBL>
T sum_m(F1 t, F2 s, I n, TBL& tbl, T id) {
    T e0 = zeroOf(id);
    if (n < 1) return e0;
    if (tbl.count(n)) return tbl[n];
    auto r = castOf(e0, t(n)), p1 = castOf(e0, s(1) - s(0));
    I q = sqrtT(n);
    for (I k = 2; k <= n / q; k++) {
        r -= castOf(e0, s(k) - s(k - 1)) * sum_m<T>(t, s, n / k, tbl, id);
    }
    for (I m = 1; m < q; m++) {
        r -= castOf(e0, s(n / m) - s(n / (m + 1))) * sum_m<T>(t, s, m, tbl, id);
    }
    return tbl[n] = r / p1;
}

/**
 * Calculates `M(n)` in `O(n^(3/4))` or `O(n^(2/3))`.
 *
 * Where:
 *   `t` is an arbitrary function such that:
 *       t(n) = Sum[M(n/k), {k, 1, n}]
 *
 * Same as `sum_m(t, s, n, tbl)` with `p(n) = 1`, `s(n) = n`.
 *
 * @param t - function as defined above; accessed via () operator
 * @param n - argument at which to evaluate `M`
 * @param tbl - table to store the calculated values; accessed via [] operator
 */
template<typename T, typename I, typename F, typename TBL>
T sum_m(F t, I n, TBL& tbl, T id) {
    T e0 = zeroOf(id);
    if (n < 1) return e0;
    if (tbl.count(n)) return tbl[n];
    auto r = castOf(e0, t(n));
    I q = sqrtT(n);
    for (I k = 2; k <= n / q; k++) {
        r -= sum_m<T>(t, n / k, tbl, id);
    }
    for (I m = 1; m < q; m++) {
        r -= sum_m<T>(t, m, tbl, id) * castOf(e0, (n / m) - (n / (m + 1)));
    }
    return tbl[n] = r;
}

/**
 * Sieves Mertens up to `n` in `O(n log log n)`.
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
 * Mertens function: `Sum[moebius_mu(k), {k, 1, n}]` in `O(n^(3/4))` or `O(n^(2/3))`.
 *
 * Note, to achieve the better `O(n^(2/3))` complexity, `tbl` values
 * smaller than `O(n^(2/3))` have to be precomputed with sieve in advance.
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
 * MertensOdd function: `Sum[moebius_mu(k), {k, 1, n, 2}]` in `O(n^(3/4))` or `O(n^(2/3))`.
 *
 * Note, to achieve the better `O(n^(2/3))` complexity, `tbl` values
 * smaller than `O(n^(2/3))` have to be precomputed with sieve in advance.
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
 * MertensEven function: `Sum[moebius_mu(k), {k, 2, n, 2}]` in `O(n^(3/4))` or `O(n^(2/3))`.
 *
 * Note, to achieve the better `O(n^(2/3))` complexity, `tbl` values
 * smaller than `O(n^(2/3))` have to be precomputed with sieve in advance.
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
 * Sieves the count of square-free numbers up to n.
 *
 * Complexity: O(n log log n)
 *
 * @param sqfc - table to store the calculated values; accessed via [] operator
 * @param n - count is calculated up to `n` (exclusive)
 * @param pa - table of all `m` prime numbers up to `sqrt(n)`
 */
template<typename T, typename TBL>
void sieve_sqfree_count(TBL& sqfc, int n, const int* pa, int m, T id = T(1)) {
    auto e0 = zeroOf(id);
    sqfc[0] = 0;
    for (int i = 1; i < n; i++) {
        sqfc[i] = id;
    }
    for (int i = 0; i < m; i++) {
        int64_t p2 = isq(pa[i]);
        for (int64_t j = p2; j < n; j += p2) {
            sqfc[j] = e0;
        }
    }
    for (int i = 1; i < n; i++) {
        sqfc[i] += sqfc[i - 1];
    }
}

/**
 * Computes the count of square-free numbers up to n.
 *
 * Complexity:
 *   O(n^1/2) - if only sqfree_count(n) is needed, no preprocessing required
 *   O(n^3/5) - if sqfree_count(n/i) for all i are needed, assuming tbl preprocessed up to O(n^3/5)
 *   O(n^2/3) - if sqfree_count(n/i) for all i are needed, if no preprocessing is done
 *
 * @param n - argument at which to evaluate quare_free_count
 * @param tbl - table to store the calculated values; accessed via [] operator
 */
template<typename T, typename I, typename TBL>
T sqfree_count(I n, TBL& tbl, T id) {
    T e0 = zeroOf(id);
    if (n < 1) return e0;
    if (tbl.count(n)) return tbl[n];
    auto r = castOf(e0, n);
    I q = cbrtT(n);
    for (I m = 1; m < q; m++) {
        r -= castOf(e0, isqrt(n / m) - isqrt(n / (m + 1))) * sqfree_count(m, tbl, id);
    }
    for (I k = sqrtT(n / q); k > 1; k--) {
        r -= sqfree_count(n / sqT(k), tbl, id);
    }
    return tbl[n] = r;
}

/**
 * A helper function for `sum_phi_D_L`.
 *
 * Denote `*` as dirichlet convolution and `.` as pointwise multiplication.
 *
 * As defined in `sum_phi_D_L`, we have:
 *   phi_D = mu * g_D
 *   g_D = 1 * phi_D
 *
 * Since a completely multiplicative function (let's call it `p`)
 * distributes pointwise multiplication over Dirichlet convolution,
 * we also have:
 *   p . g_D = p . (1 * phi_D)
 *   p . g_D = (p . 1) * (p . phi_D)
 *   p . g_D = p * (p . phi_D)
 *
 * Let's define:
 *   t' = p . g_D
 *   M' = p . phi_D
 * Then by substitution:
 *   t' = p * M'
 *   t'(n) = Sum[p(d) M'(n/d), {d|n}]
 * Which is equivalent to:
 *   t(n) = Sum[p(k) M(n/k), {k, 1, n}]
 *   which can be calculated efficiently with `sum_m`
 *
 * The above two are equivalent because:
 *   t(n) = Sum[p(k) M(n/k), {k, 1, n}]
 *   t(n) - t(n-1) = Sum[p(k) (M(n/k) - M((n-1)/k)), {k, 1, n}]
 *   But `M(n/k) - M((n-1)/k)` is precisely 1 when `k | n` and 0 otherwise
 *   t(n) - t(n-1) = Sum[p(k) (M(n/k) - M((n-1)/k)), {k|n}]
 *   t(n) - t(n-1) = Sum[p(k) (M(n/k) - M(n/k-1)), {k|n}]
 *   t'(n) = Sum[p(k) M'(n/k), {k|n}]
 */
template<typename T, typename I>
std::vector<T> sum_g_L(const polynom<T>& g, int L, const std::vector<I>& vn, int U) {
    T e1 = identityOf(g[0]), e0 = zeroOf(e1);

    // initialize polynomials
    auto p = powT(polynom<T>{ e0, e1 }, L);
    auto s = polynom_sum(p);
    auto t = polynom_sum(p * g);

    //// wrapping functions that evaluate polynomials
    auto _g = [&](I n){ return g(castOf(e0, n)); };
    auto _p = [&](I n){ return p(castOf(e0, n)); };
    auto _s = [&](I n){ return s(castOf(e0, n)); };
    auto _t = [&](I n){ return t(castOf(e0, n)); };

    // preprocess `phi_D = mu * g_D` up to `U`
    I n = *std::max_element(vn.begin(), vn.end());
    if (U <= 0) U = (int)isq(icbrt(n)); // TODO: cast
    altruct::container::sqrt_map<I, T> mm(U, n);
    moebius_transform(mm, _g, U);
    // preprocess `Sum[p(k) * phi_D[k], {k, 1, n}]` up to `U`
    mm[0] = e0;
    for (int k = 1; k < U; k++) {
        mm[k] = mm[k - 1] + _p(k) * mm[k];
    }

    // calculate the values of interest with `sum_m`
    std::vector<T> v;
    for (auto k : vn) {
        mm.reset_max(k);
        v.push_back(sum_m<T>(_t, _s, k, mm, e1));
    }
    return v;
}

/**
 * Calculates `Sum[k^L euler_phi_D(k), { k, 1, n }]` in `O(n^(2/3))`.
 *
 * Where:
 *   euler_phi_D(n) = Sum[[GCD(a1, a2, ..., aD, n) == 1], {a1 <= a2 <= ... <= aD <= n}]
 *                  = Sum[mu(n/d) * binomial(D + d - 1, D), {d|n}]
 *                  = Sum[mu(n/d) * g_phi_D(d), {d|n}]
 *   euler_phi_0(n) = [n == 1]
 *   euler_phi_1(n) = euler_totient(n)
 *   euler_phi_2(n) = euler_totient_2d(n)
 *   euler_phi_3(n) = euler_totient_3d(n)
 *   ...
 *
 * @param D - totient dimension parameter; should be a small constant
 * @param L - exponent in `k^L`; should be a small constant
 * @param vn - arguments at which to evaluate the sum
 * @param U - sieving bound, if 0 is given, `n^2/3` is used for max `n` in `vn`
 * @param id - multiplicative identity in T
 */
template<typename T>
std::vector<T> sum_phi_D_L(int D, int L, const std::vector<int64_t>& vn, int U, T id = T(1)) {
    polynom<T> g_phi_D{ id };
    for (int i = 0; i < D; i++) {
        g_phi_D *= polynom<T>{ castOf(id, i), id } / castOf(id, i + 1);
    }
    return sum_g_L(g_phi_D, L, vn, U);
}
template<typename T>
T sum_phi_D_L(int D, int L, int64_t n, int U, T id = T(1)) {
    return sum_phi_D_L(D, L, std::vector<int64_t>{ n }, U, id).back();
}
/**
 * Calculates `Sum[euler_phi(k), { k, 1, n }]` in `O(n^(2/3))`.
 * More efficient than `sum_phi_D_L` for D=1, L=0.
 *
 * @param n - argument up to which to perform computation
 * @param id - multiplicative identity in T
 * @param phi - table of euler_phi up to `n^(2/3)`
 */
template<typename T, typename I>
altruct::container::sqrt_map<I, T> sum_phi(I n, T id = T(1), int* phi = nullptr) {
    auto idn = [&](I n) { return castOf(id, n); };
    auto tri = [&](I n) { auto r = castOf(id, n); return r * (r + 1) / 2; };
    int U = (int)isq(icbrt(n)); // TODO: cast
    altruct::container::sqrt_map<I, T> mm(U, n); mm[0] = zeroOf(id);
    if (phi) {
        for (int k = 1; k < U; k++) mm[k] = mm[k - 1] + phi[k];
    } else {
        moebius_transform(mm, idn, U);
        for (int k = 1; k < U; k++) mm[k] = mm[k - 1] + mm[k];
    }
    sum_m<T>(tri, idn, n, mm, id);
    return mm;
}

/**
 * Divisor Sigma k (Sum of k-th powers of divisors) up to `n`
 *
 * A faster implementation than the `O(n log n)` ones in `primes.h`.
 *
 * Complexity: O(n log log n)
 *
 * @param ds - table to store the result; accessed via [] operator
 * @param k - dvisors are taken to the k-th power
 * @param n - sum is calculated up to `n` (exclusive)
 * @param pa - table of all `m` prime numbers up to `n`
 */
template<typename T, typename TBL>
void divisor_sigma(TBL& ds, int k, int n, int* pa, int m, T id = T(1)) {
    // s_k = Id_k * 1
    auto _c1 = [&](int n){ return id; };
    if (k == 0) {
        auto _id_0 = [&](int n){ return id; };
        dirichlet_convolution_multiplicative(ds, _id_0, _c1, n, pa, m);
    } else if (k == 1) {
        auto _id_1 = [&](int n){ return castOf(id, n); };
        dirichlet_convolution_multiplicative(ds, _id_1, _c1, n, pa, m);
    } else if (k == 2) {
        auto _id_2 = [&](int n){ return sqT(castOf(id, n)); };
        dirichlet_convolution_multiplicative(ds, _id_2, _c1, n, pa, m);
    } else {
        auto _id_k = [&](int n){ return powT(castOf(id, n), k); };
        dirichlet_convolution_multiplicative(ds, _id_k, _c1, n, pa, m);
    }
}

/**
 * Calculates `S(n)` in `O(n^(3/4) / log(n))`.
 *
 * Where:
 *   `S` is a partial sum of `f`:
 *       S(n) = Sum[f(k), {k, 1, n}]
 *   `f` is a multiplicative function such that:
 *       f(p^e) = f(f_pe1, p, e)
 *       `f_pe1 = f(p^(e-1))` is provided for convenience
 *   `s1` is a partial sum of `f` over primes:
 *       s1(n) = Sum[f(p), {prime p <= n}]
 *   `bpf` is the biggest_prime_factor function
 *
 * Each `k` is represented as `k = t * bpf(k)` where `t = k / bpf(k)`.
 * The function performs DFS traversal over all values of `t` (classes of k).
 * For each `t`, aggregated contribution of `bpf(k)` over `primes >= bpf(t)` is added.
 *
 * Recursion parameters:
 *   bpf_t_val = bpf(t)
 *   bpf_t_exp = exponent of bpf(t) in t
 *   tb = bpf_t_val ^ bpf_t_exp
 *   f_tb = f(tb)
 *
 * @param s1 - map of values of `s1` as defined above at points `floor(n / k)` for `1 <= k <= n`.
 * @param f - function as defined above; accessed via () operator at values e >= 1
 * @param n - computes the sum up to `n` (inclusive)
 * @param pa - table of all `m` prime numbers up to `sqrt(n)`
 */
template<typename T, typename F1>
T sum_multiplicative_34(const altruct::container::sqrt_map<int64_t, T>& s1, const F1& f, int64_t n, const int* pa, int m, T f_tb = T(1), int bpf_t_val = 1, int bpf_t_exp = 0) {
    T id = identityOf(f_tb);

    // add contribution of `k = tb * bpf(t)`
    T ret = bpf_t_exp ? f(f_tb, bpf_t_val, bpf_t_exp + 1) : id;

    // add contribution of `k = tb * p` for all `p > bpf(t)`, meaning `p` doesn't appear in `t`
    ret += bpf_t_exp ? f_tb * (s1[n] - s1[bpf_t_val]) : s1[n];

    // recurse over values of `t`
    for (int i = 0; i < m; ++i) {
        int p = pa[i];
        int e = 0;
        T f_pe = id;
        int64_t n_next = n / p;
        int bpf_t_val_next = bpf_t_exp ? bpf_t_val : p;
        if (n_next < bpf_t_val_next) break;
        while (n_next >= bpf_t_val_next) {
            e += 1;
            f_pe = f(f_pe, p, e);
            T f_tb_next = bpf_t_exp ? f_tb : f_pe;
            int bpf_t_exp_next = bpf_t_exp ? bpf_t_exp : e;
            T f_rec = sum_multiplicative_34(s1, f, n_next, pa, i, f_tb_next, bpf_t_val_next, bpf_t_exp_next);
            ret += bpf_t_exp ? f_pe * f_rec : f_rec;
            n_next /= p;
        }
    }
    return ret;
}

template<typename T, typename S1, typename F1>
T sum_multiplicative_34(const S1& s1, const F1& f, int64_t n, const int* pa, int m, T id = T(1)) {
    auto s1_tbl = make_sqrt_map<int64_t, T>(s1, n);
    return sum_multiplicative_34<T, F1>(s1_tbl, f, n, pa, m, id, 1, 0);
}

/**
 * Traverses all the numbers up to n whose smallest prime factor is at least p_k.
 *
 * Complexity: There are O(u * n / log(n)) such numbers when `p_k ^ u > n`.
 * See: https://en.wikipedia.org/wiki/Buchstab_function.
 *
 * @param f - multiplicative function as defined in `sum_multiplicative`
 * @param n - traverses numbers up to `n` (inclusive)
 * @param k - index of the prime p_k
 * @param pa - table of all `psz` prime numbers up to `n / p_k`
 * @param visitor - `void visitor(int64_t m, T f_m)`, where `f_m = f(m)`
 */
template<typename T, typename F, typename V>
void traverse_rough_numbers(const F& f, int64_t n, int k, const int* pa, int psz, const V& visitor, int64_t m = 1, T f_m = T(1)) {
    auto pa1 = pa - 1;
    int p = pa1[k];
    int e = 0;
    T f_pe = identityOf(f_m);
    int64_t npe = n;
    int64_t mpe = m;
    while (npe >= p) {
        e += 1;
        f_pe = f(f_pe, p, e);
        npe /= p;
        mpe *= p;
        T f_mpe = f_m * f_pe;
        visitor(mpe, f_mpe);
        for (int j = k + 1; j <= psz && pa1[j] <= npe; j++) {
            traverse_rough_numbers(f, npe, j, pa, psz, visitor, mpe, f_mpe);
        }
    }
}

/**
 * Calculates `F_k(m)` defined as:
 *   F_k(m) = Sum[f(p_k^e) * F_k+1(m / p_k^e), {e >= 0}]
 *
 * Complexity: O(log(m) / log(p_k))
 */
template<typename T, typename F, typename FK>
T calc_F_k(int p_k, int64_t m, T id, const F& f, const FK& F_k1) {
    T s = F_k1(m);
    T f_pe = id;
    int64_t mpe = m;
    for (int e = 1; (mpe /= p_k) > 0; e++) {
        f_pe = f(f_pe, p_k, e);
        s += f_pe * F_k1(mpe);
    }
    return s;
}

/**
 * Calculates `S(m)` for each `m` in `[1, n]` in `O(n^(2/3))`.
 * More efficient than `sum_multiplicative_34`.
 *
 * Where:
 *   `S` is a partial sum of `f`:
 *       S(n) = Sum[f(k), {k, 1, n}]
 *   `f` is a multiplicative function such that:
 *       f(p^e) = f(f_pe1, p, e)
 *       `f_pe1 = f(p^(e-1))` is provided for convenience
 *   `s1` is a partial sum of `f` over primes:
 *       s1(n) = Sum[f(p), {prime p <= n}]
 *   `lpf` is the lowest_prime_factor function
 *
 * @param s1 - map of values of `s1` as defined above at points `floor(n / k)` for `1 <= k <= n`.
 * @param f - function as defined above; accessed via () operator at values e >= 1
 * @param n - computes the sum up to `n` (inclusive)
 * @param pa - table of all `psz` prime numbers up to `sqrt(n)`, or `sqrt(n log n)` (it's faster)
 */
template<typename T, typename S1, typename F>
altruct::container::sqrt_map<int64_t, T> sum_multiplicative(const S1& s1, const F& f, int64_t n, const int* pa, int psz, T id = T(1)) {
    // p_k - k-th prime number (p_1 = 2, p_2 = 3, p_3 = 5, p_4 = 7, ...)
    // pi(n) - number of primes up to n (inclusive); pi(p_k) == k
    // lpf(i) - lowest prime factor of i
    // F_prime(n) = Sum[f(p), {2 <= p <= n, prime p}]
    // F_k(n) = Sum[f(i), {1 <= i <= n, lpf(i) >= p_k}]

    T zero = zeroOf(id);
    auto pa1 = pa - 1; // 1-based indexing, pa1[k] = p_k

    int q = isqrt(n);   // n^(1/2)
    int c = icbrt(n);   // n^(1/3)
    int d = c;          // ~ n^(1/3)
    int64_t nd = n / d; // ~ n^(2/3)
    int h = int(nd / pa1[psz]) + 1; // h <= n^(1/6)
    int nq = int(n / (q + 1));      // ~ n^(1/2)
    int tsz = q + 1 + nq;           // ~ n^(1/2) * 2

    altruct::container::sqrt_map<int64_t, T> F_prime(q, n);
    altruct::container::sqrt_map<int64_t, T> F_k1(q, n); // F_(k+1)
    altruct::container::sqrt_map<int64_t, T> F_k(q, n);
    altruct::math::fenwick_tree<T, std::plus<T>> Ft(tsz - d + 1, std::plus<T>(), zero);

    if (n == 1) { F_k[1] = id; return F_k; }

    // step 1: build F_prime table by evaluating s1
    for (int i = 1; i <= q; i++) {
        F_prime[i] = s1(i);
        F_prime[n / i] = s1(n / i);
    }
    int last_k;

    // step 2: calculate F_k for k = pi(n^(1/3)) + 1; O(n^(2/3) / log(n))
    {
        int k; for (k = 1; k < psz && pa1[k] <= c; k++);
        int p_k = pa1[k];
        int64_t p_k2 = isq(p_k);
        // m < p_k; only 1 is included in the sum
        F_k[0] = zero;
        for (int m = 1; m < p_k; m++) {
            F_k[m] = id;
        }
        // p_k <= m < p_k^2; also includes primes >= p_k in the sum
        T b = id - F_prime[p_k - 1];
        for (int64_t m = pa1[k]; m <= q; m++) {
            F_k[m] = b + F_prime[m];
        }
        for (int i = nq; i >= 1; i--) {
            int64_t m = n / i; if (m >= p_k2) break;
            F_k[m] = b + F_prime[m];
        }
        // p_k^2 <= m < n <= p_k^3; also includes semiprimes with factors >= p_k in the sum
        for (int i = c; i >= 1; i--) {
            int64_t m = n / i; if (m < p_k2) continue;
            T s2 = zero;
            for (int j = k; j <= psz && isq(pa1[j]) <= m; j++) {
                int p_j = pa1[j];
                T f_p = f(id, p_j, 1);
                T f_p2 = f(f_p, p_j, 2);
                s2 += f_p2 + f_p * (F_prime[m / p_j] - F_prime[p_j]);
            }
            F_k[m] = b + F_prime[m] + s2;
        }
        last_k = k;
    }

    // step 3: calculate F_k for k = {pi(n^(1/3)), ..., pi(h) + 1}; O(n^(2/3))
    {
        auto get_Ft_k1 = [&](int64_t m) {
            if (m >= nd) return F_k1[m];
            return Ft.get_sum((m <= q) ? m : tsz - n / m, zero);
        };
        auto update_Ft_k = [&](int64_t m, T f_m) {
            Ft.add((m <= q) ? m : tsz - n / m, f_m);
        };
        // store sums to the fenwick tree
        for (int64_t m = 1; m <= q; m++) Ft.add(m, F_k[m]), Ft.add(m + 1, -F_k[m]);
        for (int i = nq; i > d; i--) Ft.add(tsz - i, F_k[n / i]), Ft.add(tsz - i + 1, -F_k[n / i]);
        // main loop
        for (int k = last_k - 1; pa1[k] > h; k--) {
            int p_k = pa1[k];
            // m >= nd; compute directly
            F_k1.swap(F_k);
            for (int i = d; i >= 1; i--) {
                int64_t m = n / i;
                F_k[m] = calc_F_k(p_k, m, id, f, get_Ft_k1);
            }
            // m < nd; iterate over rough numbers and maintain sums in the fenwick tree
            traverse_rough_numbers(f, nd - 1, k, pa, psz, update_Ft_k, 1, id);
            last_k = k;
        }
        // read sums from the fenwick tree
        for (int64_t m = 0; m <= q; m++) F_k[m] = get_Ft_k1(m);
        for (int i = nq; i > d; i--) F_k[n / i] = get_Ft_k1(n / i);
    }

    // step 4: calculate F_k for k = {pi(h), ..., 1}; O(n^(2/3) / log(n))
    {
        auto get_F_k1 = [&](int64_t m) { return F_k1[m]; };
        for (int k = last_k - 1; k >= 1; k--) {
            int p_k = pa1[k];
            F_k1.swap(F_k);
            for (int64_t m = 1; m <= q; m++) {
                F_k[m] = calc_F_k(p_k, m, id, f, get_F_k1);
            }
            for (int i = nq; i >= 1; i--) {
                int64_t m = n / i;
                F_k[m] = calc_F_k(p_k, m, id, f, get_F_k1);
            }
            last_k = k;
        }
    }

    return F_k;
}

} // math
} // altruct
