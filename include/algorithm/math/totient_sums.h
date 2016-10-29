#pragma once

#include <vector>

#include "algorithm/math/base.h"
#include "algorithm/math/primes.h"
#include "algorithm/math/polynoms.h"
#include "algorithm/math/sums.h"
#include "structure/container/sqrt_map.h"

namespace altruct {
namespace math {

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
 * @param n - bound up to which to sieve
 * @param t, p - functions as defined above
 * @param M - table to store the calculated values
 */
template<typename T, typename F1, typename F2, typename MAP>
void sieve_m(int n, F1 t, F2 p, MAP& M) {
	T p1 = p(1);
	T ip1 = identityOf(p1) / p1;
	M[0] = t(0); // should be zero
	for (int i = 1; i < n; i++) {
		M[i] = t(i) - t(i - 1);
	}
	for (int d = 1; d < n; d++) {
		M[d] *= ip1;
		for (int j = 2, i = d * 2; i < n; i += d, j++) {
			M[i] -= p(j) * M[d];
		}
		M[d] += M[d - 1];
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
 */
template<typename T, typename F1, typename MAP>
void sieve_m(int n, F1 t, MAP& M) {
	M[0] = t(0); // should be zero
	for (int i = 1; i < n; i++) {
		M[i] = t(i) - t(i - 1);
	}
	for (int d = 1; d < n; d++) {
		for (int i = d * 2; i < n; i += d) {
			M[i] -= M[d];
		}
		M[d] += M[d - 1];
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
 * Note that the above relation for `t` also holds if:
 *   M(n) = Sum[p(k) * f(k), {k, 1, n}]
 *   t(n) = Sum[p(k) * g(k), {k, 1, n}]
 *   `p` is a completely-multiplicative function:
 *       p(n * m) = p(n) * p(m)
 *   `g` and `f` are Moebius transforms of each other:
 *       g(n) = Sum[f(d), {d|n}]
 *       f(n) = Sum[mu(n/d) * g(d), {d|n}]
 * This allows us to compute `M` in sublinear time given
 * that we can efficiently compute `t` and `s`.
 *
 * Note, to achieve the better `O(n^(2/3))` complexity, `tbl` values
 * smaller than `O(n^(2/3))` have to be precomputed with sieve in advance.
 * More precisely, given the precomputed values up to `U`, the complexity
 * is `O(n / sqrt(U))`. I.e. if `U` is `(n/v)^(2/3)` for some `v`, the
 * complexity is `O(n^(2/3) * v^(1/3))`. Sieving is then usually `O(n v(n))`.
 * The most common scenarios are:
 *   Sieving up to `U` | Otpimal value for `U`     | Calculating `M(n)`
 *    O(U)             |  O(n^(2/3))               |  O(n^(2/3))
 *    O(U log log U)   |  O((n / log log n)^(2/3)) |  O(n^2/3 (log log n) ^ 1/3)
 *    O(U log U)       |  O((n / log n)^(2/3))     |  O(n^2/3 (log n) ^ 1/3)
 * Sieving can always be done in `O(n log n)` or better; See `sieve_m`.
 *
 * @param n - argument at which to evaluate `M`
 * @param t, s - functions as defined above
 * @param tbl - table to store the calculated values
 */
template<typename T, typename I, typename F1, typename F2, typename MAP>
T sum_m(I n, F1 t, F2 s, MAP& tbl) {
	if (n < 1) return zeroT<T>::of(t(1));
	if (tbl.count(n)) return tbl[n];
	T r = t(n), p1 = s(1) - s(0);
	I q = sqrtT(n);
	for (I k = 2; k <= n / q; k++) {
		r -= (s(k) - s(k - 1)) * sum_m<T>(n / k, t, s, tbl);
	}
	for (I m = 1; m < q; m++) {
		r -= (s(n / m) - s(n / (m + 1))) * sum_m<T>(m, t, s, tbl);
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
 * Same as `sum_m(n, t, s, tbl)` with `p(n) = 1`, `s(n) = n`.
 */
 template<typename T, typename I, typename F, typename MAP>
T sum_m(I n, F t, MAP& tbl) {
	if (n < 1) return zeroT<T>::of(t(1));
	if (tbl.count(n)) return tbl[n];
	T r = t(n);
	I q = sqrtT(n);
	for (I k = 2; k <= n / q; k++) {
		r -= sum_m<T>(n / k, t, tbl);
	}
	for (I m = 1; m < q; m++) {
		r -= T((n / m) - (n / (m + 1))) * sum_m<T>(m, t, tbl);
	}
	return tbl[n] = r;
}

/**
 * Mertens function: `Sum[moebius_mu(k), {k, 1, n}]` in `O(n^(3/4))` or `O(n^(2/3))`.
 *
 * Note, to achieve the better `O(n^(2/3))` complexity, `tbl` values
 * smaller than `O(n^(2/3))` have to be precomputed with sieve in advance.
 *
 * @param tbl - table to store the calculated values
 */
template<typename T, typename I, typename MAP>
T mertens(I n, MAP& tbl, T id = T(1)) {
	// p = 1, f = mu, g = delta, t = 1
	return sum_m<T>(n, [&](I k){ return id; }, tbl);
}

/**
 * A helper function for `calc_sum_phi_D_L`.
 */
template<typename T, typename CAST_T>
std::vector<T> sum_g_L(const polynom<T>& g, int L, const std::vector<int64_t>& vn, int U, T id, CAST_T castT) {
	typedef polynom<T> poly;
	T e0 = zeroT<T>::of(id);

	// initialize polynomials
	auto p = powT(poly{ e0, id }, L);
	auto s = polynom_sum(p);
	auto t = polynom_sum(p * g);

	// wrapping functions that evaluate polynomials
	auto _g = [&](int64_t n){ return g(castT(n)); };
	auto _p = [&](int64_t n){ return p(castT(n)); };
	auto _s = [&](int64_t n){ return s(castT(n)); };
	auto _t = [&](int64_t n){ return t(castT(n)); };

	// preprocess `U` values of `Sum[p(k) * f[k], {k, 1, U}]`
	int64_t n = *std::max_element(vn.begin(), vn.end());
	if (U <= 0) U = (int)isq(icbrt(n));
	altruct::container::sqrt_map<int64_t, T> msf(U, n);
	moebius_transform(msf, U, _g);
	for (int k = 1; k < U; k++) {
		msf[k] = msf[k - 1] + _p(k) * msf[k];
	}

	// calculate the values of interest
	std::vector<T> v;
	for (auto k : vn) {
		msf.reset_max(k);
		v.push_back(sum_m<T>(k, _t, _s, msf));
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
 * @param vn - values of `n` to calculate
 * @param U - sieving bound, if 0 is given, `n^2/3` is used for max `n` in `vn`
 * @param id - multiplicative identity in T
 * @param castT - casts int64_t to T
 */
template<typename T, typename CAST_T>
std::vector<T> sum_phi_D_L(int D, int L, const std::vector<int64_t>& vn, int U, T id, CAST_T castT) {
	typedef polynom<T> poly;
	T e0 = zeroT<T>::of(id);
	poly g_phi_D{ id };
	for (int i = 0; i < D; i++) {
		g_phi_D *= poly{ id * i, id } / (id * (i + 1));
	}
	return sum_g_L(g_phi_D, L, vn, U, id, castT);
}
template<typename T, typename CAST_T>
T sum_phi_D_L(int D, int L, int64_t n, int U, T id, CAST_T castT) {
	return sum_phi_D_L(D, L, std::vector<int64_t>{ n }, U, id, castT).back();
}

/**
 * Calculates sum of primes: `Sum[p(k), {k, 1, n}]` in `O(n^(5/7))`.
 *
 * @param p - array of prime numbers up to `sqrt(n)` inclusive
 */
template<typename T, typename I>
T sum_primes(I n, const int* p, T id = T(1)) {
	if (n < 1) return zeroT<T>::of(id);
	// Initially, we start with the sum of all numbers:
	// d(i) =  Sum[k, {2 <= k <= i}]
	// After each round j, all multiples of a prime p(j) get eliminated:
	// d(i) = Sum[k, {2 <= k <= i, spf(k) > p(j) || is_prime(k)}]
	// spf(k) = smallest prime factor of k
	I q = sqrtT(n);
	container::sqrt_map<I, T> d(q, n);
	for (int l = 1; l <= q; l++) {
		I i = n / l;
		d[i] = id * i * (i + 1) / 2 - 1;
	}
	for (int i = n / q - 1; i >= 1; i--) {
		d[i] = id * i * (i + 1) / 2 - 1;
	}
	for (int j = 0; p[j] && p[j] <= q; j++) {
		I pj = p[j]; I p2 = sqT(pj);
		I l_max = std::min(q, n / p2);
		for (I l = 1; l <= l_max; l++) {
			I i = n / l;
			d[i] -= (d[i / pj] - d[pj - 1]) * pj;
		}
		for (I i = n / q - 1; i >= p2; i--) {
			d[i] -= (d[i / pj] - d[pj - 1]) * pj;
		}
	}
	return d[n];
}

} // math
} // altruct
