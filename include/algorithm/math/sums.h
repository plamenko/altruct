#pragma once

#include "algorithm/math/base.h"
#include "algorithm/math/recurrence.h"
#include "structure/container/sqrt_map.h"

namespace altruct {
namespace math {

/**
 * Calculates `Sum[f[k], {k, a, b}]` in `O(b - a)`.
 *
 * @param f - `f[n]`
 */
template<typename T, typename I, typename F>
T sum(F f, I a, I b, T zero = T(0)) {
	T r = zero;
	for (I k = b; k >= a; k--) {
		r += f(k);
	}
	return r;
}

/**
 * Calculates sum of powers: `Sum[k^p, {k, 1, n}]` in `O(p^2)`.
 *
 * Note: for powers bigger than 3, T must be a field.
 * Example fields that are suitable: double, fraction, modulo.
 * 
 * @param B - vector of Bernoulli numbers of the second kind up to `p`.
 */
template<typename T, typename I>
T sum_pow(int p, I n, const std::vector<T>& B) {
	T e0 = zeroT<T>::of(B[0]), e1 = identityT<T>::of(B[0]);
	if (p == 0) return e1 * n;
	if (p == 1) return e1 * n * (n + 1) / 2;
	if (p == 2) return e1 * n * (n + 1) * (2 * n + 1) / 6;
	if (p == 3) return sqT<T>(sum_pow<T, I>(1, n, B));
	//Faulhaber's formula
	T r = e0;
	T n_k = e1;
	T bin = e1;
	for (int k = 0; k <= p; k++) {
		n_k *= n;
		bin *= p - k + 1;
		bin /= k + 1;
		r += bin * B[p - k] * n_k;
	}
	return r / (p + 1);
}

/**
 * Calculates `Sum[k^p, {k, 1, n}]` in `O(p^2)`.
 *
 * Note: for powers bigger than 3, T must be a field.
 * Example fields that are suitable: double, fraction, modulo.
 */
template<typename T, typename I>
T sum_pow(int p, I n, T id = T(1)) {
	static std::vector<T> B;
	if (B.size() < p + 1) {
		int sz = (int)(B.size() + B.size() / 2);
		B = bernoulli_b<T>(max(p, sz), id);
	}
	return sum_pow(p, n, B);
}

/**
 * Calculates `Sum[f[n/k], {k, 1, n}]` in `O(sqrt(n))`.
 *
 * @param f - `f[n]`
 */
template<typename T, typename I, typename F>
T sum_sqrt(F f, I n, T zero = T(0)) {
	if (n < 1) return zero;
	I q = sqrtT(n);
	T r = zero;
	for (I k = 1; k <= n / q; k++) {
		r += f(n / k);
	}
	for (I m = 1; m < q; m++) {
		r += f(m) * ((n / m) - (n / (m + 1)));
	}
	return r;
}

/**
 * Calculates `Sum[f[k] * g[n/k], {k, 1, n}]` in `O(sqrt(n))`.
 *
 * @param sf - `sf[n] = Sum[f[k], {k, 1, n}]`
 */
template<typename T, typename I, typename F1, typename F2>
T sum_sqrt2m(F1 sf, F2 g, I n, T zero = T(0)) {
	if (n < 1) return zero;
	I q = sqrtT(n);
	T r = zero;
	T sf0 = sf(n);
	for (I k = 1; k <= n / q; k++) {
		r += (sf(k) - sf(k - 1)) * g(n / k); // f(k) * g(n / k)
	}
	for (I m = 1; m < q; m++) {
		T sf1 = sf(n / (m + 1));
		r += (sf0 - sf1) * g(m);
		sf0 = sf1;
	}
	return r;
}

/**
 * Calculates `Sum[f[k] * g[n/k], {k, 1, n}]` in `O(sqrt(n))`.
 *
 * @param sf - `sf[n] = Sum[f[k], {k, 1, n}]`
 */
template<typename T, typename I, typename F1, typename F2, typename F3>
T sum_sqrt2m(F1 f, F2 sf, F3 g, I n, T zero = T(0)) {
	if (n < 1) return zero;
	I q = sqrtT(n);
	T r = zero;
	T sf0 = sf(n);
	for (I k = 1; k <= n / q; k++) {
		r += f(k) * g(n / k);
	}
	for (I m = 1; m < q; m++) {
		T sf1 = sf(n / (m + 1));
		r += (sf0 - sf1) * g(m);
		sf0 = sf1;
	}
	return r;
}

/**
 * Calculates `Sum[f[k, n/k], {k, 1, n}]` in `O(sqrt(n))`.
 *
 * @param sf - `sf[n, m] = Sum[f[k, m], {k, 1, n}]`
 */
template<typename T, typename I, typename F>
T sum_sqrt2(F sf, I n, T zero = T(0)) {
	if (n < 1) return zero;
	I q = sqrtT(n);
	T r = zero;
	T sf0 = sf(n, 1);
	for (I k = 1; k <= n / q; k++) {
		r += sf(k, n / k) - sf(k - 1, n / k); // f(k, n / k)
	}
	for (I m = 1; m < q; m++) {
		r += sf(n / m, m) - sf(n / (m + 1), m);
	}
	return r;
}

/**
 * Calculates `sum_m(n) = Sum[p(k) * f(k), {k, 1, n}]` in `O(n^(3/4))` or `O(n^(2/3))`.
 *
 * Where:
 *   `p` is a completely-multiplicative function:
 *       p(n * m) = p(n) * p(m)
 *   `g` and `f` are Moebius transforms of each other:
 *       g(n) = Sum[f(d), {d|n}]
 *       f(n) = Sum[mu(n/d) * g(d), {d|n}]
 *   sp(n) = Sum[p(k), {k, 1, n}]
 *   st(n) = Sum[p(k) * g(k), {k, 1, n}]
 *
 * Then the following relation holds and is used for calculation:
 *   st(n) = Sum[p(k) sum_m(n/k), {k, 1, n}]
 *
 * Note, to achieve the better `O(n^(2/3))` complexity, `tbl` values
 * smaller than `O(n^(2/3))` have to be precomputed with sieve in advance.
 *
 * @param tbl - table to store the calculated values
 */
template<typename T, typename I, typename F1, typename F2, typename MAP>
T sum_m(I n, F1 st, F2 sp, MAP& tbl) {
	if (n < 1) return zeroT<T>::of(st(1));
	if (tbl.count(n)) return tbl[n];
	T r = st(n);
	I q = sqrtT(n);
	for (I k = 2; k <= n / q; k++) {
		r -= (sp(k) - sp(k - 1)) * sum_m<T>(n / k, st, sp, tbl);
	}
	for (I m = 1; m < q; m++) {
		r -= (sp(n / m) - sp(n / (m + 1))) * sum_m<T>(m, st, sp, tbl);
	}
	return tbl[n] = r;
}

/**
 * Same as `sum_m(n, st, sp, tbl)` with `p(n) = 1`, `sp(n) = n`.
 */
template<typename T, typename I, typename F, typename MAP>
T sum_m(I n, F st, MAP& tbl) {
	if (n < 1) return zeroT<T>::of(st(1));
	if (tbl.count(n)) return tbl[n];
	T r = st(n);
	I q = sqrtT(n);
	for (I k = 2; k <= n / q; k++) {
		r -= sum_m<T>(n / k, st, tbl);
	}
	for (I m = 1; m < q; m++) {
		r -= T((n / m) - (n / (m + 1))) * sum_m<T>(m, st, tbl);
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
	return sum_m<T>(n, [&](I k){ return id; }, tbl);
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
		I l_max = min(q, n / p2);
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
