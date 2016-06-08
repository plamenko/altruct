#pragma once

#include "algorithm/math/base.h"
#include "algorithm/math/recurrence.h"

namespace altruct {
namespace math {

/**
 * Calculates `Sum[f[k], {k, a, b}]` in `O(b - a)`.
 *
 * @param f - `f[n]`
 */
template<typename T, typename I, typename F>
T sum(const F& f, I a, I b, T zero = T(0)) {
	T r = zero;
	for (I k = a; k <= b; k++) {
		r += f(k);
	}
	return r;
}

/**
 * Calculates `Sum[k^p, {k, 1, n}]` in `O(p^2)`.
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
	return sum_pow(p, n, bernoulli_b<T>(p, id));
}

/**
 * Calculates `Sum[f[n/k], {k, 1, n}]` in `O(sqrt(n))`.
 *
 * @param f - `f[n]`
 */
template<typename T, typename I, typename F>
T sum_sqrt(const F& f, I n, T zero = T(0)) {
	I q = sqrtT(n);
	T r = zero;
	if (n < 1) return r;
	for (I k = n / q; k > 0; k--) {
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
T sum_sqrt2m(const F1& sf, const F2& g, I n, T zero = T(0)) {
	I q = sqrtT(n);
	T r = zero;
	if (n < 1) return r;
	for (I k = n / q; k > 0; k--) {
		r += (sf(k) - sf(k - 1)) * g(n / k); // f(k) * g(n / k)
	}
	T sf0 = sf(n);
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
T sum_sqrt2m(const F1& f, const F2& sf, const F3& g, I n, T zero = T(0)) {
	I q = sqrtT(n);
	T r = zero;
	if (n < 1) return r;
	for (I k = n / q; k > 0; k--) {
		r += f(k) * g(n / k);
	}
	T sf0 = sf(n);
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
T sum_sqrt2(const F& sf, I n, T zero = T(0)) {
	I q = sqrtT(n);
	T r = zero;
	if (n < 1) return r;
	for (I k = n / q; k > 0; k--) {
		r += sf(k, n / k) - sf(k - 1, n / k); // f(k, n / k)
	}
	for (I m = 1; m < q; m++) {
		r += sf(n / m, m) - sf(n / (m + 1), m);
	}
	return r;
}

} // math
} // altruct
