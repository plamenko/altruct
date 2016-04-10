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
 */
template<typename T, typename I>
T sum_pow(int p, I n, T zero = T(0)) {
	T e0 = zero, e1 = identityT<T>::of(zero);
	if (p == 0) return e1 * n;
	if (p == 1) return e1 * n * (n + 1) / 2;
	if (p == 2) return e1 * n * (n + 1) * (2 * n + 1) / 6;
	if (p == 3) return sqT<T>(sum_pow<T, I>(1, n));
	//Faulhaber's formula
	std::vector<T> B = bernoulli_b<T>(p, e1);
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
		r += ((n / m) - (n / (m + 1))) * f(m);
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
