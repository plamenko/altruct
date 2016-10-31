#pragma once

#include "algorithm/math/base.h"
#include "algorithm/math/recurrence.h"
#include "structure/math/polynom.h"

namespace altruct {
namespace math {

/**
 * Calculates `Sum[(a * k + b) / q, {k, 0, n - 1}]` in `O(log min(q, n))`.
 * 
 * Note: `a` and `b` must be non-negative integers, `q` must be a positive integer.
 */
template<typename I>
I sum_ratio(I a, I b, I q, I n) {
	I s = 0;
	while (n > 0) {
		I n1 = n - 1;
		s += (b / q) * n + (a / q) * n * n1 / 2;
		b %= q, a %= q; if (a == 0) break;
		n = (a * n1 + b) / q;
		b = (q - 1) - b, std::swap(a, q);
		s += n * n1;
		s = -s;
	}
	// we could keep track of the sign above, but the sum can't be negative
	return s < 0 ? -s : s;
}

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
	if (p == 2) return e1 * n * (n + 1) * (n * 2 + 1) / 6;
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
		B = bernoulli_b<T>(std::max(p, sz), id);
	}
	return sum_pow(p, n, B);
}

/**
 * Calculates `Sum[k^m x^k, {k, 1, n}]` in `O(m^2)`.
 *
 * Note: `x != 1` must hold. For `x == 1` use `sum_pow`.
 */
template<typename T, typename I>
T sum_powx(int m, T x, I n) {
	T T0 = zeroT<T>::of(x), T1 = identityT<T>::of(x);
	auto Tn = T1 * n;
	polynom<T> p{ 0, 1 }, q{ 1 }, z{ 0, -1, 1 };
	for (int k = 1; k <= m; k++) {
		auto Tk = T1 * k;
		p = z * p.derivative() - (polynom<T>{Tn, Tk - Tn}) * p;
		q = z * q.derivative() - (polynom<T>{T0, Tk - T0}) * q;
	}
	return (powT(x, n) * p(x) - q(x)) / (powT(x - T1, m + 1)) - (m == 0 ? T1 : T0);
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

} // math
} // altruct
