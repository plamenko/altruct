#pragma once

#include "base.h"
#include "structure/math/polynom.h"
#include "structure/math/modulo.h"

#include <vector>

namespace altruct {
namespace math {

/**
 * n-th element of a linear recurrence
 *
 * f[i] = f_init[i]
 * f[n+1] = sum{f[n-i] * f_coeff[i]}
 *
 * @param f_coeff - coefficients
 * @param f_init - initial values
 */
template<typename T>
T linear_recurrence(const std::vector<T> &f_coeff, const std::vector<T> &f_init, long long n) {
	// characteristic polynomial: p(x) = 0
	int L = (int) f_coeff.size();
	polynom<T> p(vector<T>(L + 1));
	p[L] = 1;
	for (int i = 0; i < L; i++)
		p[L - 1 - i] = -f_coeff[i];
	// x^n % p(x)
	typedef modulo<polynom<T>, 1> polymod;
	polymod::M = p;
	polynom<T> x = { 0, 1 };
	polymod xn = powT<polymod>(x, n);
	// f[n]
	T r = 0;
	for (int i = 0; i < L; i++) {
		r += xn.v[i] * f_init[i];
	}
	return r;
}

/**
 * n-th Fibonacci number
 */
template<typename T>
T fibonacci(long long n) {
	return linear_recurrence<T>({ 1, 1 }, {0, 1}, n);
}

/**
 * n-th element of the Lucas L sequence
 */
template<typename T>
T lucas_l(long long n) {
	return linear_recurrence<T>({ 1, 1 }, { 2, 1 }, n);
}

/**
 * n-th element of the Lucas U(p, q) sequence
 */
template<typename T>
T lucas_u(T p, T q, long long n) {
	return linear_recurrence<T>({ p, -q }, { 0, 1 }, n);
}

/**
 * n-th element of the Lucas V(p, q) sequence
 */
template<typename T>
T lucas_v(T p, T q, long long n) {
	return linear_recurrence<T>({ p, -q }, { 2, p }, n);
}

/*
 * First B0...Bn Bernoulli numbers of the second kind
 */
template<typename T>
std::vector<T> bernoulli_b(int n) {
	std::vector<T> a(n+1), b(n+1);
	for (int m = 0; m <= n; m++) {
		a[m] = T(1) / (m + 1);
		for (int j = m; j >= 1; j--) {
			a[j - 1] = (a[j - 1] - a[j]) * j;
		}
		b[m] = a[0];
	}
	return b;
}

/**
 * Finds the characteristic polynomial of a linearly recurrent sequence
 *
 * Works for an arbitrary field. The field requirement means that all non-zero
 * elements need to have a multiplicative inverse.
 * 
 * @param a - the first 2n elements (or more) of the sequence,
 *            where n is the degree of the polynomial
 */
template<typename T>
polynom<T> berlekamp_massey(std::vector<T> a) {
	int n = (int)a.size() / 2;
	int m = 2 * n - 1;
	polynom<T> R, R0, R1, V, V0, V1, Q;
	R0[m + 1] = 1;
	for (int i = 0; i <= m; i++) {
		R1[i] = a[m - i];
	}
	V0[0] = 0;
	V1[0] = 1;
	while (R1.deg() >= n) {
		Q = R0 / R1;
		R = R0 % R1;
		V = V0 - Q * V1;
		V0 = V1; V1 = V;
		R0 = R1; R1 = R;
	}
	return V1 / V1[V1.deg()];
}

} // math
} // altruct
