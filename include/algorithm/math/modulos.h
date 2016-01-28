#pragma once

#include "structure/math/modulo.h"
#include "structure/math/quadratic.h"

#include <algorithm>

namespace altruct {
namespace math {

/**
 * Jacobi symbol
 * 
 * For prime `m`, this is equivalent to Legendre symbol:
 *    0 - if `n` is `0` mod `m`
 *   +1 - if `n` is a quadratic residue mod `m`
 *   -1 - if `n` is a quadratic nonresidue mod `m`
 *
 * For composite m
 */
template <typename I>
int jacobi(I n, I m) {
	int e, j = 1;
	for (;;) {
		if (m == 1) return j;
		n %= m;
		if (n == 0) return 0;
		for (e = 0; (n % 2) == 0; e++) n /= 2;
		if ((e % 2) == 1 && (m % 8) == 3) j = -j;
		if ((e % 2) == 1 && (m % 8) == 5) j = -j;
		if ((n % 4) == 3 && (m % 4) == 3) j = -j;
		std::swap(n, m);
	}
}

/**
 * Square root of `n` modulo `p`
 */
template <typename T, int ID>
modulo<T, ID> sqrt_cipolla(modulo<T, ID> n) {
	typedef modulo<T, ID> mod;
	const T& p = mod::M;
	// find a quadratic nonresidue `d` modulo `p`
	mod a = 0, d = 0;
	do {
		a += 1;
		d = a * a - n;
	} while (powT(d, (p - 1) / 2) == 1); // jacobi(d, p) == 1
	// r = (a + sqrt(d)) ^ ((p + 1) / 2)
	typedef quadratic<mod, -1> quad;
	quad::D = d;
	quad r = powT(quad(a, 1), (p + 1) / 2);
	return r.a;
}

/**
 * Square root of `n` modulo `p`
 */
template <typename I>
I sqrt_cipolla(I n, I p) {
	typedef modulo<I, -1> mod;
	mod::M = p;
	return sqrt_cipolla(mod(n)).v;
}

/**
 * Square root of `n` modulo `p^k`
 */
template <typename I>
I sqrt_hensel_lift(I n, int p, int k) {
	typedef modulo<I, -1> mod;
	mod::M = p; int p_i = p, p_k = powT(p, k);
	// f(r) == r^2 - n; f'(r) == 2r;
	mod r = sqrt_cipolla(mod(n));
	for (int i = 1; i < k; i *= 2) {
		int phi = p_i / p * (p - 1); // euler_phi(p_i)
		mod u = powT(r * 2, -1 + phi); // f'(r) ^-1
		mod::M = p_i = (i * 2 < k) ? p_i * p_i : p_k; // lift modulus
		mod v = r * r - n; // f(r)
		r -= u * v;
	}
	return r.v;
}

} // math
} // altruct
