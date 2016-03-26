#pragma once

#include "structure/math/modulo.h"
#include "structure/math/quadratic.h"
#include "structure/math/prime_holder.h"

#include <vector>
#include <algorithm>

namespace altruct {
namespace math {

/**
 * Chinese Remainder
 *
 * Calculates `a` and `n` so that:
 *   n = lcm(n1, n2)
 *   a % n1 == a1
 *   a % n2 == a2
 *   0 <= a < n
 * `n1` and `n2` don't have to be coprime.
 * Correctly handles 64bit result for long long type.
 */
template<typename T>
void chinese_remainder(T &a, T&n, T a1, T n1, T a2, T n2) {
	T e0 = zeroT<T>::of(a);
	T ni1, ni2; T g = gcd_ex(n1, n2, &ni1, &ni2);
	if ((a2 - a1) % g != e0) { n = e0, a = e0; return; }
	T t1 = modulo_multiply(a1, ni2, n1);
	T t2 = modulo_multiply(a2, ni1, n2);
	n1 /= g; n2 /= g; n = n1 * n2 * g;
	a = modulo_multiply(t1, n2, n) + modulo_multiply(t2, n1, n);
	return modulo_normalize(a, n);
}
template<typename T>
void chinese_remainder(T &ar, T&nr, T a, T n) {
	chinese_remainder(ar, nr, ar, nr, a, n);
}

/**
 * Jacobi symbol
 * 
 * For prime `m`, this is equivalent to Legendre symbol:
 *    0 - if `n` is `0` mod `m`
 *   +1 - if `n` is a quadratic residue mod `m`
 *   -1 - if `n` is a quadratic nonresidue mod `m`
 *
 * For composite `m`, result of `+1` only means that `n` is
 * a quadratic nonresidue for an even number (zero or more)
 * of prime factors of `m`. In order for `n` to be a quadratic
 * residue, it has to be a residue for each prime factor of `m`.
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
 * Square root of `y.v` modulo prime `y.M`
 *
 * @param M - the modulo<I, ...> type
 */
template <typename M>
M sqrt_cipolla(const M& y) {
	M e0 = zeroT<M>::of(y), e1 = identityT<M>::of(y);
	// find a quadratic nonresidue `d` modulo `p`
	M a = e0, d = e0;
	do {
		a += 1, d = a * a - y;
	} while (powT(d, (y.M - 1) / 2) == 1); // jacobi(d, p) == 1
	// r = (a + sqrt(d)) ^ ((p + 1) / 2)
	return powT(quadraticX<M>(a, e1, d), (y.M + 1) / 2).a;
}

/**
 * Square root of `y` modulo prime `p`
 */
template <typename I>
I sqrt_cipolla(const I& y, const I& p) {
	return sqrt_cipolla(moduloX<I>(y, p)).v;
}

/**
 * Square root of `y` modulo prime power `p^k`
 */
template <typename I>
I sqrt_hensel_lift(const I& y, const I& p, I k) {
	typedef moduloX<I> modx;
	// f(r) == r^2 - y; f'(r) == 2r;
	modx r = sqrt_cipolla(modx(y, p));
	for (I i = 1; i < k; i *= 2) {
		I phi = r.M / p * (p - 1); // euler_phi(r.M)
		modx u = powT(r * 2, phi - 1); // f'(r) ^-1
		r.M = (i * 2 < k) ? r.M * r.M : powT(p, k); // lift modulus
		modx v = r * r - y; // f(r)
		r -= v * u;
	}
	return r.v;
}

/**
 * Primitive root modulo `m`
 *
 * `m` must be 2, 4, p^k or 2p^k.
 *
 * @param m - modulus
 * @param phi - `euler_phi(m)`; number of coprimes with `m` up to `m`
 * @param phi_factors - unique prime factors of `phi`
 */
template<typename I>
I primitive_root(I m, I phi, const std::vector<I> &phi_factors) {
	typedef moduloX<I> modx;
	for (I g = 1; g < m; g += 1) {
		if (gcd(g, m) > 1) continue;
		bool primitive = true;
		for (const I& p : phi_factors) {
			primitive &= (powT(modx(g, m), phi / p) != 1);
		}
		if (primitive) return g;
	}
	return 0;
}

/**
 * Primitive root modulo `m`
 *
 * `m` must be 2, 4, p^k or 2p^k.
 */
int primitive_root(int m, prime_holder& prim);

/**
 * `k-th` roots of unity modulo `m`
 *
 * `m` must be 2, 4, p^k or 2p^k.
 *
 * @param m - modulus
 * #param k - k-th root
 * @param lam - `carmichael_lambda(m)`
 * @param phi - `euler_phi(m)`
 * @param phi_factors - unique prime factors of `phi`
 */
template<typename I>
std::vector<I> kth_roots(I m, I k, I lam, I phi, const std::vector<I> &phi_factors) {
	typedef moduloX<I> modx;
	I g = primitive_root(m, phi, phi_factors);
	I l = gcd(lam, k);
	modx w = powT(modx(g, m), lam / l);
	modx r = identityT<modx>::of(w);
	std::vector<I> vr;
	for (I j = 0; j < l; j++) {
		vr.push_back(r.v);
		r *= w;
	}
	sort(vr.begin(), vr.end());
	return vr;
}

/**
 * `k-th` roots of unity modulo `m`
 *
 * `m` must be 2, 4, p^k or 2p^k.
 */
std::vector<int> kth_roots(int m, int k, prime_holder& prim);

/**
 * Builds the look-up table for `factorial_mod_p`
 */
template<typename P>
void factorial_mod_p_table(P p, P* table) {
	table[0] = 1;
	for (P i = 1; i < p; i++) {
		table[i] = modulo_multiply(table[i - 1], i, p);
	}
}

/**
 * Factorial of n modulo p
 *
 * `(n! / p^e) % p`, where `p^e` is the highest power of `p` dividing `n`.
 * Note, before modulo operation gets applied, all factors `p` are removed.
 *
 * Complexity: O(p + log_p n)
 *
 * @param I - type of the number `n`
 * @param P - type of the prime moduli
 * @param n - number to take factorial of
 * @param p - prime moduli
 * @param table - look-up table of `k! % p` for all `k < p`
 * @param e_out - if given, the highest exponent of `p` that divides `n!`
 *                will be stored in it
 */
template<typename I, typename P>
P factorial_mod_p(I n, P p, P* table, I *e_out = nullptr) {
	P r = 1;
	I e = 0;
	while (n > 1) {
		r = modulo_multiply(r, table[n % p], p);
		n /= p;
		e += n;
		if (n % 2 == 1) r = -r;
	}
	modulo_normalize(r, p);
	if (e_out) *e_out = e;
	return r;
}

} // math
} // altruct
