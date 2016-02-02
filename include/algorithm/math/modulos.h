#pragma once

#include "structure/math/modulo.h"
#include "structure/math/quadratic.h"
#include "structure/math/prime_holder.h"

#include <vector>
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
 * Square root of `y` modulo `p`
 */
template <typename I, int ID>
modulo<I, ID> sqrt_cipolla(modulo<I, ID> y) {
	typedef modulo<I, ID> mod;
	const I& p = mod::M;
	// find a quadratic nonresidue `d` modulo `p`
	mod a = 0, d = 0;
	do {
		a += 1;
		d = a * a - y;
	} while (powT(d, (p - 1) / 2) == 1); // jacobi(d, p) == 1
	// r = (a + sqrt(d)) ^ ((p + 1) / 2)
	typedef quadratic<mod, 0> quad;
	quad::D = d;
	quad r = powT(quad(a, 1), (p + 1) / 2);
	return r.a;
}

/**
 * Square root of `y` modulo `p`
 */
template <typename I>
I sqrt_cipolla(I y, I p) {
	typedef modulo<I, 1> mod;
	mod::M = p;
	return sqrt_cipolla(mod(y)).v;
}

/**
 * Square root of `y` modulo `p^k`
 */
template <typename I>
I sqrt_hensel_lift(I y, int p, int k) {
	typedef modulo<I, 1> mod;
	mod::M = p; int p_i = p, p_k = powT(p, k);
	// f(r) == r^2 - y; f'(r) == 2r;
	mod r = sqrt_cipolla(mod(y));
	for (int i = 1; i < k; i *= 2) {
		int phi = p_i / p * (p - 1); // euler_phi(p_i)
		mod u = powT(r * 2, -1 + phi); // f'(r) ^-1
		mod::M = p_i = (i * 2 < k) ? p_i * p_i : p_k; // lift modulus
		mod v = r * r - y; // f(r)
		r -= u * v;
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
	typedef modulo<I, 1> mod;
	mod::M = m;
	for (I g = 1; g < m; g += 1) {
		if (gcd(g, m) > 1) continue;
		bool primitive = true;
		for (const I& p : phi_factors) {
			primitive &= (powT<mod>(g, phi / p) != 1);
		}
		if (primitive) return g;
	}
	return 0;
}

/**
 * Primitive root modulo `m`
 */
int primitive_root(int m, prime_holder& prim) {
	int phi = euler_phi(prim.factor_integer(m));
	auto phi_factors = prime_factors(prim.factor_integer(phi));
	return primitive_root(m, phi, phi_factors);
}

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
std::vector<I> kth_roots(I m, int k, I lam, I phi, const std::vector<I> &phi_factors) {
	typedef modulo<I, 1> mod;
	mod::M = m;
	I g = primitive_root(m, phi, phi_factors);
	I l = gcd(lam, k);
	mod w = powT<mod>(g, lam / l);
	mod r = 1;
	std::vector<int> vr;
	for (int j = 0; j < l; j++) {
		vr.push_back(r.v);
		r *= w;
	}
	sort(vr.begin(), vr.end());
	return vr;
}

/**
 * `k-th` roots of unity modulo `m`
 */
std::vector<int> kth_roots(int m, int k, prime_holder& prim) {
	auto vf = prim.factor_integer(m);
	int lam = carmichael_lambda(vf);
	int phi = euler_phi(vf);
	auto phi_factors = prime_factors(prim.factor_integer(phi));
	return kth_roots(m, k, lam, phi, phi_factors);
}


/**
 * Builds the look-up table for `factorial_mod_p`
 */
template<typename M>
void factorial_mod_p_table(M p, M* table) {
	table[0] = 1;
	for (M i = 1; i < p; i++) {
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
 * @param M - type of the prime moduli
 * @param n - number to take factorial of
 * @param p - prime moduli
 * @param table - look-up table of `k! % p` for all `k < p`
 * @param e_out - if given, the highest exponent of `p` that divides `n!`
 *                will be stored in it
 */
template<typename I, typename M>
M factorial_mod_p(I n, M p, M* table, I *e_out = nullptr) {
	M r = 1;
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
