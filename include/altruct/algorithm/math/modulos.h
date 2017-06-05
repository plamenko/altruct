#pragma once

#include "altruct/structure/math/modulo.h"
#include "altruct/structure/math/quadratic.h"
#include "altruct/structure/math/prime_holder.h"

#include <set>
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
void chinese_remainder(T* a, T* n, T a1, T n1, T a2, T n2) {
	T e0 = zeroT<T>::of(*a);
	T ni1, ni2; T g = gcd_ex(n1, n2, &ni1, &ni2);
	if ((a2 - a1) % g != e0) { *n = e0, *a = e0; return; }
	modulo_normalize(&ni2, n1);
	modulo_normalize(&ni1, n2);
	T t1 = modulo_mul(a1, ni2, n1);
	T t2 = modulo_mul(a2, ni1, n2);
	n1 /= g; n2 /= g; *n = n1 * n2 * g;
	*a = modulo_mul(t1, n2, *n) + modulo_mul(t2, n1, *n);
	return modulo_normalize(a, *n);
}
template<typename T>
void chinese_remainder(T* ar, T* nr, T a, T n) {
	chinese_remainder(ar, nr, *ar, *nr, a, n);
}
template<typename T>
T chinese_remainder(T a1, T n1, T a2, T n2) {
	T a, n; chinese_remainder(&a, &n, a1, n1, a2, n2);
	return a;
}

/**
 * Calculates coefficients of the mixed radix equation using Garner Algorithm.
 *
 * Let `u` be given as a set of congruence equations: `u == a[i] (mod p[i])`.
 * Let `q[i]` be the product of the first `i` elements of `p`. `q[0] = 1`.
 * Then `a` can be represented in the mixed radix form: `u = Sum[x[i] * q[i]]`.
 * Note: `u` is unique modulo `q[n]` as per Chinese Remainder Theorem.
 *
 * @param vap - vector of (remainder, modulus) pairs, moduli shoulde be pairwise relatively prime
 * @return vx - vector of (coeeficient, modulus) pairs, moduli are same as in `vap`
 */
template<typename V>
V garner(const V& vap) {
	V vx(vap.size());
	for (int i = 0; i < vx.size(); i++) {
		auto y = vap[i];
		for (int j = 0; j < i; j++) {
			y -= vx[j];
			y /= vap[j].M();
		}
		vx[i] = y;
	}
	return vx;
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
 * Square root of `y.v` modulo odd prime `y.M`
 *
 * @param M - the modulo<I, ...> type
 * @param y - integer such that: y = x^2 (mod p), (y, p) = 1
 */
template <typename M>
M sqrt_cipolla(const M& y) {
	M e0 = zeroT<M>::of(y), e1 = identityT<M>::of(y);
	// find a quadratic nonresidue `d` modulo `p`
	M a = e0, d = e0;
	do {
		a += 1, d = a * a - y;
	} while (powT(d, (y.M() - 1) / 2) == 1); // jacobi(d, p) == 1
	// r = (a + sqrt(d)) ^ ((p + 1) / 2)
	return powT(quadraticX<M>(a, e1, d), (y.M() + 1) / 2).a;
}

/**
 * Square root of `y` modulo odd prime `p`
 *
 * @param y - integer such that: y = x^2 (mod p), (y, p) = 1
 */
template <typename I>
I sqrt_cipolla(const I& y, const I& p) {
	return sqrt_cipolla(moduloX<I>(y, p)).v;
}

/**
 * Square root of `y` modulo odd prime power `p^k`
 *
 * @param y - integer such that: y = x^2 (mod p^k), (y, p) = 1
 */
template <typename I>
I sqrt_hensel_lift(const I& y, const I& p, I k) {
	typedef moduloX<I> modx;
	// f(r) == r^2 - y; f'(r) == 2r;
	modx r = sqrt_cipolla(modx(y, p));
	for (I i = 1; i < k; i *= 2) {
		I phi = r.M() / p * (p - 1); // euler_phi(r.M)
		modx u = powT(r * 2, phi - 1); // f'(r) ^-1
		r.M() = (i * 2 < k) ? r.M() * r.M() : powT(p, k); // lift modulus
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
 * @return - root or 0 if there is none
 */
template<typename I>
I primitive_root(I m, I phi, const std::vector<I> &phi_factors) {
	typedef moduloX<I> modx;
	for (I g = 1; g < m; g += 1) {
		if (gcd(g, m) > 1) continue;
		bool primitive = true;
		for (const I& p : phi_factors) {
			primitive &= (powT(modx(g, m), I(phi / p)) != 1);
		}
		if (primitive) return g;
	}
	return 0;
}

/**
 * Primitive root of unity modulo `m`
 *
 * @param m - modulus
 * @param lam - `carmichael_lambda(m)`;
 * @param lam_factors - unique prime factors of `lam`
 * @return - root or 0 if there is none
 */
template<typename I>
I primitive_root_of_unity(I m, I lam, const std::vector<I> &lam_factors) {
	return primitive_root(m, lam, lam_factors);
}

/**
 * Primitive root modulo `m`
 *
 * `m` must be 2, 4, p^k or 2p^k.
 */
int primitive_root(int m, prime_holder& prim);

/**
 * Primitive root of unity modulo `m`
 */
int primitive_root_of_unity(int m, prime_holder& prim);

/**
 * `k-th` roots of unity modulo `m`
 *
 * `m` must be 2, 4, p^k or 2p^k.
 *
 * @param m - modulus
 * @param k - k-th root
 * @param lam - `carmichael_lambda(m)`
 * @param g - primitive root of unity modulo `m`
 */
template<typename I>
std::set<I> kth_roots_of_unity(I m, I k, I lam, I g) {
	typedef moduloX<I> modx;
	I d = gcd(k, lam);
	modx w = powT(modx(g, m), I(lam / d));
	modx r = identityT<modx>::of(w);
	std::set<I> sr;
	for (I j = 0; j < d; j++) {
		sr.insert(r.v);
		r *= w;
	}
	return sr;
}

/**
 * `k-th` roots of unity modulo `m`
 */
std::set<int> kth_roots_of_unity(int m, int k, prime_holder& prim);

/**
* `k-th` roots of `n` modulo `m`
*
* `m` must be 2, 4, p^k or 2p^k.
*
* @param m - modulus
* @param k - k-th root
* @param phi - `euler_phi(m)`
* @param g - primitive root modulo `m`
* @param l - `l` such that `g ^ l == n  (mod m)`
*/
template<typename I>
std::set<I> kth_roots(I m, I k, I phi, I g, I l) {
	typedef moduloX<I> modx;
	I d = gcd(k, phi);
	if (d == 0 || l % d != 0) return{};
	phi /= d; l /= d; k /= d;
	// g^(l/k) == n^(1/k)
	modx h = modx(l, phi) / k;
	modx r = powT(modx(g, m), h.v);
	modx w = powT(modx(g, m), phi);
	std::set<I> sr;
	for (int i = 0; i < d; i++) {
		sr.insert(r.v);
		r *= w;
	}
	return sr;
}

} // math
} // altruct
