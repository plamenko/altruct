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
 * Dirichlet convolution of `f` and `g` up to `n` in `O(n log n)`.
 *
 * Calculates `h` where `h[n] = Sum[f(n/d) * g(d), {d|n}]`
 *
 * Where:
 *   `f` and `g` are arbitrary arithmetic functions
 *
 * @param h - table to store the result
 * @param f, g - functions as defined above
 */
template<typename T, typename F1, typename F2, typename TBL>
void dirichlet_convolution(TBL& h, F1 f, F2 g, int n) {
	T e0 = zeroOf(f(0));
	for (int i = 0; i < n; i++) {
		h[i] = e0;
	}
	for (int d = 1; d < n; d++) {
		for (int e = 1, i = d; i < n; i += d, e++) {
			h[i] += f(d) * g(e);
		}
	}
}

/**
 * Dirichlet inverse of `f` up to `n` in `O(n log n)`.
 *
 * Calculates `f_inv` such that: `f * f_inv = e`
 *
 * Where:
 *   `f` is an arbitrary arithmetic function such that:
 *       f(1) != 0, and f(1) is invertible
 *   `e` is the dirichlet multiplicative identity: `e(n) = [n == 1]`
 *
 * @param f_inv - table to store the result
 * @param f - function as defined above
 */
template<typename T, typename F1, typename TBL>
void dirichlet_inverse(TBL& f_inv, F1 f, int n) {
	T f1 = f(1);
	T e0 = zeroOf(f1), e1 = identityOf(f1);
	T if1 = e1 / f1;
	for (int i = 0; i < n; i++) {
		f_inv[i] = e0;
	}
	f_inv[1] = e1;
	for (int d = 1; d < n; d++) {
		f_inv[d] *= if1;
		for (int j = 2, i = d * 2; i < n; i += d, j++) {
			f_inv[i] -= f(j) * f_inv[d];
		}
	}
}

/**
 * Calculates all the values of a multiplicative function `f` up to `n`,
 * from the values at prime powers, in `O(n log log n)`.
 *
 * @param f - table of values of `f`
 *            must be set to the actual value for prime powers and 1 elsewhere
 * @param pf - table of prime factors up to `n`; pf[k] = some_prime_factor_of(k)
 */
template<typename TBL>
void calc_multiplicative(TBL& f, int n, int* pf) {
	auto e1 = f[1];
	for (int p = 2; p < n; p++) {
		if (pf[p] != p) continue; // not a prime
		for (int64_t qq = p; qq < n; qq *= p) {
			for (int q = (int)qq, l = 2, m = 2 * q; m < n; m += q, l++) {
				if (l % p != 0) f[m] *= f[q];
			}
		}
	}
}

/**
 * Dirichlet convolution of `f` and `g` up to `n` in `O(n log log n)`.
 *
 * Calculates `h` where `h[n] = Sum[f(n/d) * g(d), {d|n}]`
 *
 * Where:
 *   `f` and `g` are arbitrary arithmetic functions such that
 *   `h = f * g` is a multiplicative function
 *
 * Note that only `h` needs to be multiplicative!
 *
 * @param h - table to store the result
 * @param f, g - functions as defined above
 * @param pf - table of prime factors up to `n`; pf[k] = some_prime_factor_of(k)
 */
template<typename T, typename F1, typename F2, typename TBL>
void dirichlet_convolution_multiplicative(TBL& h, F1 f, F2 g, int n, int *pf) {
	auto e1 = identityOf(f(1)), e0 = zeroOf(f(1));
	for (int i = 1; i < n; i++) {
		h[i] = e1;
	}
	int q[32]; T fq[32], gq[32];
	for (int p = 2; p < n; p++) {
		if (pf[p] != p) continue; // not a prime
		int m = 0;
		for (int64_t qq = 1; qq < n; qq *= p) {
			fq[m] = f((int)qq);
			gq[m] = g((int)qq);
			q[m++] = (int)qq;
		}
		for (int k = 0; k < m; k++) {
			auto hq_k = e0;
			for (int i = 0; i <= k; i++) {
				hq_k += fq[k - i] * gq[i];
			}
			h[q[k]] = hq_k;
		}
	}
	calc_multiplicative(h, n, pf);
}

/**
 * Dirichlet inverse of `f` up to `n` in `O(n log log n)`.
 *
 * Calculates `f_inv` such that: `f * f_inv = e`
 *
 * Where:
 *   `f` is a multiplicative function, which means:
 *       its inverse `f_inv` is also multiplicative
 *   `e` is the dirichlet multiplicative identity: `e(n) = [n == 1]`
 *
 * @param f_inv - table to store the result
 * @param f - function as defined above
 * @param pf - table of prime factors up to `n`; pf[k] = some_prime_factor_of(k)
 */
template<typename T, typename F1, typename TBL>
void dirichlet_inverse_multiplicative(TBL& f_inv, F1 f, int n, int* pf) {
	auto e1 = f(1), e0 = zeroOf(f(1));
	for (int i = 1; i < n; i++) {
		f_inv[i] = e1;
	}
	int q[32]; T fq[32], hq[32];
	for (int p = 2; p < n; p++) {
		if (pf[p] != p) continue; // not a prime
		int m = 0;
		for (int64_t qq = 1; qq < n; qq *= p) {
			fq[m] = f((int)qq);
			q[m++] = (int)qq;
		}
		hq[0] = e1;
		for (int k = 1; k < m; k++) {
			hq[k] = e0;
			for (int i = 0; i < k; i++) {
				hq[k] -= fq[k - i] * hq[i];
			}
			f_inv[q[k]] = hq[k];
		}
	}
	calc_multiplicative(f_inv, n, pf);
}

/**
 * Calculates all the values of a completely multiplicative function `f` up to `n`,
 * from the values at primes, in `O(n)`.
 *
 * param f - table of values of `f`
 *           must be set to the actual value for primes and 1 elsewhere
 * @param pf - table of prime factors up to `n`; pf[k] = some_prime_factor_of(k)
 */
template<typename TBL>
void calc_completely_multiplicative(TBL& f, int n, int* pf) {
	auto e1 = f[1];
	for (int i = 2; i < n; i++) {
		int p = pf[i];
		if (pf[i] != i) f[i] = f[i / p] * f[p];
	}
}

/**
 * Dirichlet convolution of `f` and `g` up to `n` in `O(n)`.
 *
 * Calculates `h` where `h[n] = Sum[f(n/d) * g(d), {d|n}]`
 *
 * Where:
 *   `f` and `g` are arbitrary arithmetic functions such that
 *   `h = f * g` is a completely multiplicative function
 *
 * Note that only `h` needs to be completely multiplicative!
 * For example, `f(n) = mu(n)` and `g(n) = sigma1(n)` are not
 * completely multiplicative, but its convolution `h(n) = n` is!
 *
 * @param h - table to store the result
 * @param f, g - functions as defined above
 * @param pf - table of prime factors up to `n`; pf[k] = some_prime_factor_of(k)
 */
template<typename T, typename F1, typename F2, typename TBL>
void dirichlet_convolution_completely_multiplicative(TBL& h, F1 f, F2 g, int n, int *pf) {
	auto e1 = identityOf(f(1)), f1 = f(1), g1 = g(1);
	for (int i = 1; i < n; i++) {
		h[i] = e1;
	}
	for (int p = 2; p < n; p++) {
		if (pf[p] == p) h[p] = f(p) * g1 + g(p) * f1;
	}
	calc_completely_multiplicative(h, n, pf);
}

/**
 * Dirichlet inverse of `f` up to `n` in `O(n)`.
 *
 * Calculates `f_inv` such that: `f * f_inv = e`
 *
 * Where:
 *   `f` is an arbitrary arithmetic function such that
 *   `f_inv = f^-1` is a completely multiplicative function
 *   `e` is the dirichlet multiplicative identity: `e(n) = [n == 1]`
 *
 * Note that only `f_inv` needs to be completely multiplicative!
 * For example, `f(n) = n mu(n)` is not completely multiplicative,
 * but its inverse `f_inv(n) = n` is!
 *
 * @param f_inv - table to store the result
 * @param f - function as defined above
 * @param pf - table of prime factors up to `n`; pf[k] = some_prime_factor_of(k)
 */
template<typename T, typename F1, typename TBL>
void dirichlet_inverse_completely_multiplicative(TBL& f_inv, F1 f, int n, int* pf) {
	auto e1 = f(1);
	for (int i = 1; i < n; i++) {
		f_inv[i] = e1;
	}
	for (int p = 2; p < n; p++) {
		if (pf[p] == p) f_inv[p] = -f(p);
	}
	calc_completely_multiplicative(f_inv, n, pf);
}

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
 * @param M - table to store the calculated values
 * @param t, p - functions as defined above
 * @param n - bound up to which to sieve
 */
template<typename T, typename F1, typename F2, typename TBL>
void sieve_m(TBL& M, F1 t, F2 p, int n) {
	T p1 = p(1);
	T ip1 = identityOf(p1) / p1;
	M[1] = t(1);
	for (int i = 2; i < n; i++) {
		M[i] = t(i) - t(i - 1);
	}
	for (int d = 1; d < n; d++) {
		M[d] *= ip1;
		for (int j = 2, i = d * 2; i < n; i += d, j++) {
			M[i] -= p(j) * M[d];
		}
		if (d > 1) M[d] += M[d - 1];
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
template<typename T, typename F1, typename TBL>
void sieve_m(TBL& M, F1 t, int n) {
	M[1] = t(1);
	for (int i = 2; i < n; i++) {
		M[i] = t(i) - t(i - 1);
	}
	for (int d = 1; d < n; d++) {
		for (int i = d * 2; i < n; i += d) {
			M[i] -= M[d];
		}
		if (d > 1) M[d] += M[d - 1];
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
 * @param t, s - functions as defined above
 * @param n - argument at which to evaluate `M`
 * @param tbl - table to store the calculated values
 */
template<typename T, typename I, typename F1, typename F2, typename TBL>
T sum_m(F1 t, F2 s, I n, TBL& tbl) {
	if (n < 1) return zeroT<T>::of(t(1));
	if (tbl.count(n)) return tbl[n];
	T r = t(n), p1 = s(1) - s(0);
	I q = sqrtT(n);
	for (I k = 2; k <= n / q; k++) {
		r -= (s(k) - s(k - 1)) * sum_m<T>(t, s, n / k, tbl);
	}
	for (I m = 1; m < q; m++) {
		r -= (s(n / m) - s(n / (m + 1))) * sum_m<T>(t, s, m, tbl);
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
 * Same as `sum_m(t, s, n, tbl)` with `p(n) = 1`, `s(n) = n`.
 */
template<typename T, typename I, typename F, typename TBL>
T sum_m(F t, I n, TBL& tbl) {
	if (n < 1) return zeroT<T>::of(t(1));
	if (tbl.count(n)) return tbl[n];
	T r = t(n);
	I q = sqrtT(n);
	for (I k = 2; k <= n / q; k++) {
		r -= sum_m<T>(t, n / k, tbl);
	}
	for (I m = 1; m < q; m++) {
		r -= T((n / m) - (n / (m + 1))) * sum_m<T>(t, m, tbl);
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
template<typename T, typename I, typename TBL>
T mertens(I n, TBL& tbl, T id = T(1)) {
	// p = 1, f = mu, g = delta, t = 1
	return sum_m<T>([&](I k){ return id; }, n, tbl);
}

/**
 * A helper function for `sum_phi_D_L`.
 *
 * Denote `*` as dirichlet convolution and `.` as pointwise multiplication.
 *
 * As defined in `sum_phi_D_L`, we have:
 *   phi_D = mu * g_D
 *   g_D = 1 * phi_D
 *
 * Since a completely multiplicative function (let's call it `p`)
 * distributes pointwise multiplication over Dirichlet convolution,
 * we also have:
 *   p . g_D = p . (1 * phi_D)
 *   p . g_D = (p . 1) * (p . phi_D)
 *   p . g_D = p * (p . phi_D)
 *
 * Let's define:
 *   t' = p . g_D
 *   M' = p . phi_D
 * Then by substitution:
 *   t' = p * M'
 *   t'(n) = Sum[p(d) M'(n/d), {d|n}]
 * Which is equivalent to:
 *   t(n) = Sum[p(k) M(n/k), {k, 1, n}]
 *   which can be calculated efficiently with `sum_m`
 *
 * The above two are equivalent because:
 *   t(n) = Sum[p(k) M(n/k), {k, 1, n}]
 *   t(n) - t(n-1) = Sum[p(k) (M(n/k) - M((n-1)/k)), {k, 1, n}]
 *   But `M(n/k) - M((n-1)/k)` is precisely 1 when `k | n` and 0 otherwise
 *   t(n) - t(n-1) = Sum[p(k) (M(n/k) - M((n-1)/k)), {k|n}]
 *   t(n) - t(n-1) = Sum[p(k) (M(n/k) - M(n/k-1)), {k|n}]
 *   t'(n) = Sum[p(k) M'(n/k), {k|n}]
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

	// preprocess `phi_D = mu * g_D` up to `U`
	int64_t n = *std::max_element(vn.begin(), vn.end());
	if (U <= 0) U = (int)isq(icbrt(n));
	altruct::container::sqrt_map<int64_t, T> mm(U, n);
	moebius_transform(mm, U, _g);
	// preprocess `Sum[p(k) * phi_D[k], {k, 1, n}]` up to `U`
	for (int k = 1; k < U; k++) {
		mm[k] = mm[k - 1] + _p(k) * mm[k];
	}

	// calculate the values of interest with `sum_m`
	std::vector<T> v;
	for (auto k : vn) {
		mm.reset_max(k);
		v.push_back(sum_m<T>(_t, _s, k, mm));
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
