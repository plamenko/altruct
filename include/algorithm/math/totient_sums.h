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
 * A helper function for `calc_sum_phi_D_L`.
 */
template<typename T, typename CAST_T>
std::vector<T> sum_g_L(const polynom<T>& g, int L, const std::vector<int64_t>& vn, T id, CAST_T castT) {
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

	// preprocess `U = n^(2/3)` values of `Sum[p(k) * f[k], {k, 1, U}]`
	int64_t n = *std::max_element(vn.begin(), vn.end());
	int U = (int)isq(icbrt(n));
	altruct::container::sqrt_map<int64_t, T> msf(U, n);
	moebius_transform(msf, U, _g);
	for (int k = 1; k < U; k++) {
		msf[k] = _p(k) * msf[k] + msf[k - 1];
	}

	// calculate the values of interest
	std::vector<T> v;
	for (auto k : vn) {
		msf.reset_max(k);
		v.push_back(sum_m<T>(k, _t, _s, msf));
	}
	return v;
}

/**
 * Calculates `Sum[k^L euler_phi_D(k), { k, 1, n }]` in `O(n^(2/3))`.
 *
 * Where:
 *   euler_phi_D(n) = Sum[[GCD(a1, a2, ..., aD, n) == 1], {a1 <= a2 <= ... <= aD <= n}]
 *                  = Sum[mu(n/d) * binomial(D + d - 1, D), {d|x}]
 *   euler_phi_0(n) = [n == 1]
 *   euler_phi_1(n) = euler_totient(n)
 *   euler_phi_2(n) = euler_totient_2d(n)
 *   euler_phi_3(n) = euler_totient_3d(n)
 *   ...
 */
template<typename T, typename CAST_T>
std::vector<T> sum_phi_D_L(int D, int L, const std::vector<int64_t>& vn, T id, CAST_T castT) {
	typedef polynom<T> poly;
	T e0 = zeroT<T>::of(id);
	// g_phi_D(n) = binomial(D + n - 1, D)
	poly g_phi_D{ id };
	for (int i = 0; i < D; i++) {
		g_phi_D *= poly{ id * i, id } / (id * (i + 1));
	}
	return sum_g_L(g_phi_D, L, vn, id, castT);
}
template<typename T, typename CAST_T>
T sum_phi_D_L(int D, int L, int64_t n, T id, CAST_T castT) {
	return sum_phi_D_L(D, L, std::vector<int64_t>{ n }, id, castT).back();
}

} // math
} // altruct
