#pragma once

#include "base.h"
#include <unordered_map>

namespace altruct {
namespace math {

/**
 * PrimePi - number of primes up to `m` in `O(m^(2/3))`.
 * Meissel-Lehmer algorithm;  O(m^(2/3)) space.
 *
 * @param pi - array of prime_pi values up to `m^(2/3)`
 * @param p - array of primes up to `m^(1/2)`
 * @param pi_size - value up to which prime_pi is precomputed
 *
 * @deprecated - `prime_pi_sqrt` implementation is simpler and faster in practice
 *               this is because this implementation requires too much memory
 */
template<typename I, typename J>
I prime_pi_PHI_key(I m, J n) {
	// this works for I=i64, m < 2^48 (10^14), n < 2^16
	// this works for I=i32, m < 2^24 (10^7), n < 2^8
	return (m << 16) + n;
}
template<typename I, typename J>
std::unordered_map<I, I>& prime_pi_PHI_tbl() {
	static std::unordered_map<I, I> tbl;
	return tbl;
}
template<typename I, typename J>
I prime_pi_PHI(I m, J n, const int* p) {
	if (m == 0 || n == 0) return m;
	auto t = prime_pi_PHI_key(m, n);
	auto& tbl = prime_pi_PHI_tbl<I, J>();
	if (tbl.count(t)) return tbl[t];
	I r = prime_pi_PHI(m, n - 1, p) - prime_pi_PHI(m / p[n - 1], n - 1, p);
	return tbl[t] = r;
}
template<typename I, typename J>
I prime_pi_P2(I m, J n, const J* pi, const int* p) {
	I r = 0;
	for (J k = n; (I)p[k] * p[k] <= m; k++)
		r += pi[m / p[k]] - k;
	return r;
}
template<typename I>
std::unordered_map<I, I>& prime_pi_tbl() {
	static std::unordered_map<I, I> tbl;
	return tbl;
}
template<typename I, typename J>
I prime_pi_deprecated(I m, const J* pi, const int* p, int pi_size = 0) {
	if (m < 2) return 0;
	if (m < pi_size) return pi[m];
	auto& tbl = prime_pi_tbl<I>();
	if (tbl.count(m)) return tbl[m];
	J y = (J)cbrtT(m) + 1;
	J n = pi[y];
	return tbl[m] = prime_pi_PHI(m, n, p) - prime_pi_P2(m, n, pi, p) + n - 1;
}

} // math
} // altruct
