#pragma once

#include "algorithm/math/base.h"
#include "structure/container/sqrt_map.h"

namespace altruct {
namespace math {

/**
 * Calculates `PrimeSum[n / k]` for each `k` in `[1, n]` in `O(n^(5/7))`.
 *
 * Where:
 *  `PrimeSum[n] := Sum[If[IsPrime[k], k, 0], {k, 1, n}]`
 *
 * Note, there is only `O(sqrt n)` different values and the result is given as `sqrt_map`.
 *
 * @param castT - casts I to T
 */
template<typename T, typename I, typename CAST_T, typename J = int>
container::sqrt_map<I, T> prime_sum_sqrt(I n, T id, CAST_T castT) {
	// Initially, we start with the sum of all numbers:
	// s[i] = Sum[k, {2 <= k <= i}]
	// After each round j, all multiples of a prime p[j] get eliminated:
	// s[i] = Sum[k, {2 <= k <= i, SmallestPrimeFactorOf[k] > p[j] || IsPrime[k]}]
	J q = J(sqrtT(n)) + 1;
	container::sqrt_map<I, T> s(q - 1, n);
	for (J i = 1; i < q; i++) {
		T t = castT(i);
		s[i] = id * t * (t + 1) / 2 - 1;
	}
	for (J k = J(n / q); k >= 1; k--) {
		I i = n / k;
		T t = castT(i);
		s[i] = id * t * (t + 1) / 2 - 1;
	}
	for (J p = 2; p < q; p++) {
		if (s[p - 1] == s[p]) continue;
		T t = s[p - 1];
		I p2 = sqT<I>(p);
		J k_max = J(std::min(n / q, n / p2));
		for (J k = 1; k <= k_max; k++) {
			//I i = n / k;
			//s.el(i) -= (s.el(i / p) - t) * p;
			I j = n / (I(k) * p);
			s.hi(k) -= (s.el(j) - t) * p;
		}
		for (J i = q - 1; i >= p2; i--) {
			s.lo(i) -= (s.lo(i / p) - t) * p;
		}
	}
	return s;
}

/**
 * Calculates `PrimeSum[n]` in `O(n^(5/7))`.
 *
 * @param castT - casts I to T
 */
template<typename T, typename I, typename CAST_T, typename J = int>
T prime_sum(I n, T id, CAST_T castT) {
	return (n < 1) ? zeroT<T>::of(id) : prime_sum_sqrt<T, I, CAST_T, J>(n, id, castT)[n];
}

/**
 * Calculates `PrimePi[n / k]` for each `k` in `[1, n]` in `O(n^(5/7))`.
 *
 * Where:
 *  `PrimePi[n] := Sum[If[IsPrime[k], 1, 0], {k, 1, n}]`
 *
 * Note, there is only `O(sqrt n)` different values and the result is given as `sqrt_map`.
 */
template<typename I = int64_t, typename J = int>
container::sqrt_map<I, I> prime_pi_sqrt(I n) {
	J q = J(sqrtT(n)) + 1;
	container::sqrt_map<I, I> pi(q - 1, n);
	for (J i = 1; i < q; i++) {
		pi[i] = i - 1;
	}
	for (J k = J(n / q); k >= 1; k--) {
		I i = n / k;
		pi[i] = i - 1;
	}
	for (J p = 2; p < q; p++) {
		if (pi[p - 1] == pi[p]) continue;
		I t = pi[p - 1];
		I p2 = sqT<I>(p);
		J k_max = J(std::min(n / q, n / p2));
		for (J k = 1; k <= k_max; k++) {
			//I i = n / k;
			//pi.el(i) -= pi.el(i / p) - t;
			I j = n / (I(k) * p);
			pi.hi(k) -= pi.el(j) - t;
		}
		for (J i = q - 1; i >= p2; i--) {
			pi.lo(i) -= pi.lo(i / p) - t;
		}
	}
	return pi;
}

/**
 * Calculates `PrimePi[n]` in `O(n^(5/7))`.
 */
template<typename I = int64_t, typename J = int>
I prime_pi(I n) {
	return (n < 1) ? 0 : prime_pi_sqrt<I, J>(n)[n];
}

} // math
} // altruct
