#pragma once

#include "base.h"

#include <vector>

namespace altruct {
namespace math {

/**
 * Calculates the exponent of the prime `p` in `n!`.
 */
template <typename I>
I factorial_prime_exponent(I p, I n) {
	I e = 0;
	while ((n /= p) > 0) e += n;
	return e;
}

/**
 * Calculates the exponent of the prime `p` in `binomial(n, k)`.
 */
template <typename I>
I binomial_prime_exponent(I p, I n, I k) {
	return factorial_prime_exponent(p, n) - factorial_prime_exponent(p, n - k) - factorial_prime_exponent(p, k);
}

/**
 * Calculates the exponent of the prime `p` in `multinomial(k1, ..., kl)`.
 */
template <typename I, typename IT>
I multinomial_prime_exponent(I p, IT k_begin, IT k_end) {
	I e = 0;
	I n = 0;
	for (IT k_it = k_begin; k_it != k_end; ++k_it) {
		e += factorial_prime_exponent(p, *k_it);
		n += *k_it;
	}
	return factorial_prime_exponent(p, n) - e;
}

/**
 * Calculates the exponent of the prime `p` in `multinomial(k1, ..., kl)`.
 */
template <typename I, typename C>
I multinomial_prime_exponent(I p, C k_container) {
	return multinomial_prime_exponent(p, k_container.begin(), k_container.end());
}

/**
 * Calculates the multinomial coefficient based on the elements.
 *
 * Note, the elements must be in sorted order.
 *
 * E.g. for elements {a, a, b, b, b, c}, this calculates
 * `multinomial(2, 3, 1) = (2+3+1)!/(2!3!1!)`.
 */
template<typename T, typename It>
T elements_multinomial(It begin, It end, T id = T(1)) {
	if (begin == end) return id;
	T n = id, l = id;
	T f = id, d = id;
	for (It prev = begin++; begin != end; ++prev, ++begin) {
		n += id;
		l = (*prev == *begin) ? l + id : id;
		f *= n;
		d *= l;
	}
	return f / d;
}

/**
 * Calculates the multinomial coefficient based on the elements.
 */
template <typename T, typename C>
T elements_multinomial(C container, T id = T(1)) {
	return elements_multinomial(container.begin(), container.end(), id);
}

} // math
} // altruct
