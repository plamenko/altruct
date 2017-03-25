#pragma once

#include "base.h"

#include <vector>

namespace altruct {
namespace math {

/**
 * Calculates the exponent of the prime `p` in `n!`.
 */
template <typename I>
I factorial_prime_exponent(I n, I p) {
	I e = 0;
	while ((n /= p) > 0) e += n;
	return e;
}

/**
 * Calculates the exponent of the prime `p` in `binomial(n, k)`.
 */
template <typename I>
I binomial_prime_exponent(I n, I k, I p) {
	return factorial_prime_exponent(n, p) - factorial_prime_exponent(n - k, p) - factorial_prime_exponent(k, p);
}

/**
 * Calculates the exponent of the prime `p` in `multinomial(k1, ..., kl)`.
 */
template <typename I, typename IT>
I multinomial_prime_exponent(IT k_begin, IT k_end, I p) {
	I e = 0;
	I n = 0;
	for (IT k_it = k_begin; k_it != k_end; ++k_it) {
		e += factorial_prime_exponent(*k_it, p);
		n += *k_it;
	}
	return factorial_prime_exponent(n, p) - e;
}

/**
 * Calculates the exponent of the prime `p` in `multinomial(k1, ..., kl)`.
 */
template <typename I, typename C>
I multinomial_prime_exponent(C k_container, I p) {
	return multinomial_prime_exponent(k_container.begin(), k_container.end(), p);
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
