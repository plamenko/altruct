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

} // math
} // altruct
