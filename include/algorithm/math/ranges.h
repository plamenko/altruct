#pragma once

#include <iterator>
#include <type_traits>

#include "algorithm/math/base.h"

namespace altruct {
namespace math {

/**
 * Builds range look-up table up to `n`.
 */
template<typename It, typename T = typename std::iterator_traits<It>::value_type>
void range(It begin, It end, T step = T(1)) {
	T v = zeroT<T>::of(step);
	for (It it = begin; it != end; ++it) {
		*it = v; v += step;
	}
}

/**
 * Builds the powers of `base` look-up table up to `n`.
 */
template<typename It, typename T = typename std::iterator_traits<It>::value_type>
void powers(It begin, It end, T base) {
	T v = identityT<T>::of(base);
	for (It it = begin; it != end; ++it) {
		*it = v; v *= base;
	}
}

/**
 * Builds the factorial look-up table up to `n`.
 */
template<typename It, typename T = typename std::iterator_traits<It>::value_type>
void factorials(It begin, It end, T id = T(1)) {
	T v = id, i = id;
	for (It it = begin; it != end; ++it) {
		*it = v; v *= i; i += id;
	}
}

/**
 * Builds the inverse factorial look-up table up to `n`.
 */
template<typename It, typename T = typename std::iterator_traits<It>::value_type>
void inv_factorials(It begin, It end, T id = T(1)) {
	T fact = id, i = id;
	for (It it = begin; it != end; ++it) {
		fact *= i; i += id;
	}
	T ifact = id / fact;
	for (It it = end; it != begin; ) {
		i -= id; ifact *= i; *--it = ifact;
	}
}

/**
 * Inverts the table elements with respect to multiplication.
 * Zero values are left zeros.
 */
template<typename It, typename T = typename std::iterator_traits<It>::value_type>
void invert(It begin, It end, T id = T(1)) {
	T e0 = zeroT<T>::of(id);
	for (It it = begin; it != end; ++it) {
		if (*it != e0) *it = id / *it;
	}
}

/**
 * Negates the table elements.
 */
template<typename It, typename T = typename std::iterator_traits<It>::value_type>
void negate(It begin, It end) {
	for (It it = begin; it != end; ++it) {
		*it = -*it;
	}
}

/**
 * Alternates the sign of the table elements.
 * I.e. negates odd elements.
 */
template<typename It, typename T = typename std::iterator_traits<It>::value_type>
void alternate(It begin, It end) {
	int s = 1;
	for (It it = begin; it != end; ++it, s = -s) {
		if (s < 0) *it = -*it;
	}
}

/**
 * Accumulates the table elements.
 */
template<typename It, typename T = typename std::iterator_traits<It>::value_type>
void accumulate(It begin, It end) {
	if (begin == end) return;
	It prev = begin++;
	for (It it = begin; it != end; ++it) {
		*it += *prev++;
	}
}

/**
 * Differences between the table elements.
 */
template<typename It, typename T = typename std::iterator_traits<It>::value_type>
void differences(It begin, It end) {
	if (begin == end) return;
	It prev = --end;
	for (It it = end; it != begin; --it) {
		*it -= *--prev;
	}
}

} // math
} // altruct
