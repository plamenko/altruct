#pragma once

#include <cstdint>
#include <limits>
#include <functional>

namespace altruct {
namespace random {

/**
 * Maps unsigned integer to a double in the `[0, 1]` range, both inclusive.
 *
 * The value is uniformly mapped to the `[0, 1]` range.
 */
template<typename U>
double integer_to_double_0_1(U val) {
	return val / (double)std::numeric_limits<U>::max();
}

/**
 * Maps unsigned integer to an integer in the `[min, max]` range, both inclusive.
 *
 * The value is almost uniformly mapped to the given range.
 * In case the range of `val` is not a proper multiple of the range `[min, max]`,
 * the lower values will be hit more often. See `next_uniform` to avoid this.
 */
template<typename U>
U integer_to_range(U val, U min, U max) {
	U width = max - min + 1;
	return (width == 0) ? val : min + val % width;
}

/**
 * Biggest multiple of `width` that is less than the unsigned integer size.
 *
 * Note: integer size being `2^L`, where `L` is the number of bits.
 * E.g. integer size for uint64_t type is 2^64.
 */
template<typename U>
U biggest_multiple(U width) {
	// 2^L - ((2^L - w) % w)
	return (width == 0) ? 0 : -(U(-width) % width);
}

/**
 * Uniformly selects an integer from the `[min, max]` range, both inclusive.
 *
 * Note, this assumes `next` uniformly provides an unsigned integer from the
 * whole range of `U`.
 *
 * When using modulo operation to reduce an integer value to a range of a
 * certain width, in case the width does not divide the integer range, there
 * will be a slight bias towards the lower numbers of the target range. The
 * bigger the width, the stronger the bias.
 * This bias can be avoided by accepting only values smaller than some multiple
 * of width. We choose the bigest multiple of width that fits the integer size.
 * This multiple is always bigger than half of the integer range, which means
 * that in the worst case, there is a less than 50% chance of not getting a
 * random value smaller than it. This in turn means that the expected number of
 * iterations is only 2. To avoid theoretical infinite loop however, at most
 * 20 iterations are performed. This means that there is still a small chance
 * (2^-20) of a biased selection. This seems like a reasonable trade-off.
 */
template<typename U>
U uniform_next(std::function<U(void)> next, U min, U max) {
	U multiple = biggest_multiple<U>(max - min + 1);
	U val;
	int iter = 0;
	do {
		val = next();
	} while (multiple != 0 && val >= multiple && ++iter < 20);
	return integer_to_range<U>(val, min, max);
}

} // random
} // altruct
