#include "algorithm/math/base.h"

#include <cmath>
#include <stdint.h>

namespace altruct {
namespace math {

// integer square & square root

int64_t isq(int64_t x) {
	return x * x;
}

int32_t isqrt(int64_t x) {
	if (x < 0) return -isqrtc(-x);
	int32_t q = int32_t(floor(sqrt(double(x))));
	while (isq(q) > x) q--;
	while (isq(q + 1) <= x) q++;
	return q;
}

int32_t isqrtc(int64_t x) {
	if (x < 0) return -isqrt(-x);
	int32_t q = int32_t(ceil(sqrt(double(x))));
	while (isq(q) < x) q++;
	while (q > 0 && isq(q - 1) >= x) q--;
	return q;
}

bool is_square(int64_t x) {
	return (isq(isqrt(x)) == x);
}

// integer cube & cube root

int64_t icb(int64_t x) {
	return x * x * x;
}

int32_t icbrt(int64_t x) {
	if (x < 0) return -icbrtc(-x);
	int32_t q = int32_t(floor(pow(double(x), 1.0 / 3)));
	while (icb(q) > x) q--;
	while (icb(q + 1) <= x) q++;
	return q;
}

int32_t icbrtc(int64_t x) {
	if (x < 0) return -icbrt(-x);
	int32_t q = int32_t(ceil(pow(double(x), 1.0 / 3)));
	while (icb(q) < x) q++;
	while (icb(q - 1) >= x) q--;
	return q;
}

bool is_cube(int64_t x) {
	return (icb(icbrt(x)) == x);
}

} // math
} // altruct
