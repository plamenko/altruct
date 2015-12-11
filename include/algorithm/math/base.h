#pragma once

#include <cmath>
#include <stdint.h>

namespace altruct {
namespace math {

// generic exponentiation by squaring
template<typename T>
T powT(T x, int64_t y) {
	T r = 1;
	for (; y > 0; y >>= 1) {
		if (y & 1) r *= x;
		x *= x;
	}
	return r;
}

// generic squaring
template<typename T> T sqT(T x) {
	return x * x;
}

// integer square & square root
int64_t isq(int64_t x);
int32_t isqrt(int64_t x);
int32_t isqrtc(int64_t x);
bool is_square(int64_t x);

// integer cube & cube root
int64_t icb(int64_t x);
int32_t icbrt(int64_t x);
int32_t icbrtc(int64_t x);
bool is_cube(int64_t x);

// integer floor & ceil division
int64_t div_floor(int64_t a, int64_t b);
int64_t div_ceil(int64_t a, int64_t b);

} // math
} // altruct
