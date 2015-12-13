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

// generic GCD
// note: for integral types and negative input the result might be of incorrect sign!
template<typename T>
T gcd(T a, T b) {
	while (a != 0) { T r = b % a; b = a; a = r; }
	return b;
}

// generic extended GCD
// note: for integral types and negative input the result might be of incorrect sign!
template<typename T>
T gcd_ex(const T& a, const T& b, T *x = 0, T *y = 0) {
	T r, q, g = a, h = b;
	T xo = 0, xn = 1;
//	T yo = 1, yn = 0;
	while (h != 0) {
		q = g / h;
		r = g  - q * h;  g  = h;  h  = r;
		r = xn - q * xo; xn = xo; xo = r;
//		r = yn - q * yo; yn = yo; yo = r;
//		T gn = a * xn + b * yn;
	}
	if (x) *x = xn;
//	if (y) *y = yn;
	if (y) *y = (b != 0) ? (g - a * xn) / b : 0;
	return g;
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
