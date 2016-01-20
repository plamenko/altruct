#pragma once

#include <cmath>
#include <stdint.h>

namespace altruct {
namespace math {

/**
 * Exponentiation by squaring
 *
 * @param x - base
 * @param y - exponent
 * @param id - multiplicative identity
 * @return x^y
 */
template<typename T>
T powT(T x, int64_t y, T id = T(1)) {
	T r = id;
	for (; y > 0; y >>= 1) {
		if (y & 1) r *= x;
		x *= x;
	}
	return r;
}

/**
 * Greatest Common Divisor
 *
 * Note: for integral types and negative input the result might be of incorrect sign!
 */
template<typename T>
T gcd(T a, T b) {
	while (a != 0) { T r = b % a; b = a; a = r; }
	return b;
}

/**
 * Extended Greatest Common Divisor
 *
 * Calculates `x`, `y` and `g` so that: `a * x + b * y = g`.
 * Note: for integral types and negative input the result might be of incorrect sign!
 *
 * @param a - first operand
 * @param b - second operand
 * @param x - output argument `x`
 * @param y - output argument `y`
 * @return g
 */
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

/**
 * Least Common Multiple
 *
 * Note: for integral types and negative input the result might be of incorrect sign!
 */
template<typename T>
T lcm(const T& a, const T& b) {
	return a * (b / gcd(a, b));
}

/**
 * Chinese Remainder
 *
 * Calculates `a` and `n` so that:
 *   n = lcm(n1, n2)
 *   a % n1 == a1
 *   a % n2 == a2
 *   0 <= a < n
 * `n1` and `n2` don't have to be coprime.
 * Correctly handles 64bit result for long long type.
 */
template<typename T>
void chinese_remainder(T &a, T&n, T a1, T n1, T a2, T n2) {
	T ni1, ni2; T g = gcd_ex(n1, n2, &ni1, &ni2);
	if ((a2 - a1) % g != 0) { n = 0, a = 0; return; }
	T t1 = (a1 * ni2) % n1;
	T t2 = (a2 * ni1) % n2;
	n1 /= g; n2 /= g; n = n1 * n2 * g;
	a = (t1 * n2 + t2 * n1) % n; if (a < 0) a += n;
}

/**
 * Square
 */
template<typename T> T sqT(T x) {
	return x * x;
}

/**
 * Integer square & square root
 */
int64_t isq(int64_t x);
int32_t isqrt(int64_t x);
int32_t isqrtc(int64_t x);
bool is_square(int64_t x);

/**
 * Integer cube & cube root
 */
int64_t icb(int64_t x);
int32_t icbrt(int64_t x);
int32_t icbrtc(int64_t x);
bool is_cube(int64_t x);

/**
 * Integer floor & ceil division
 */
int64_t div_floor(int64_t a, int64_t b);
int64_t div_ceil(int64_t a, int64_t b);

} // math
} // altruct
