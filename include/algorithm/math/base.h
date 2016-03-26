#pragma once

#include <cmath>
#include <stdint.h>

namespace altruct {
namespace math {

/**
 * Gives the multiplicative identity element for the element `x`.
 *
 * For example:
 * if `x` is a `5x5` matrix, `e` is an identity matrix of rank `5`.
 * If `x` is an integer modulo M, `e` is `1 (mod M)`.
 * If `x` is an integer, `e` is simply 1.
 *
 * Note: implemented as a class because C++ doesn't allow
 * partial template specializations for functions.
 *
 * @param x - the element to provide the multiplicative identity for
 * @return e - such that `e * x = x * e = x`
 */
template<typename T>
struct identityT {
	static T of(const T& x) {
		return T(1);
	}
};

/**
 * Gives the additive identity element (multiplicative zero) for the element `x`.
 *
 * For example:
 * if `x` is a `5x5` matrix, `e` is a `5x5` zero matrix.
 * If `x` is an integer modulo M, `e` is `0 (mod M)`.
 * If `x` is an integer, `e` is simply 0.
 *
 * Note: implemented as a class because C++ doesn't allow
 * partial template specializations for functions.
 *
 * @param x - the element to provide the additive identity for
 * @return e - such that `e + x = x + e = x`
 */
template<typename T>
struct zeroT {
	static T of(const T& x) {
		return T(0);
	}
};

/**
 * Absolute value
 */
template <typename T>
T absT(const T &x) {
	T e0 = zeroT<T>::of(x);
	return (x < e0) ? -x : x;
}

/**
 * Exponentiation by squaring
 *
 * @param x - base
 * @param y - exponent
 * @return x^y
 */
template<typename T, typename I>
T powT(T x, I y) {
	T e1 = identityT<T>::of(x);
	T r = e1;
	for (; y > 0; y /= 2) {
		if (y % 2) r *= x;
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
	T e0 = zeroT<T>::of(a);
	while (a != e0) { T r = b % a; b = a; a = r; }
	return b;
}

/**
 * Extended Greatest Common Divisor
 *
 * Calculates `x`, `y` and `g` so that: `a * x + b * y = g`.
 * Note: for integral types and negative input the result might be of incorrect sign!
 *
 * @param a - the first operand
 * @param b - the second operand
 * @param x - the output argument `x`
 * @param y - the output argument `y`
 * @return g - the greatest common divisor of `a` and `b`
 */
template<typename T>
T gcd_ex(const T& a, const T& b, T *x = 0, T *y = 0) {
	T e0 = zeroT<T>::of(a), e1 = identityT<T>::of(a);
	T r, q, g = a, h = b;
	T xo = e0, xn = e1;
//	T yo = e1, yn = e0;
	while (h != e0) {
		q = g / h;
		r = g  - q * h;  g  = h;  h  = r;
		r = xn - q * xo; xn = xo; xo = r;
//		r = yn - q * yo; yn = yo; yo = r;
//		T gn = a * xn + b * yn;
	}
	if (x) *x = xn;
//	if (y) *y = yn;
	if (y) *y = (b != e0) ? (g - a * xn) / b : e0;
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
template<typename I>
I div_floor(I a, I b) {
	if (b < 0) a = -a, b = -b;
	return (a < 0) ? (a + 1) / b - 1 : a / b;
}

template<typename I>
I div_ceil(I a, I b) {
	if (b < 0) a = -a, b = -b;
	return (a > 0) ? (a - 1) / b + 1 : a / b;
}

/**
 * Multiple of `a`, greater than or equal to `b`.
 */
template<typename I>
I multiple(I a, I b) {
	return div_ceil(b, a) * a;
}

} // math
} // altruct
