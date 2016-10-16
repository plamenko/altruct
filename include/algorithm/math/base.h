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
template<typename T>
T identityOf(const T& x) { return identityT<T>::of(x); }

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
template<typename T>
T zeroOf(const T& x) { return zeroT<T>::of(x); }

/**
 * Absolute value.
 */
template <typename T>
T absT(const T &x) {
	T e0 = zeroT<T>::of(x);
	return (x < e0) ? -x : x;
}

/**
 * Minimum.
 */
template <typename T>
T minT(const T &x, const T &y) {
	return (x < y) ? x : y;
}

/**
 * Maximum.
 */
template <typename T>
T maxT(const T &x, const T &y) {
	return (x < y) ? y : x;
}

/**
 * Exponentiation by squaring.
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
		if (y % 2 != 0) r *= x;
		x *= x;
	}
	return r;
}

/**
 * Greatest Common Divisor.
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
 * Extended Greatest Common Divisor.
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
	T yo = e1, yn = e0;
	while (h != e0) {
		q = g / h;
		r = g  - q * h;  g  = h;  h  = r;
		r = xn - q * xo; xn = xo; xo = r;
		r = yn - q * yo; yn = yo; yo = r;
//		T gn = a * xn + b * yn;
	}
	if (x) *x = xn;
	if (y) *y = yn;
//	if (y) *y = (b != e0) ? (g - a * xn) / b : e0;
	return g;
}

/**
 * Maximal divisor `g` of `b`, such that `squarefree_kernel(g)` divides `a`.
 *
 * Formula: `g = gcd_max(a, b) = gcd(a ^ inf, b)`
 * The following holds: `gcd(a, b / g) = 1`
 * In comparison, regular `gcd` function returns the
 * maximal divisor `g` of `b`, such that `g` divides `a`.
 */
template<typename T>
T gcd_max(T a, T b) {
	T e0 = zeroT<T>::of(b);
	if (b == e0) return a;
	T go = e0, g = 1;
	while (go != g) {
		go = g; g = gcd(g * a, b);
	}
	return g;
}

/**
 * Least Common Multiple.
 *
 * Note: for integral types and negative input the result might be of incorrect sign!
 */
template<typename T>
T lcm(const T& a, const T& b) {
	return a * (b / gcd(a, b));
}

/**
 * Integer square & square root.
 */
int64_t isq(int64_t x);
int32_t isqrt(int64_t x);
int32_t isqrtc(int64_t x);

/**
 * Integer cube & cube root.
 */
int64_t icb(int64_t x);
int32_t icbrt(int64_t x);
int32_t icbrtc(int64_t x);

/**
 * Square of `x`.
 */
template<typename T>
T sqT(T x) {
	return x * x;
}

/**
 * Square root of `x`, rounded towards 0.
 *
 * Note: for negative `x` result is `-sqrtT(-x)`.
 * Base implementation uses the Newton-Raphson method in O(log3(n)).
 */
template<typename T>
T sqrtT(T x, T eps = T(1)) {
	if (x < 0) return -sqrtT<T>(-x, eps);
	if (x == 0) return 0;
	if (x == 1) return 1;
	T q1 = x / 2;
	T q2 = x / q1;
	while (absT(T(q1 - q2)) > eps) {
		q1 = (q1 + q2) / 2;
		q2 = x / q1;
	}
	return minT(q1, q2);
}
template<> inline float sqrtT(float x, float) { return sqrt(x); }
template<> inline double sqrtT(double x, double) { return sqrt(x); }
template<> inline long double sqrtT(long double x, long double) { return sqrt(x); }
template<> inline int8_t sqrtT(int8_t x, int8_t) { return isqrt(x); }
template<> inline uint8_t sqrtT(uint8_t x, uint8_t) { return isqrt(x); }
template<> inline int16_t sqrtT(int16_t x, int16_t) { return isqrt(x); }
template<> inline uint16_t sqrtT(uint16_t x, uint16_t) { return isqrt(x); }
template<> inline int32_t sqrtT(int32_t x, int32_t) { return isqrt(x); }
template<> inline uint32_t sqrtT(uint32_t x, uint32_t) { return isqrt(x); }
template<> inline int64_t sqrtT(int64_t x, int64_t) { return isqrt(x); }
template<> inline uint64_t sqrtT(uint64_t x, uint64_t) { return isqrt(x); }

/**
 * Tests whether `x` is a square.
 */
template<typename I>
bool is_square(I x) {
	return (sqT(sqrtT(x)) == x);
}

/**
 * Cube of `x`.
 */
template<typename T>
T cbT(T x) {
	return x * x * x;
}

/**
 * Cube root of `x`, rounded towards 0.
 *
 * Note: for negative `x` result is `-cbrtT(-x)`.
 * Base implementation uses the Newton-Raphson method in O(log3(n)).
 */
template<typename T>
T cbrtT(T x, T eps = T(1)) {
	if (x < 0) return -cbrtT<T>(-x, eps);
	if (x == 0) return 0;
	if (x == 1) return 1;
	T r0 = 0;
	T r1 = sqrtT(x, eps);
	T r2 = x / sqT(r1);
	while (r1 != r0 && absT(T(r1 - r2)) > eps) {
		r0 = r1;
		r1 = (r1 + r1 + r2) / 3;
		r2 = x / sqT(r1);
	}
	return minT(r1, r2);
}
template<> inline float cbrtT(float x, float) { return cbrt(x); }
template<> inline double cbrtT(double x, double) { return cbrt(x); }
template<> inline int16_t cbrtT(int16_t x, int16_t) { return icbrt(x); }
template<> inline uint16_t cbrtT(uint16_t x, uint16_t) { return icbrt(x); }
template<> inline int32_t cbrtT(int32_t x, int32_t) { return icbrt(x); }
template<> inline uint32_t cbrtT(uint32_t x, uint32_t) { return icbrt(x); }
template<> inline int64_t cbrtT(int64_t x, int64_t) { return icbrt(x); }
template<> inline uint64_t cbrtT(uint64_t x, uint64_t) { return icbrt(x); }

/**
 * Tests whether `x` is a cube.
 */
template<typename I>
bool is_cube(I x) {
	return (cbT(cbrtT(x)) == x);
}

/**
 * Integer floor & ceil division.
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

template<typename I>
I div_round(I a, I b) {
	if (b < 0) a = -a, b = -b;
	return (a > 0) ? (a + b / 2) / b : (a - b / 2) / b;
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
