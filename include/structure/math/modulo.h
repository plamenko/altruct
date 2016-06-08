#pragma once

#include "algorithm/math/base.h"

#include <type_traits>

namespace altruct {
namespace math {

// modulo normalization
template<typename T>
void modulo_normalize(T* v, const T& M) { *v %= M; }
// integral type specializations
template<typename I>
void modulo_normalize_int(I* v, const I& M) {
	if (*v < 0) *v += M;
	if (*v >= M) *v -= M;
	if (*v < 0 || *v >= M) *v %= M;
	if (*v < 0) *v += M;
}
inline void modulo_normalize(long long* v, long long M) { modulo_normalize_int(v, M); }
inline void modulo_normalize(int* v, int M) { modulo_normalize_int(v, M); }

// modulo multiplication
template<typename T>
T modulo_multiply(const T& x, const T& y, const T& M) { return (x * y) % M; }
// integral type specializations
template<typename I>
I modulo_long_multiply_int(I x, I y, I M) {
	modulo_normalize_int(&x, M);
	modulo_normalize_int(&y, M);
	I r = 0;
	for (; y > 0; y >>= 1) {
		if (y & 1) r += x, modulo_normalize_int(&r, M);
		x += x, modulo_normalize_int(&x, M);
	}
	return r;
}
inline long long modulo_multiply(long long x, long long y, long long M) {
	modulo_normalize_int(&x, M);
	modulo_normalize_int(&y, M);
	bool fit = (x < (1LL << 31) && y < (1LL << 31));
	return fit ? (x * y) % M : modulo_long_multiply_int(x, y, M);
}
inline int modulo_multiply(int x, int y, int M) {
	modulo_normalize_int(&x, M);
	modulo_normalize_int(&y, M);
	return ((long long)x * y) % M;
}

// modulo inversion
template<typename T> T modulo_inverse(const T& v, const T& M) {
	T vi; T g = gcd_ex(v, M, &vi);
	if (g != 1) vi /= g;
	return vi;
}
// integral type specializations
template<typename I> I modulo_inverse_int(I v, I M) {
	modulo_normalize_int(&v, M);
	I vi; gcd_ex(v, M, &vi);
	modulo_normalize_int(&vi, M);
	return vi;
}
inline long long modulo_inverse(long long v, long long M) { return modulo_inverse_int(v, M); }
inline int modulo_inverse(int v, int M) { return modulo_inverse_int(v, M); }

// modulo division
template<typename T> T modulo_divide(const T& x, const T& y, const T& M) {
	return modulo_multiply(x, modulo_inverse(y, M), M);
}
// integral type specializations
template<typename I>
I modulo_divide_int(I x, I y, I M) {
	modulo_normalize_int(&x, M);
	modulo_normalize_int(&y, M);
	if (y != 0 && x % y == 0) return x / y; // fast path if y divides x
	I yi; I g = gcd_ex(y, M, &yi);
	if (g != 1) { // uh oh, y and M are not coprime, try common gcd
		g = gcd(x, g);
		x /= g; y /= g;
		gcd_ex(y, M / g, &yi); // k
		// if (k != 1), there is no result, or more precisely,
		// the result will be `k` times bigger than it should.
		modulo_normalize_int(&yi, M / g);
	}
	return modulo_multiply(x, yi, M);
}
inline long long modulo_divide(long long x, long long y, long long M) { return modulo_divide_int(x, y, M); }
inline int modulo_divide(int x, int y, int M) { return modulo_divide_int(x, y, M); }

template<typename T, int ID, bool STATIC = true>
struct modulo_members;

template<typename T, int ID>
struct modulo_members<T, ID, false> {
	T M;
	modulo_members(const T& M = 0) : M(M) {}
};

template<typename T, int ID>
struct modulo_members<T, ID, true> {
	static T M;
	modulo_members(const T& M = 0) {}
};

/**
 * Modulo M arithmetics
 *
 * modulo<int, 3> - Z/3Z
 *
 * @param T - the underlying type
 * @param ID - ID of the modulo type (useful when STATIC = true)
 * @param STATIC - whether M is a static or instance member
 */
template<typename T, int ID, bool STATIC = true>
class modulo : public modulo_members<T, ID, STATIC> {
	typedef modulo_members<T, ID, STATIC> my_modulo_members;
public:
	T v;

	// construct from int, but only if T is not integral to avoid constructor clashing
	template <typename I = T, typename = std::enable_if_t<!std::is_integral<I>::value>>
	modulo(int v) : v(v), my_modulo_members(1) { if (STATIC) normalize(); }
	modulo(const T& v = 0) : v(v), my_modulo_members(1) { if (STATIC) normalize(); }
	modulo(const T& v, const T& M) : v(v), my_modulo_members(M) { normalize(); }
	modulo(const modulo& rhs) : v(rhs.v), my_modulo_members(rhs.M) {}

	void normalize() { modulo_normalize(&v, this->M); }
	
	bool operator == (const modulo &rhs) const { return (v == rhs.v); }
	bool operator != (const modulo &rhs) const { return (v != rhs.v); }
	bool operator <  (const modulo &rhs) const { return (v <  rhs.v); }
	bool operator >  (const modulo &rhs) const { return (v >  rhs.v); }
	bool operator <= (const modulo &rhs) const { return (v <= rhs.v); }
	bool operator >= (const modulo &rhs) const { return (v >= rhs.v); }
	
	modulo  operator +  (const modulo &rhs) const { modulo t(*this); t += rhs; return t; }
	modulo  operator -  (const modulo &rhs) const { modulo t(*this); t -= rhs; return t; }
	modulo  operator -  ()                  const { modulo t(-v, this->M);     return t; }
	modulo  operator *  (const modulo &rhs) const { modulo t(*this); t *= rhs; return t; }
	modulo  operator /  (const modulo &rhs) const { modulo t(*this); t /= rhs; return t; }
	modulo  operator %  (const modulo &rhs) const { modulo t(*this); t %= rhs; return t; }

	modulo& operator += (const modulo &rhs) { v += rhs.v;     normalize(); return *this; }
	modulo& operator -= (const modulo &rhs) { v -= rhs.v;     normalize(); return *this; }
	modulo& operator *= (const modulo &rhs) { v = modulo_multiply(v, rhs.v, this->M); return *this; }
	modulo& operator /= (const modulo &rhs) { v = modulo_divide(v, rhs.v, this->M); return *this; }
	modulo& operator %= (const modulo &rhs) { v %= rhs.v;     normalize(); return *this; }
};

template<typename T>
using moduloX = modulo<T, 0, false>;

template<typename T, int ID>
T modulo_members<T, ID, true>::M = (T)ID;

template<typename T, int ID, bool STATIC>
struct identityT<modulo<T, ID, STATIC>> {
	static modulo<T, ID, STATIC> of(const modulo<T, ID, STATIC>& x) {
		return modulo<T, ID, STATIC>(identityT<T>::of(x.v), x.M);
	}
};

template<typename T, int ID, bool STATIC>
struct zeroT<modulo<T, ID, STATIC>> {
	static modulo<T, ID, STATIC> of(const modulo<T, ID, STATIC>& x) {
		return modulo<T, ID, STATIC>(zeroT<T>::of(x.v), x.M);
	}
};

template<typename T>
T modT(T v, const T& m) { modulo_normalize(&v, m); return v; }

template<typename T, typename I>
T modulo_power(const T& x, const I& y, const T& M) {
	return powT(moduloX<T>(x, M), y).v;
}

} // math
} // altruct
