#pragma once

#include "algorithm/math/base.h"

namespace altruct {
namespace math {

// modulo normalization
template<typename T>
void modulo_normalize(T& v, const T& M) { v %= M; }
// integral type specializations
template<typename I>
void modulo_normalize_int(I& v, I M) { if (v < 0 || v >= M) v %= M; if (v < 0) v += M; }
void modulo_normalize(long long& v, long long M);
void modulo_normalize(int& v, int M);

// modulo multiplication
template<typename T>
T modulo_multiply(const T& x, const T& y, const T& M) { return (x * y) % M; }
// integral type specializations
template<typename I>
I modulo_long_multiply_int(I x, I y, I M) {
	modulo_normalize_int(x, M);
	modulo_normalize_int(y, M);
	I r = 0;
	for (; y > 0; y >>= 1) {
		if (y & 1) r += x, modulo_normalize_int(r, M);
		x += x, modulo_normalize_int(x, M);
	}
	return r;
}
long long modulo_multiply(long long x, long long y, long long M);
int modulo_multiply(int x, int y, int M);

// modulo inversion
template<typename T> T modulo_inverse(const T& v, const T& M) {
	T vi; T g = gcd_ex(v, M, &vi);
	if (g != 1) vi /= g;
	return vi;
}
// integral type specializations
template<typename I> I modulo_inverse_int(I v, I M) {
	modulo_normalize_int(v, M);
	I vi; gcd_ex(v, M, &vi);
	modulo_normalize_int(vi, M);
	return vi;
}
long long modulo_inverse(long long v, long long M);
int modulo_inverse(int v, int M);

// modulo division
template<typename T> T modulo_divide(const T& x, const T& y, const T& M) {
	return modulo_multiply(x, modulo_inverse(y, M), M);
}
// integral type specializations
template<typename I>
I modulo_divide_int(I x, I y, I M) {
	modulo_normalize_int(x, M);
	modulo_normalize_int(y, M);
	return (x % y == 0) ? x / y : modulo_multiply(x, modulo_inverse(y, M), M);
}
long long modulo_divide(long long x, long long y, long long M);
int modulo_divide(int x, int y, int M);

// T must be constructable from int
template<typename T = long long, int ID = 1000000007>
class modulo {
public:
	static T M;

	T v;

	modulo(const modulo &rhs) : v(rhs.v) {}
	modulo(T v = 0) : v(v) { normalize(); }

	void normalize() { modulo_normalize(v, M); }
	
	bool operator == (const modulo &rhs) const { return (v == rhs.v); }
	bool operator != (const modulo &rhs) const { return (v != rhs.v); }
	bool operator <  (const modulo &rhs) const { return (v <  rhs.v); }
	bool operator >  (const modulo &rhs) const { return (v >  rhs.v); }
	bool operator <= (const modulo &rhs) const { return (v <= rhs.v); }
	bool operator >= (const modulo &rhs) const { return (v >= rhs.v); }
	
	modulo  operator +  (const modulo &rhs) const { modulo t(*this); t += rhs; return t; }
	modulo  operator -  (const modulo &rhs) const { modulo t(*this); t -= rhs; return t; }
	modulo  operator -  ()                  const { modulo t(-v);              return t; }
	modulo  operator *  (const modulo &rhs) const { modulo t(*this); t *= rhs; return t; }
	modulo  operator /  (const modulo &rhs) const { modulo t(*this); t /= rhs; return t; }
	modulo  operator %  (const modulo &rhs) const { modulo t(*this); t %= rhs; return t; }

	modulo& operator += (const modulo &rhs) { v += rhs.v;     normalize(); return *this; }
	modulo& operator -= (const modulo &rhs) { v -= rhs.v;     normalize(); return *this; }
	modulo& operator *= (const modulo &rhs) { v = modulo_multiply(v, rhs.v, M); return *this; }
	modulo& operator /= (const modulo &rhs) { v = modulo_divide(v, rhs.v, M); return *this; }
	modulo& operator %= (const modulo &rhs) { v %= rhs.v;     normalize(); return *this; }
};

template<typename T, int ID>
T modulo<T, ID>::M = (T) ID;

} // math
} // altruct
