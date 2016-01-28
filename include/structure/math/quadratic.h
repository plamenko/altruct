#pragma once

#include "algorithm/math/base.h"

namespace altruct {
namespace math {

// a + b * sqrt(D)
template<typename T = int, int ID = -1>
class quadratic {
public:
	static T D;

	T a, b;

	quadratic(const quadratic &rhs) : a(rhs.a), b(rhs.b) {}
	quadratic(T a = 0, T b = 0) : a(a), b(b) {}

	bool operator == (const quadratic &rhs) const { return (a == rhs.a && b == rhs.b); }
	bool operator != (const quadratic &rhs) const { return (a != rhs.a || b != rhs.b); }
	bool operator <  (const quadratic &rhs) const { return (a < rhs.a || (a == rhs.a && b < rhs.b)); }
	bool operator >  (const quadratic &rhs) const { return (a > rhs.a || (a == rhs.a && b > rhs.b)); }
	bool operator <= (const quadratic &rhs) const { return (a < rhs.a || (a == rhs.a && b <= rhs.b)); }
	bool operator >= (const quadratic &rhs) const { return (a > rhs.a || (a == rhs.a && b >= rhs.b)); }
	
	quadratic  operator +  (const quadratic &rhs) const { quadratic t(*this); t += rhs; return t; }
	quadratic  operator -  (const quadratic &rhs) const { quadratic t(*this); t -= rhs; return t; }
	quadratic  operator -  ()                     const { quadratic t(-a, -b);          return t; }
	quadratic  operator *  (const quadratic &rhs) const { quadratic t(*this); t *= rhs; return t; }
	quadratic  operator /  (const quadratic &rhs) const { quadratic t(*this); t /= rhs; return t; }
	quadratic  operator %  (const quadratic &rhs) const { quadratic t(*this); t %= rhs; return t; }

	quadratic  operator *  (const T &rhs) const { quadratic t(*this); t *= rhs; return t; }
	quadratic  operator /  (const T &rhs) const { quadratic t(*this); t /= rhs; return t; }
	
	quadratic& operator += (const quadratic &rhs) { a += rhs.a; b += rhs.b; return *this; }
	quadratic& operator -= (const quadratic &rhs) { a -= rhs.a; b -= rhs.b; return *this; }
	quadratic& operator *= (const quadratic &rhs) { return *this = quadratic(a * rhs.a + b * rhs.b * D, a * rhs.b + b * rhs.a); }
	quadratic& operator /= (const quadratic &rhs) { return *this = *this * rhs.conjugate() / rhs.norm(); }
	quadratic& operator %= (const quadratic &rhs) { return *this = *this - rhs * (*this / rhs); }
	
	quadratic& operator *= (const T &rhs) { a *= rhs; b *= rhs; return *this; }
	quadratic& operator /= (const T &rhs) { a /= rhs; b /= rhs; return *this; }

	quadratic conjugate() const { return quadratic(a, -b); }
	T norm() const { return a * a - b * b * D; }
};

template<typename T, int ID>
T quadratic<T, ID>::D = (T)ID;

} // math
} // altruct
