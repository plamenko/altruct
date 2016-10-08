#pragma once

#include "algorithm/math/base.h"

namespace altruct {
namespace math {

struct vector2d_empty_data {};
template <typename T, typename PAYLOAD = vector2d_empty_data>
class vector2d {
public:
	typedef T value_type;
	static T EPS;
	T x, y;
	PAYLOAD data;

	vector2d(const T& x = 0, const T& y = 0, const PAYLOAD& data = PAYLOAD()) : x(x), y(y), data(data) {}
	
	bool near        (const vector2d& v) const { return absT(x - v.x) <= EPS && absT(y - v.y) <= EPS; }
	bool operator == (const vector2d& v) const { return (x == v.x) && (y == v.y); }
	bool operator != (const vector2d& v) const { return !(*this == v); }
	bool operator <  (const vector2d& v) const { return (x != v.x) ? (x < v.x) : (y < v.y); }
	bool operator >  (const vector2d& v) const { return (v < *this); }
	bool operator <= (const vector2d& v) const { return !(v < *this); }
	bool operator >= (const vector2d& v) const { return !(*this < v); }

	vector2d& operator += (const vector2d& v) { x += v.x; y += v.y; return *this; }
	vector2d& operator -= (const vector2d& v) { x -= v.x; y -= v.y; return *this; }
	vector2d& operator *= (const vector2d& v) { x *= v.x; y *= v.y; return *this; } // (element by element)
	vector2d& operator /= (const vector2d& v) { x /= v.x; y /= v.y; return *this; } // (element by element)
	
	vector2d& operator *= (const T& t) { x *= t; y *= t; return (*this); }
	vector2d& operator /= (const T& t) { x /= t; y /= t; return (*this); }
	
	vector2d  operator +  (const vector2d& v) const { return vector2d(x + v.x, y + v.y); }
	vector2d  operator -  (const vector2d& v) const { return vector2d(x - v.x, y - v.y); }
	vector2d  operator -  ()                  const { return vector2d(-x, -y); }
	vector2d  operator *  (const vector2d& v) const { return vector2d(x * v.x, y * v.y); } // (element by element)
	vector2d  operator /  (const vector2d& v) const { return vector2d(x / v.x, y / v.y); } // (element by element)
	
	vector2d  operator *  (const T& t) const { return vector2d(x * t, y * t); }
	vector2d  operator /  (const T& t) const { return vector2d(x / t, y / t); }
	
	T         operator &  (const vector2d& v) const { return (x * v.x + y * v.y); } // dot product
	T         operator ^  (const vector2d& v) const { return (x * v.y - y * v.x); } // cross product
	
	T         dot         (const vector2d& v1, const vector2d& v2) const { return ((v1 - *this) & (v2 - *this)); } // dot product
	T         cross       (const vector2d& v1, const vector2d& v2) const { return ((v1 - *this) ^ (v2 - *this)); } // cross product
	
	vector2d  unit        (const vector2d& v0 = vector2d()) const { T d = abs1(); return (d > EPS) ? *this / d : v0; }
	vector2d  rot         (const vector2d& r)               const { return vector2d(x * r.x - y * r.y, y * r.x + x * r.y); }
	vector2d  irot        (const vector2d& r)               const { return vector2d(x * r.x + y * r.y, y * r.x - x * r.y); }
	vector2d  perp        ()                                const { return vector2d(-y, x); }
	
	T         abs1        () const { return (sqrtT(abs2())); }
	T         abs2        () const { return (x * x + y * y); }
	T         diff2       () const { return (x * x - y * y); }
};

// comparison tolerance
template <typename T, typename PAYLOAD> T vector2d<T, PAYLOAD>::EPS = 0;

} // math
} // altruct
