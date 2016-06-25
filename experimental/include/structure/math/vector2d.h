#pragma once

#include "algorithm/math/base.h"

namespace altruct {
namespace math {

template <typename T>
class vector2d {
public:
	static T EPS;
	T x, y;
	
	vector2d(const T& x = 0, const T& y = 0) : x(x), y(y) {}
	
	bool near(const vector2d<T>& v) const { return absT(x - v.x) <= EPS && absT(y - v.y) <= EPS; }
	bool operator == (const vector2d<T>& v) const { return (x == v.x) && (y == v.y); }
	bool operator != (const vector2d<T>& v) const { return !(*this == v); }
	bool operator <  (const vector2d<T>& v) const { return (x != v.x) ? (x < v.x) : (y < v.y); }
	bool operator >  (const vector2d<T>& v) const { return (v < *this); }
	bool operator <= (const vector2d<T>& v) const { return !(v < *this); }
	bool operator >= (const vector2d<T>& v) const { return !(*this < v); }

	vector2d<T>& operator += (const vector2d<T>& v) { x += v.x; y += v.y; return *this; }
	vector2d<T>& operator -= (const vector2d<T>& v) { x -= v.x; y -= v.y; return *this; }
	vector2d<T>& operator *= (const vector2d<T>& v) { x *= v.x; y *= v.y; return *this; } // (element by element)
	vector2d<T>& operator /= (const vector2d<T>& v) { x /= v.x; y /= v.y; return *this; } // (element by element)
	
	vector2d<T>& operator *= (const T& t) { x *= t; y *= t; return (*this); }
	vector2d<T>& operator /= (const T& t) { x /= t; y /= t; return (*this); }
	
	vector2d<T>  operator +  (const vector2d<T>& v) const { return vector2d(x + v.x, y + v.y); }
	vector2d<T>  operator -  (const vector2d<T>& v) const { return vector2d(x - v.x, y - v.y); }
	vector2d<T>  operator -  ()                     const { return vector2d(     -x,      -y); }
	vector2d<T>  operator *  (const vector2d<T>& v) const { return vector2d(x * v.x, y * v.y); } // (element by element)
	vector2d<T>  operator /  (const vector2d<T>& v) const { return vector2d(x / v.x, y / v.y); } // (element by element)
	
	vector2d<T>  operator *  (const T& t) const { return vector2d(x * t, y * t); }
	vector2d<T>  operator /  (const T& t) const { return vector2d(x / t, y / t); }
	
	T            operator &  (const vector2d<T>& v) const { return (x * v.x + y * v.y); } // dot product
	T            operator ^  (const vector2d<T>& v) const { return (x * v.y - y * v.x); } // cross product
	
	T            dot         (const vector2d<T>& v1, const vector2d<T>& v2) const { return ((v1 - *this) ^ (v2 - *this)); } // dot product
	T            cross       (const vector2d<T>& v1, const vector2d<T>& v2) const { return ((v1 - *this) ^ (v2 - *this)); } // cross product
	
	vector2d<T>  unit        (const vector2d<T>& v0 = vector2d()) const { T d = abs1(); return (d > EPS) ? *this / d : v0; }
	vector2d<T>  rot         (const vector2d<T>& r)               const { return vector2d(x * r.x - y * r.y, y * r.x + x * r.y); }
	vector2d<T>  perp        ()                                   const { return vector2d(-y, x); }
	
	T            abs1        () const { return (sqrtT(abs2())); }
	T            abs2        () const { return (x * x + y * y); }
	T            diff2       () const { return (x * x - y * y); }
};

// comparison tolerance
template <typename T> T vector2d<T>::EPS = 0;

} // math
} // altruct
