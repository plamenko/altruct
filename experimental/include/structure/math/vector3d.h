#pragma once

#include "algorithm/math/base.h"

namespace altruct {
namespace math {

template <typename T>
class vector3d {
public:
	static T EPS;
	T x, y, z;
	
	vector3d(const T& x = 0, const T& y = 0, const T& z = 0) : x(x), y(y), z(z) {}
	
	bool near(const vector3d<T>& v) const { return absT(x - v.x) <= EPS && absT(y - v.y) <= EPS && absT(z - v.z) <= EPS; }
	bool operator == (const vector3d<T>& v) const { return (x == v.x) && (y == v.y) && (z == v.z); }
	bool operator != (const vector3d<T>& v) const { return !(*this == v); }
	bool operator <  (const vector3d<T>& v) const { return (x != v.x) ? (x < v.x) : ((y != v.y) ? (y < v.y) : (z < v.z)); }
	bool operator >  (const vector3d<T>& v) const { return (v < *this); }
	bool operator <= (const vector3d<T>& v) const { return !(v < *this); }
	bool operator >= (const vector3d<T>& v) const { return !(*this < v); }

	vector3d<T>& operator += (const vector3d<T>& v) { x += v.x; y += v.y; z += v.z; return *this; }
	vector3d<T>& operator -= (const vector3d<T>& v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
	vector3d<T>& operator *= (const vector3d<T>& v) { x *= v.x; y *= v.y; z *= v.z; return *this; } // element by element
	vector3d<T>& operator /= (const vector3d<T>& v) { x /= v.x; y /= v.y; z /= v.z; return *this; } // element by element
	
	vector3d<T>& operator *= (const T& t) { x *= t; y *= t; z *= t; return *this; }
	vector3d<T>& operator /= (const T& t) { x /= t; y /= t; z /= t; return *this; }
	
	vector3d<T>  operator +  (const vector3d<T>& v) const { return vector3d(x + v.x, y + v.y, z + v.z); }
	vector3d<T>  operator -  (const vector3d<T>& v) const { return vector3d(x - v.x, y - v.y, z - v.z); }
	vector3d<T>  operator -  ()                     const { return vector3d(     -x,      -y,      -z); }
	vector3d<T>  operator *  (const vector3d<T>& v) const { return vector3d(x * v.x, y * v.y, z * v.z); } // element by element
	vector3d<T>  operator /  (const vector3d<T>& v) const { return vector3d(x / v.x, y / v.y, z / v.z); } // element by element

	vector3d<T>  operator *  (const T& t) const { return vector3d(x * t, y * t, z * t); }
	vector3d<T>  operator /  (const T& t) const { return vector3d(x / t, y / t, z / t); }

	vector3d<T>& operator ^= (const vector3d<T>& v) { return *this = *this ^ v; } // cross product
	
	T            operator &  (const vector3d<T>& v) const { return (x * v.x + y * v.y + z * v.z); } // dot product
	vector3d<T>  operator ^  (const vector3d<T>& v) const { return vector3d(y * v.z - v.y * z, v.x * z - x * v.z, x * v.y - v.x * y); } // cross product
	
	T            dot         (const vector3d<T>& v1, const vector3d<T>& v2) const { return ((v1 - *this) & (v2 - *this)); } // dot product
	vector3d<T>  cross       (const vector3d<T>& v1, const vector3d<T>& v2) const { return ((v1 - *this) ^ (v2 - *this)); } // cross product
	
	vector3d<T>  unit        (const vector3d<T>& v0 = vector3d()) const { T d = abs1(); return ((d > EPS) ? *this / d : v0); }
	
	T            abs1        () const { return (sqrtT(abs2())); }
	T            abs2        () const { return (x * x + y * y + z * z); }
};

// comparison tolerance
template <typename T> T vector3d<T>::EPS = 0;

} // math
} // altruct
