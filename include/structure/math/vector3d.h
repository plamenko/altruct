#pragma once

#include "algorithm/math/base.h"

namespace altruct {
namespace math {

struct vector3d_empty_data {};
template <typename T, typename PAYLOAD = vector3d_empty_data>
class vector3d {
public:
	typedef T value_type;
	static T EPS;
	T x, y, z;
	PAYLOAD data;

	vector3d(const T& x = 0, const T& y = 0, const T& z = 0, const PAYLOAD& data = PAYLOAD()) : x(x), y(y), z(z), data(data) {}
	
	bool near        (const vector3d& v) const { return absT(x - v.x) <= EPS && absT(y - v.y) <= EPS && absT(z - v.z) <= EPS; }
	bool operator == (const vector3d& v) const { return (x == v.x) && (y == v.y) && (z == v.z); }
	bool operator != (const vector3d& v) const { return !(*this == v); }
	bool operator <  (const vector3d& v) const { return (x != v.x) ? (x < v.x) : ((y != v.y) ? (y < v.y) : (z < v.z)); }
	bool operator >  (const vector3d& v) const { return (v < *this); }
	bool operator <= (const vector3d& v) const { return !(v < *this); }
	bool operator >= (const vector3d& v) const { return !(*this < v); }

	vector3d& operator += (const vector3d& v) { x += v.x; y += v.y; z += v.z; return *this; }
	vector3d& operator -= (const vector3d& v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
	vector3d& operator *= (const vector3d& v) { x *= v.x; y *= v.y; z *= v.z; return *this; } // element by element
	vector3d& operator /= (const vector3d& v) { x /= v.x; y /= v.y; z /= v.z; return *this; } // element by element
	
	vector3d& operator *= (const T& t) { x *= t; y *= t; z *= t; return *this; }
	vector3d& operator /= (const T& t) { x /= t; y /= t; z /= t; return *this; }
	
	vector3d  operator +  (const vector3d& v) const { return vector3d(x + v.x, y + v.y, z + v.z); }
	vector3d  operator -  (const vector3d& v) const { return vector3d(x - v.x, y - v.y, z - v.z); }
	vector3d  operator -  ()                  const { return vector3d(-x, -y, -z); }
	vector3d  operator *  (const vector3d& v) const { return vector3d(x * v.x, y * v.y, z * v.z); } // element by element
	vector3d  operator /  (const vector3d& v) const { return vector3d(x / v.x, y / v.y, z / v.z); } // element by element

	vector3d  operator *  (const T& t) const { return vector3d(x * t, y * t, z * t); }
	vector3d  operator /  (const T& t) const { return vector3d(x / t, y / t, z / t); }

	vector3d& operator ^= (const vector3d& v) { return *this = *this ^ v; } // cross product
	
	T         operator &  (const vector3d& v) const { return (x * v.x + y * v.y + z * v.z); } // dot product
	vector3d  operator ^  (const vector3d& v) const { return vector3d(y * v.z - v.y * z, v.x * z - x * v.z, x * v.y - v.x * y); } // cross product
	
	T         dot         (const vector3d& v1, const vector3d& v2) const { return ((v1 - *this) & (v2 - *this)); } // dot product
	vector3d  cross       (const vector3d& v1, const vector3d& v2) const { return ((v1 - *this) ^ (v2 - *this)); } // cross product
	
	vector3d  unit        (const vector3d& v0 = vector3d()) const { T d = abs1(); return ((d > EPS) ? *this / d : v0); }
	
	T         abs1        () const { return (sqrtT(abs2())); }
	T         abs2        () const { return (x * x + y * y + z * z); }
};

// comparison tolerance
template <typename T, typename PAYLOAD> T vector3d<T, PAYLOAD>::EPS = 0;

} // math
} // altruct
