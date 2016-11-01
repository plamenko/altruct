#pragma once

#include "algorithm/math/base.h"

namespace altruct {
namespace math {

namespace quadratic_storage {
	enum type { INSTANCE, STATIC, CONSTANT };
}

template<typename T, int ID, int STORAGE_TYPE>
struct quadratic_members;

template<typename T, int ID>
struct quadratic_members<T, ID, quadratic_storage::INSTANCE> {
	T _D;
	quadratic_members(const T& _D = 0) : _D(_D) {}
	const T& D() const { return _D; }
	T& D() { return _D; }
};

template<typename T, int ID>
struct quadratic_members<T, ID, quadratic_storage::STATIC> {
	static T _D;
	quadratic_members(const T& _D = 0) {}
	static T& D() { return _D; }
};

template<typename T, int ID>
struct quadratic_members<T, ID, quadratic_storage::CONSTANT> {
	quadratic_members(const T& _D = 0) {}
	static T D() { return T(ID); }
};

/**
 * a + b * sqrt(D)
 *
 * quadratic<int, -1> - gaussian integers
 * quadratic<double, -1> - complex numbers
 *
 * @param T - the underlying type
 * @param ID - ID of the quadratic type (useful when STATIC = true)
 * @param STATIC - whether D is a static or instance member
 */
template<typename T, int ID, int STORAGE_TYPE = quadratic_storage::STATIC>
class quadratic : public quadratic_members<T, ID, STORAGE_TYPE> {
	typedef quadratic_members<T, ID, STORAGE_TYPE> my_quadratic_members;
public:
	T a, b;

	// construct from int, but only if T is not integral to avoid constructor clashing
	template <typename I = T, typename = std::enable_if_t<!std::is_integral<I>::value>>
	quadratic(int a) : my_quadratic_members(0), a(a), b(0) {}
	quadratic(const T& a = 0, const T& b = 0, const T& D = 0) : my_quadratic_members(D), a(a), b(b) {}
	quadratic(const quadratic& rhs) : my_quadratic_members(rhs.D()), a(rhs.a), b(rhs.b) {}

	bool operator == (const quadratic& rhs) const { return (a == rhs.a && b == rhs.b); }
	bool operator != (const quadratic& rhs) const { return (a != rhs.a || b != rhs.b); }
	bool operator <  (const quadratic& rhs) const { return (a < rhs.a || (a == rhs.a && b < rhs.b)); }
	bool operator >  (const quadratic& rhs) const { return (a > rhs.a || (a == rhs.a && b > rhs.b)); }
	bool operator <= (const quadratic& rhs) const { return (a < rhs.a || (a == rhs.a && b <= rhs.b)); }
	bool operator >= (const quadratic& rhs) const { return (a > rhs.a || (a == rhs.a && b >= rhs.b)); }

	quadratic  operator +  (const quadratic& rhs) const { quadratic t(*this); t += rhs; return t; }
	quadratic  operator -  (const quadratic& rhs) const { quadratic t(*this); t -= rhs; return t; }
	quadratic  operator -  ()                     const { quadratic t(-a, -b, this->D()); return t; }
	quadratic  operator *  (const quadratic& rhs) const { quadratic t(*this); t *= rhs; return t; }
	quadratic  operator /  (const quadratic& rhs) const { quadratic t(*this); t /= rhs; return t; }
	quadratic  operator %  (const quadratic& rhs) const { quadratic t(*this); t %= rhs; return t; }

	quadratic  operator *  (const T& rhs) const { quadratic t(*this); t *= rhs; return t; }
	quadratic  operator /  (const T& rhs) const { quadratic t(*this); t /= rhs; return t; }

	quadratic& operator += (const quadratic& rhs) { a += rhs.a; b += rhs.b; return *this; }
	quadratic& operator -= (const quadratic& rhs) { a -= rhs.a; b -= rhs.b; return *this; }
	quadratic& operator *= (const quadratic& rhs) { return *this = quadratic(a * rhs.a + b * rhs.b * this->D(), a * rhs.b + b * rhs.a, this->D()); }
	quadratic& operator /= (const quadratic& rhs) { return *this = *this * rhs.conjugate() / rhs.norm(); }
	quadratic& operator %= (const quadratic& rhs) { return *this = *this - rhs * (*this / rhs); }

	quadratic& operator *= (const T& rhs) { a *= rhs; b *= rhs; return *this; }
	quadratic& operator /= (const T& rhs) { a /= rhs; b /= rhs; return *this; }

	quadratic conjugate() const { return quadratic(a, -b, this->D()); }
	T norm() const { return a * a - b * b * this->D(); }
};

template<typename T>
using quadraticX = quadratic<T, 0, quadratic_storage::INSTANCE>;

template<typename T, int ID>
T quadratic_members<T, ID, quadratic_storage::STATIC>::_D = T(ID);

template<typename T, int ID, int STORAGE_TYPE>
struct identityT<quadratic<T, ID, STORAGE_TYPE>> {
	static quadratic<T, ID, STORAGE_TYPE> of(const quadratic<T, ID, STORAGE_TYPE>& x) {
		return quadratic<T, ID, STORAGE_TYPE>(identityT<T>::of(x.a), zeroT<T>::of(x.b), x.D());
	}
};

template<typename T, int ID, int STORAGE_TYPE>
struct zeroT<quadratic<T, ID, STORAGE_TYPE>> {
	static quadratic<T, ID, STORAGE_TYPE> of(const quadratic<T, ID, STORAGE_TYPE>& x) {
		return quadratic<T, ID, STORAGE_TYPE>(zeroT<T>::of(x.a), zeroT<T>::of(x.b), x.D());
	}
};

} // math
} // altruct
