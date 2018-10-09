#pragma once

#include "altruct/algorithm/math/base.h"

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
 * @param ID - ID of the quadratic type (useful with quadratic_storage::CONSTANT)
 * @param STORAGE_TYPE - whether D is a constant, static or instance member
 *   See `STORAGE_TYPE` in the `modulo` class for details.
 */
template<typename T, int ID, int STORAGE_TYPE = quadratic_storage::STATIC>
class quadratic : public quadratic_members<T, ID, STORAGE_TYPE> {
    typedef quadratic_members<T, ID, STORAGE_TYPE> my_quadratic_members;
public:
    T a, b;

    // construct from int, but only if T is not integral to avoid constructor clashing
    template <typename I = T, typename = std::enable_if_t<!std::is_integral<I>::value>>
    quadratic(int a) : my_quadratic_members(0), a(castOf(this->D(), a)), b(zeroOf(this->D())) {}
    quadratic() : my_quadratic_members(-1), a(zeroOf(this->D())), b(zeroOf(this->D())) {}
    quadratic(const T& a) : my_quadratic_members(-1), a(a), b(zeroOf(a)) {}
    quadratic(const T& a, const T& b, const T& D = -1) : my_quadratic_members(D), a(a), b(b) {}
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
T quadratic_members<T, ID, quadratic_storage::STATIC>::_D = castOf<T>(ID);

template<typename T, int ID, int STORAGE_TYPE, typename I>
struct castT<quadratic<T, ID, STORAGE_TYPE>, I> {
    typedef quadratic<T, ID, STORAGE_TYPE> quad;
    static quad of(const I& x) {
        return quad(castOf<T>(x));
    }
    static quad of(const quad& ref, const I& x) {
        return quad(castOf(ref.a, x), castOf(ref.b, zeroOf(x)), ref.D());
    }
};
template<typename T, int ID, int STORAGE_TYPE>
struct castT<quadratic<T, ID, STORAGE_TYPE>, quadratic<T, ID, STORAGE_TYPE>> : nopCastT<quadratic<T, ID, STORAGE_TYPE>>{};

template<typename T, int ID, int STORAGE_TYPE>
struct identityT<quadratic<T, ID, STORAGE_TYPE>> {
    typedef quadratic<T, ID, STORAGE_TYPE> quad;
    static quad of(const quad& x) {
        return quad(identityOf(x.a), zeroOf(x.b), x.D());
    }
};

template<typename T, int ID, int STORAGE_TYPE>
struct zeroT<quadratic<T, ID, STORAGE_TYPE>> {
    typedef quadratic<T, ID, STORAGE_TYPE> quad;
    static quad of(const quad& x) {
        return quad(zeroOf(x.a), zeroOf(x.b), x.D());
    }
};

template<typename T, int ID, int STORAGE_TYPE>
struct conjugateT<quadratic<T, ID, STORAGE_TYPE>> {
    typedef quadratic<T, ID, STORAGE_TYPE> quad;
    static quad of(const quad& x) {
        return x.conjugate();
    }
};

} // math
} // altruct
