#pragma once

#include "altruct/algorithm/math/base.h"

#include <type_traits>

namespace altruct {
namespace math {

// modulo normalization
template<typename T>
void modulo_normalize(T* v, const T& M) { *v %= M; }
// integral type specializations
template<typename I>
void modulo_normalize_int(I* v, I M) { *v %= M; if (*v < 0) *v += M; }
inline void modulo_normalize(int64_t* v, int64_t M) { modulo_normalize_int(v, M); }
inline void modulo_normalize(int32_t* v, int32_t M) { modulo_normalize_int(v, M); }

// modulo addition
template<typename T>
T modulo_add(const T& x, const T& y, const T& M) { return (x + y) % M; }

// modulo subtraction
template<typename T>
T modulo_sub(const T& x, const T& y, const T& M) { return (x - y) % M; }
// integral type specializations, input is assumed to be normalized
template<typename I>
I modulo_sub_int(I x, I y, I M) { return (x + (M - y)) % M; }
inline int64_t modulo_sub(int64_t x, int64_t y, int64_t M) { return modulo_sub_int(x, y, M); }
inline int32_t modulo_sub(int32_t x, int32_t y, int32_t M) { return modulo_sub_int(x, y, M); }

// modulo multiplication
template<typename T>
T modulo_mul(const T& x, const T& y, const T& M) { return (x * y) % M; }
// integral type specializations, input is assumed to be normalized
template<typename I>
I modulo_mul_int_long(I x, I y, I M) {
    I r = 0;
    for (; y > 0; y >>= 1) {
        if (y & 1) r += x, r %= M;
        x += x, x %= M;
    }
    return r;
}
inline int64_t modulo_mul(int64_t x, int64_t y, int64_t M) {
    return ((x | y) >> 31 == 0) ? (x * y) % M : modulo_mul_int_long(x, y, M);
}
inline int32_t modulo_mul(int32_t x, int32_t y, int32_t M) {
    return (int64_t(x) * y) % M;
}

// modulo inversion
template<typename T> T modulo_inv(const T& v, const T& M) {
    T vi; T g = gcd_ex(v, M, &vi);
    if (g != 1) vi /= g;
    return vi;
}
// integral type specializations, input is assumed to be normalized
template<typename I> I modulo_inv_int(I v, I M) {
    I vi; gcd_ex(v, M, &vi);
    modulo_normalize_int(&vi, M);
    return vi;
}
inline int64_t modulo_inv(int64_t v, int64_t M) { return modulo_inv_int(v, M); }
inline int32_t modulo_inv(int32_t v, int32_t M) { return modulo_inv_int(v, M); }

// modulo division
template<typename T> T modulo_div(const T& x, const T& y, const T& M) {
    return modulo_mul(x, modulo_inv(y, M), M);
}
// integral type specializations, input is assumed to be normalized
template<typename I>
I modulo_div_int(I x, I y, I M) {
    if (y != 0 && x % y == 0) return x / y; // fast path if y divides x
    I yi; I g = gcd_ex(y, M, &yi);
    if (g != 1) { // uh oh, y and M are not coprime, try common gcd
        g = altruct::math::gcd(x, g);
        x /= g; y /= g;
        gcd_ex(y, I(M / g), &yi); // == k
        // if (k != 1), there is no result, or more precisely,
        // the result will be `k` times bigger than it should.
    }
    modulo_normalize_int(&yi, I(M / g));
    return modulo_mul(x, yi, M);
}
inline int64_t modulo_div(int64_t x, int64_t y, int64_t M) { return modulo_div_int(x, y, M); }
inline int32_t modulo_div(int32_t x, int32_t y, int32_t M) { return modulo_div_int(x, y, M); }


// modulo storage type

namespace modulo_storage {
    enum type { INSTANCE, STATIC, CONSTANT };
}

template<typename T, int ID, int STORAGE_TYPE>
struct modulo_members;

template<typename T, int ID>
struct modulo_members<T, ID, modulo_storage::INSTANCE> {
    T _M;
    modulo_members(const T& _M = 0) : _M(_M) {}
    const T& M() const { return _M; }
    T& M() { return _M; }
};

template<typename T, int ID>
struct modulo_members<T, ID, modulo_storage::STATIC> {
    static T _M;
    modulo_members(const T& _M = 0) {}
    static T& M() { return _M; }
};

template<typename T, int ID>
struct modulo_members<T, ID, modulo_storage::CONSTANT> {
    modulo_members(const T& _M = 0) {}
    static T M() { return T(ID); }
};

/**
 * Modulo M arithmetics
 *
 * modulo<int, 3> - Z/3Z
 *
 * @param T - the underlying type
 * @param ID - ID of the modulo type (useful with modulo_storage::CONSTANT)
 * @param STORAGE_TYPE - whether M is a constant, static or instance member
 *    CONSTANT - Uses the ID template argument as M which is a constant.
 *      This allows compiler to employ optimized division by constant.
 *      If your moudlo is always say 1000000007, this is the way to go.
 *      This option can only be used when the type T is int.
 *    STATIC - Has a class static member for M. Separate for each <T, ID>.
 *      This allows to avoid having the same instance of M for each instance
 *      of modulo. Useful when having a large array of instances when M gets
 *      known only at run time, or when the type T is not int.
 *    INSTANCE - Each modulo instance consists of both value v and modulus M.
 *      This is not as time and space efficient as the above two, but is the
 *      prefered option when keeping an instance of M for each instance of
 *      modulo is not a problem.
 *      Note that operations between two instances (v1, m1) and (v2, m2) with
 *      different moduli are allowed and in such case m1 gets used as modulus.
 *      This is both for performance reasons (avoids a check) and convenience
 *      as one can do `modx(v, M) * int(u)` in which case the second operand
 *      gets resolved to `modx(u, 0)` which has an invalid modulus 0.
 */
template<typename T, int ID, int STORAGE_TYPE = modulo_storage::STATIC>
class modulo : public modulo_members<T, ID, STORAGE_TYPE> {
    typedef modulo_members<T, ID, STORAGE_TYPE> my_modulo_members;
public:
    T v;

    // construct from int, but only if T is not integral to avoid constructor clashing
    template <typename I = T, typename = std::enable_if_t<!std::is_integral<I>::value>>
    modulo(int v) : my_modulo_members(1), v(v) { if (STORAGE_TYPE != modulo_storage::INSTANCE) normalize(); }
    modulo() : my_modulo_members(1), v(zeroOf(this->M())) { if (STORAGE_TYPE != modulo_storage::INSTANCE) normalize(); }
    modulo(const T& v) : my_modulo_members(1), v(v) { if (STORAGE_TYPE != modulo_storage::INSTANCE) normalize(); }
    modulo(const T& v, const T& M) : my_modulo_members(M), v(v) { normalize(); }
    modulo(const modulo& rhs) : my_modulo_members(rhs.M()), v(rhs.v) {}

    void normalize() { modulo_normalize(&v, this->M()); }

    bool operator == (const modulo &rhs) const { return (v == rhs.v); }
    bool operator != (const modulo &rhs) const { return (v != rhs.v); }
    bool operator <  (const modulo &rhs) const { return (v <  rhs.v); }
    bool operator >  (const modulo &rhs) const { return (v >  rhs.v); }
    bool operator <= (const modulo &rhs) const { return (v <= rhs.v); }
    bool operator >= (const modulo &rhs) const { return (v >= rhs.v); }

    modulo  operator +  (const modulo &rhs) const { modulo t(*this); t += rhs; return t; }
    modulo  operator -  (const modulo &rhs) const { modulo t(*this); t -= rhs; return t; }
    modulo  operator -  ()                  const { modulo t(-v, this->M());   return t; }
    modulo  operator *  (const modulo &rhs) const { modulo t(*this); t *= rhs; return t; }
    modulo  operator /  (const modulo &rhs) const { modulo t(*this); t /= rhs; return t; }
    modulo  operator %  (const modulo &rhs) const { modulo t(*this); t %= rhs; return t; }

    modulo& operator += (const modulo &rhs) { v = modulo_add(v, rhs.v, this->M()); return *this; }
    modulo& operator -= (const modulo &rhs) { v = modulo_sub(v, rhs.v, this->M()); return *this; }
    modulo& operator *= (const modulo &rhs) { v = modulo_mul(v, rhs.v, this->M()); return *this; }
    modulo& operator /= (const modulo &rhs) { v = modulo_div(v, rhs.v, this->M()); return *this; }
    modulo& operator %= (const modulo &rhs) { v %= rhs.v;                          return *this; }

    modulo inv() const { return modulo_inv(v, this->M()); }
};

template<typename T>
using moduloX = modulo<T, 0, modulo_storage::INSTANCE>;

template<typename T, int ID>
T modulo_members<T, ID, modulo_storage::STATIC>::_M = castOf<T>(ID);

template<typename T, int ID, typename I>
struct castT<modulo<T, ID, modulo_storage::INSTANCE>, I> {
    typedef modulo<T, ID, modulo_storage::INSTANCE> mod;
    static mod of(const I& x) {
        return mod(castOf<T>(x));
    }
    static mod of(const mod& ref, const I& x) {
        return mod(castOf(ref.v, x % ref.M()), ref.M());
    }
};
template<typename T, int ID>
struct castT<modulo<T, ID, modulo_storage::INSTANCE>, modulo<T, ID, modulo_storage::INSTANCE>> : nopCastT<modulo<T, ID, modulo_storage::INSTANCE>>{};

template<typename T, int ID, int STORAGE_TYPE, typename I>
struct castT<modulo<T, ID, STORAGE_TYPE>, I> {
    typedef modulo<T, ID, STORAGE_TYPE> mod;
    static mod of(const I& x) {
        return mod(castOf<T>(x % mod::M()));
    }
    static mod of(const mod& ref, const I& x) {
        return mod(castOf(ref.v, x % ref.M()), ref.M());
    }
};
template<typename T, int ID, int STORAGE_TYPE>
struct castT<modulo<T, ID, STORAGE_TYPE>, modulo<T, ID, STORAGE_TYPE>> : nopCastT<modulo<T, ID, STORAGE_TYPE>>{};

template<typename T, int ID, int STORAGE_TYPE>
struct identityT<modulo<T, ID, STORAGE_TYPE>> {
    typedef modulo<T, ID, STORAGE_TYPE> mod;
    static mod of(const mod& x) {
        return mod(identityOf(x.v), x.M());
    }
};

template<typename T, int ID, int STORAGE_TYPE>
struct zeroT<modulo<T, ID, STORAGE_TYPE>> {
    typedef modulo<T, ID, STORAGE_TYPE> mod;
    static mod of(const mod& x) {
        return mod(zeroOf(x.v), x.M());
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
