#pragma once

#include "altruct/algorithm/math/base.h"
#include "altruct/algorithm/math/intrinsic.h"

#include <type_traits>

namespace altruct {
namespace math {

// operations for non-integral types

template<typename V, typename T, typename std::enable_if_t<!std::is_integral<T>::value, bool> = true>
T modulo_normalize(const V& v, const T& M) { return castOf(M, v) % M; }
template<typename T, typename std::enable_if_t<!std::is_integral<T>::value, bool> = true>
T modulo_add(const T& x, const T& y, const T& M) { return (x + y) % M; }
template<typename T, typename std::enable_if_t<!std::is_integral<T>::value, bool> = true>
T modulo_sub(const T& x, const T& y, const T& M) { return (x - y) % M; }
template<typename T, typename std::enable_if_t<!std::is_integral<T>::value, bool> = true>
T modulo_neg(const T& v, const T& M) { return -v; }
template<typename T, typename std::enable_if_t<!std::is_integral<T>::value, bool> = true>
T modulo_mul(const T& x, const T& y, const T& M) { return (x * y) % M; }
template<typename T, typename std::enable_if_t<!std::is_integral<T>::value, bool> = true>
T modulo_gcd_ex(const T& n1, const T& n2, T& ni1, T& ni2) {
    return gcd_ex(n1, n2, &ni1, &ni2);
}
template<typename T, typename std::enable_if_t<!std::is_integral<T>::value, bool> = true>
T modulo_inv(const T& v, const T& M) {
    T vi; T g = gcd_ex(v, M, &vi);
    // `gcd_ex` produces `vi` and `Mi` (not stored) such that:
    // v * vi + M * Mi == g
    // for certain types the produced `g` might not be identity,
    // but the inverse still exists if `vi/g` and `Mi/g` exist:
    // v * vi/g + M * Mi/g == 1
    // for example, for polynomials, `g` may be a 0-degree polynomial
    // which is just a scalar, and dividing `vi` and `Mi` by `g` are
    // both well defined provided coefficients are invertible.
    // in such a case the inverse is simply `vi/g`:
    if (g != 1) vi /= g;
    return vi;
}
template<typename T, typename std::enable_if_t<!std::is_integral<T>::value, bool> = true>
T modulo_div(const T& x, const T& y, const T& M) {
    return modulo_mul(x, modulo_inv(y, M), M);
}

// operations for integral types; input is assumed to be normalized

template<typename U>
U modulo_add_uint(U x, U y, U M) { // y = M is allowed
    U r; bool of = add_overflow(x, y, &r);
    return (!of && r < M) ? r : (r - M);
}
template<typename S>
S modulo_add_int(S x, S y, S M) { // y = M is allowed
    // TODO: use add_overflow when MSVC adds intrinsics
    S r = x + y; bool of = (r < 0);
    return (!of && r < M) ? r : (r - M);
}
template<typename U, typename std::enable_if_t<std::is_unsigned<U>::value, bool> = true>
U modulo_add(U x, U y, U M) { return modulo_add_uint(x, y, M); }
template<typename S, typename std::enable_if_t<std::is_signed<S>::value, bool> = true>
S modulo_add(S x, S y, S M) { return modulo_add_int(x, y, M); }
template<typename I>
I modulo_sub_int(I x, I y, I M) { return modulo_add<I>(x, M - y, M); } // y = 0 is allowed
template<typename I, typename std::enable_if_t<std::is_integral<I>::value, bool> = true>
I modulo_sub(I x, I y, I M) { return modulo_sub_int(x, y, M); }
template<typename I>
I modulo_neg_int(I v, I M) { return (v == 0) ? v : (M - v); }
template<typename I, typename std::enable_if_t<std::is_integral<I>::value, bool> = true>
I modulo_neg(I v, I M) { return modulo_neg_int(v, M); }
template<typename I>
I modulo_mul_int_long(I x, I y, I M) {
    I r = 0;
    for (; y > 0; y >>= 1) {
        if (y & 1) r = modulo_add(r, x, M);
        x = modulo_add(x, x, M);
    }
    return r;
}
template<typename U, typename std::enable_if_t<std::is_unsigned<U>::value, bool> = true>
U modulo_mul(U x, U y, U M) {
    return (uint32_t(x) * uint32_t(y)) % uint32_t(M);
}
template<> inline uint32_t modulo_mul(uint32_t x, uint32_t y, uint32_t M) {
    return (uint64_t(x) * y) % M;
}
template<> inline uint64_t modulo_mul(uint64_t x, uint64_t y, uint64_t M) {
    return ((x >> 32) == 0 && (y >> 32) == 0) ? (x * uint32_t(y)) % M : modulo_mul_int_long(x, y, M);
}
template<typename S, typename std::enable_if_t<std::is_signed<S>::value, bool> = true>
S modulo_mul(S x, S y, S M) {
    using U = typename std::make_unsigned<S>::type;
    return S(modulo_mul(U(x), U(y), U(M)));
}
template<typename I>
I modulo_gcd_ex_int(I n1, I n2, I& ni1, I& ni2) {
    int s;
    I g = gcd_ex<I>(n1, n2, &ni1, &ni2, &s);
    if (s % 2 == 1 && ni1 != 0) ni1 += n2;
    if (s % 2 == 0 && ni2 != 0) ni2 += n1;
    return g;
}
template<typename I, typename std::enable_if_t<std::is_integral<I>::value, bool> = true>
I modulo_gcd_ex(I n1, I n2, I& ni1, I& ni2) { return modulo_gcd_ex_int(n1, n2, ni1, ni2); }
template<typename I>
I modulo_inv_int(I v, I M) {
    I vi; int s;
    I g = gcd_ex<I>(M, v, nullptr, &vi, &s);
    if (s % 2 == 0 && vi != 0) vi += M;
    return (g == 1) ? vi : 0;
}
template<typename I, typename std::enable_if_t<std::is_integral<I>::value, bool> = true>
I modulo_inv(I v, I M) { return modulo_inv_int(v, M); }
template<typename I>
I modulo_div_int(I x, I y, I M) {
    if (y != 0 && x % y == 0) return x / y; // fast path if y divides x
    I yi = modulo_inv_int(y, M); // modular inverse
    if (yi != 0) return modulo_mul(x, yi, M);
    // y and M are not coprime, try dividing by common gcd
    I g = altruct::math::gcd(altruct::math::gcd(x, y), M);
    return modulo_mul<I>(x / g, modulo_inv_int<I>(y / g, M / g), M);
}
template<typename I, typename std::enable_if_t<std::is_integral<I>::value, bool> = true>
I modulo_div(I x, I y, I M) { return modulo_div_int(x, y, M); }
// modulo_normalize for integral types
template<typename U, typename I, typename std::enable_if_t<std::is_unsigned<U>::value && std::is_integral<I>::value, bool> = true>
I modulo_normalize(U v, I M) {
    auto UM = static_cast<typename std::make_unsigned<I>::type>(M);
    if (v < UM) return static_cast<I>(v);
    return static_cast<I>(v % UM);
}
template<typename S, typename I, typename std::enable_if_t<std::is_signed<S>::value && std::is_integral<I>::value, bool> = true>
I modulo_normalize(S v, I M) {
    if (v < 0) return modulo_neg_int(modulo_normalize(-v, M), M);
    return modulo_normalize(static_cast<typename std::make_unsigned<S>::type>(v), M);
}


// modulo storage type

namespace modulo_storage {
    enum type { INSTANCE, STATIC, CONSTANT };
}

template<typename T, uint64_t ID, int STORAGE_TYPE>
struct modulo_members;

template<typename T, uint64_t ID>
struct modulo_members<T, ID, modulo_storage::INSTANCE> {
    T _M;
    modulo_members(const T& _M = T(1)) : _M(_M) {}
    const T& M() const { return _M; }
    T& M() { return _M; }
};

template<typename T, uint64_t ID>
struct modulo_members<T, ID, modulo_storage::STATIC> {
    static T _M;
    modulo_members(const T& _M = T(1)) {}
    static T& M() { return _M; }
};

template<typename T, uint64_t ID>
struct modulo_members<T, ID, modulo_storage::CONSTANT> {
    modulo_members(const T& _M = T(1)) {}
    static T M() { return castOf<T>(ID); }
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
 *      This option can only be used when the type T is integral.
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
template<typename T, uint64_t ID, int STORAGE_TYPE = modulo_storage::STATIC>
class modulo : public modulo_members<T, ID, STORAGE_TYPE> {
    typedef modulo_members<T, ID, STORAGE_TYPE> my_modulo_members;
public:
    T v;

    modulo() : my_modulo_members(), v(zeroOf(this->M())) {}
    modulo(const T& v_, const T& M_) : my_modulo_members(M_), v(modulo_normalize(v_, this->M())) {}
    modulo(const T& v_) : my_modulo_members(), v((STORAGE_TYPE != modulo_storage::INSTANCE) ? modulo_normalize(v_, this->M()) : v_) {}
    // construct from a different type I
    template<typename I, typename Enable = std::enable_if_t<!std::is_same<T, I>::value, bool>>
    modulo(const I& v_, const T& M_) : my_modulo_members(M_), v(modulo_normalize(v_, this->M())) {}
    template<typename I, typename Enable = std::enable_if_t<!std::is_same<T, I>::value && STORAGE_TYPE != modulo_storage::INSTANCE, bool>>
    modulo(const I& v_) : my_modulo_members(), v(modulo_normalize(v_, this->M())) {}

    bool operator == (const modulo &rhs) const { return (v == rhs.v); }
    bool operator != (const modulo &rhs) const { return (v != rhs.v); }
    bool operator <  (const modulo &rhs) const { return (v <  rhs.v); }
    bool operator >  (const modulo &rhs) const { return (v >  rhs.v); }
    bool operator <= (const modulo &rhs) const { return (v <= rhs.v); }
    bool operator >= (const modulo &rhs) const { return (v >= rhs.v); }

    modulo  operator +  (const modulo &rhs) const { auto t = *this; t += rhs; return t; }
    modulo  operator -  (const modulo &rhs) const { auto t = *this; t -= rhs; return t; }
    modulo  operator -  ()                  const { return neg(); }
    modulo  operator *  (const modulo &rhs) const { auto t = *this; t *= rhs; return t; }
    modulo  operator /  (const modulo &rhs) const { auto t = *this; t /= rhs; return t; }
    modulo  operator %  (const modulo &rhs) const { auto t = *this; t %= rhs; return t; }

    modulo& operator += (const modulo &rhs) { v = modulo_add(v, rhs.v, this->M()); return *this; }
    modulo& operator -= (const modulo &rhs) { v = modulo_sub(v, rhs.v, this->M()); return *this; }
    modulo& operator *= (const modulo &rhs) { v = modulo_mul(v, rhs.v, this->M()); return *this; }
    modulo& operator /= (const modulo &rhs) { v = modulo_div(v, rhs.v, this->M()); return *this; }
    modulo& operator %= (const modulo &rhs) { v %= rhs.v;                          return *this; }

    modulo neg() const { auto t = *this; t.v = modulo_neg(v, this->M()); return t; }
    modulo inv() const { auto t = *this; t.v = modulo_inv(v, this->M()); return t; }
};

template<typename T>
using moduloX = modulo<T, 0, modulo_storage::INSTANCE>;

template<typename T, uint64_t ID>
T modulo_members<T, ID, modulo_storage::STATIC>::_M = castOf<T>(ID);

template<typename T, uint64_t ID, int STORAGE_TYPE, typename I>
struct castT<modulo<T, ID, STORAGE_TYPE>, I> {
    typedef modulo<T, ID, STORAGE_TYPE> mod;
    static mod of(const I& v) {
        return mod(v);
    }
    static mod of(const mod& ref, const I& v) {
        return mod(v, ref.M());
    }
};
template<typename T, uint64_t ID, int STORAGE_TYPE>
struct castT<modulo<T, ID, STORAGE_TYPE>, modulo<T, ID, STORAGE_TYPE>> : nopCastT<modulo<T, ID, STORAGE_TYPE>>{};

template<typename T, uint64_t ID, int STORAGE_TYPE>
struct identityT<modulo<T, ID, STORAGE_TYPE>> {
    typedef modulo<T, ID, STORAGE_TYPE> mod;
    static mod of(const mod& x) {
        return mod(identityOf(x.v), x.M());
    }
};

template<typename T, uint64_t ID, int STORAGE_TYPE>
struct zeroT<modulo<T, ID, STORAGE_TYPE>> {
    typedef modulo<T, ID, STORAGE_TYPE> mod;
    static mod of(const mod& x) {
        return mod(zeroOf(x.v), x.M());
    }
};

template<typename T, uint64_t ID, int STORAGE_TYPE>
struct hasherT<modulo<T, ID, STORAGE_TYPE>> {
    typedef modulo<T, ID, STORAGE_TYPE> mod;
    size_t operator()(const mod& x) const {
        return hasherT<T>()(x.v);
    }
};

template<typename T>
T modT(T v, const T& M) { return modulo_normalize(v, M); }

template<typename T, typename I>
T modulo_power(const T& x, const I& y, const T& M) {
    return powT(moduloX<T>(x, M), y).v;
}

} // math
} // altruct
