#pragma once

#include "altruct/algorithm/math/base.h"
#include "altruct/structure/math/vector3d.h"

namespace altruct {
namespace math {

/**
 * Plane-based Geometric Algebra in 3D
 * Clifford Algebra Cl(3, 0, 1) a.k.a Geometric Algebra G(3, 0, 1) in 3D
 *
 *   orthogonal basis: {e0, e1, e2, e3}
 *     3 positive vectors: {e1, e2, e3}
 *     0 negative vectors: {}
 *     1 null     vector:  {e0}
 *
 *        e0 e1 e2 e3
 *     e0  0  0  0  0
 *     e1  0  1  0  0
 *     e2  0  0  1  0
 *     e3  0  0  0  1
 *
 *     e23 = e2e3 = -e3e2 = -e32
 *     e31 = e3e1 = -e1e3 = -e13
 *     e12 = e1e2 = -e2e1 = -e21
 *     e123 = e231 = e312 = -e321 = -e213 = -e132
 *
 *   elements:   {e0, {e1, e2, e3},   1   , {e23, e31, e12}, {e01, e02, e03}, e0123  ,  e123, {e032, e013, e021}}
 *   blades:     |    vector      | salar |   E-bivector   |   e-bivector   | pseudo |        trivector         |
 *   primitives: |    plane       |       |          line or screw          |        |          point           |
 *   operators:  |                |     rotor              |           translator    |                          |
 *
 *   dual:                     !a = {e123, triP, e0123, bie, biE, s,     e0,   v   }
 *                                 !{e0,   v,    s,     biE, bie, e0123, e123, triP}
 *   reverse:                  ~a = a.rev(); change sign of blade2 and blade3
 * > geometric product:        a * b
 * > wedge/meet/outer product: a ^ b
 * > dot/inner product:        a & b
 * > join (dual to meet):      a | b = !(!a ^ !b)
 *
 * > translator(tx, ty, tz) = 1 + e t / 2 = exp(-t e / 2)
 *              t = tx e1 + tx e2 + tx e3
 *
 * > rotor(line L, angle phi) = exp(-L phi/2)
 *
 * > dot(line L, point P) = plane through a point P and perpendicular to a line L
 *
 */
namespace pga {
template<typename T>
class zero {};

template<typename T>
class blade1 {
public:
    T e0; vector3d<T> v;

    blade1() {}
    explicit blade1(T e0) : e0(e0), v({ zeroOf(e0), zeroOf(e0) , zeroOf(e0) }) {}
    explicit blade1(vector3d<T> v) : e0(zeroOf(v.z)), v(v) {}
    blade1(T e0, vector3d<T> v) : e0(std::move(e0)), v(std::move(v)) {}

    blade1& operator += (const blade1& b) { e0 += b.e0; v += b.v; return *this; }
    blade1& operator -= (const blade1& b) { e0 -= b.e0; v -= b.v; return *this; }
    blade1& operator *= (const T& s) { e0 *= s; v *= s; return *this; }
    blade1& operator /= (const T& s) { e0 /= s; v /= s; return *this; }

    blade1  operator +  ()                const { return blade1(+e0, +v); }
    blade1  operator -  ()                const { return blade1(-e0, -v); }
    blade1  operator +  (const blade1& b) const { auto r = *this; r += b; return r; }
    blade1  operator -  (const blade1& b) const { auto r = *this; r -= b; return r; }
    blade1  operator *  (const T& s)      const { auto r = *this; r *= s; return r; }
    blade1  operator /  (const T& s)      const { auto r = *this; r /= s; return r; }

    blade1  operator ~  () const { return rev(); }

    blade1 rev() const { return blade1(e0, v); } // same
};

template<typename T>
class blade02 {
public:
    T s; vector3d<T> biE;

    blade02() {}
    explicit blade02(T s) : s(s), biE({ zeroOf(s), zeroOf(s) , zeroOf(s) }) {}
    explicit blade02(vector3d<T> biE) : s(zeroOf(biE.z)), biE(biE) {}
    blade02(T s, vector3d<T> biE) : s(std::move(s)), biE(std::move(biE)) {}

    blade02& operator += (const blade02& b) { s += b.s; biE += b.biE; return *this; }
    blade02& operator -= (const blade02& b) { s -= b.s; biE -= b.biE; return *this; }
    blade02& operator *= (const T& s) { this->s *= s; biE *= s; return *this; }
    blade02& operator /= (const T& s) { this->s /= s; biE /= s; return *this; }

    blade02  operator +  ()                 const { return blade02(+s, +biE); }
    blade02  operator -  ()                 const { return blade02(-s, -biE); }
    blade02  operator +  (const blade02& b) const { auto r = *this; r += b; return r; }
    blade02  operator -  (const blade02& b) const { auto r = *this; r -= b; return r; }
    blade02  operator *  (const T& s)       const { auto r = *this; r *= s; return r; }
    blade02  operator /  (const T& s)       const { auto r = *this; r /= s; return r; }

    blade02  operator ~  () const { return rev(); }

    blade02 rev() const { return blade02(s, -biE); } // negate blade2 part
};

template<typename T>
class blade24 {
public:
    vector3d<T> bie; T e0123;

    blade24() {}
    explicit blade24(T e0123) : bie({ zeroOf(e0123), zeroOf(e0123) , zeroOf(e0123) }), e0123(e0123) {}
    explicit blade24(vector3d<T> bie) : bie(bie), e0123(zeroOf(bie.z)) {}
    blade24(vector3d<T> bie, T e0123) : bie(std::move(bie)), e0123(std::move(e0123)) {}

    blade24& operator += (const blade24& b) { bie += b.bie; e0123 += b.e0123; return *this; }
    blade24& operator -= (const blade24& b) { bie -= b.bie; e0123 -= b.e0123; return *this; }
    blade24& operator *= (const T& s) { bie *= s; e0123 *= s; return *this; }
    blade24& operator /= (const T& s) { bie /= s; e0123 /= s; return *this; }

    blade24  operator +  ()                 const { return blade24(+bie, +e0123); }
    blade24  operator -  ()                 const { return blade24(-bie, -e0123); }
    blade24  operator +  (const blade24& b) const { auto r = *this; r += b; return r; }
    blade24  operator -  (const blade24& b) const { auto r = *this; r -= b; return r; }
    blade24  operator *  (const T& s)       const { auto r = *this; r *= s; return r; }
    blade24  operator /  (const T& s)       const { auto r = *this; r /= s; return r; }

    blade24  operator ~  () const { return rev(); }

    blade24 rev() const { return blade24(-bie, e0123); } // negate blade2 part
};

template<typename T>
class blade3 {
public:
    T e123; vector3d<T> triP;

    blade3() {}
    explicit blade3(T e123) : e123(e123), triP({ zeroOf(e123), zeroOf(e123) , zeroOf(e123) }) {}
    explicit blade3(vector3d<T> triP) : e123(zeroOf(triP.z)), triP(triP) {}
    blade3(T e123, vector3d<T> triP) : e123(std::move(e123)), triP(std::move(triP)) {}

    blade3& operator += (const blade3& b) { e123 += b.e123; triP += b.triP; return *this; }
    blade3& operator -= (const blade3& b) { e123 -= b.e123; triP -= b.triP; return *this; }
    blade3& operator *= (const T& s) { e123 *= s; triP *= s; return *this; }
    blade3& operator /= (const T& s) { e123 /= s; triP /= s; return *this; }

    blade3  operator +  ()                const { return blade3(+e123, +triP); }
    blade3  operator -  ()                const { return blade3(-e123, -triP); }
    blade3  operator +  (const blade3& b) const { auto r = *this; r += b; return r; }
    blade3  operator -  (const blade3& b) const { auto r = *this; r -= b; return r; }
    blade3  operator *  (const T& s)      const { auto r = *this; r *= s; return r; }
    blade3  operator /  (const T& s)      const { auto r = *this; r /= s; return r; }

    blade3  operator ~  () const { return rev(); }

    blade3 rev() const { return blade3(-e123, -triP); } // negate blade3 part
};

//-------------------------------------------------------------------------------

template<typename T>
class blade13 {
public:
    blade1<T> b1;
    blade3<T> b3;

    blade13() {}
    explicit blade13(blade1<T> b1) : b1(b1), b3(zeroOf(b1.e0)) {}
    explicit blade13(blade3<T> b3) : b1(zeroOf(b3.e123)), b3(b3) {}
    blade13(blade1<T> b1, blade3<T> b3) : b1(std::move(b1)), b3(std::move(b3)) {}

    blade13& operator += (const blade13& b) { b1 += b.b1; b3 += b.b3; return *this; }
    blade13& operator -= (const blade13& b) { b1 -= b.b1; b3 -= b.b3; return *this; }
    blade13& operator *= (const T& s) { b1 *= s; b3 *= s; return *this; }
    blade13& operator /= (const T& s) { b1 /= s; b3 /= s; return *this; }

    blade13  operator +  ()                 const { return blade13(+b1, +b3); }
    blade13  operator -  ()                 const { return blade13(-b1, -b3); }
    blade13  operator +  (const blade13& b) const { auto r = *this; r += b; return r; }
    blade13  operator -  (const blade13& b) const { auto r = *this; r -= b; return r; }
    blade13  operator *  (const T& s)       const { auto r = *this; r *= s; return r; }
    blade13  operator /  (const T& s)       const { auto r = *this; r /= s; return r; }

    blade13  operator ~  () const { return rev(); }

    blade13 rev() const { return blade13(b1.rev(), b3.rev()); }
};

template<typename T>
class blade024 {
public:
    blade02<T> b02;
    blade24<T> b24;

    blade024() {}
    explicit blade024(blade02<T> b02) : b02(b02), b24(zeroOf(b02.s)) {}
    explicit blade024(blade24<T> b24) : b02(zeroOf(b24.e0123)), b24(b24) {}
    blade024(blade02<T> b02, blade24<T> b24) : b02(std::move(b02)), b24(std::move(b24)) {}

    blade024& operator += (const blade024& b) { b02 += b.b02; b24 += b.b24; return *this; }
    blade024& operator -= (const blade024& b) { b02 -= b.b02; b24 -= b.b24; return *this; }
    blade024& operator *= (const T& s) { b02 *= s; b24 *= s; return *this; }
    blade024& operator /= (const T& s) { b02 /= s; b24 /= s; return *this; }

    blade024  operator +  ()                  const { return blade024(+b02, +b24); }
    blade024  operator -  ()                  const { return blade024(-b02, -b24); }
    blade024  operator +  (const blade024& b) const { auto r = *this; r += b; return r; }
    blade024  operator -  (const blade024& b) const { auto r = *this; r -= b; return r; }
    blade024  operator *  (const T& s)        const { auto r = *this; r *= s; return r; }
    blade024  operator /  (const T& s)        const { auto r = *this; r /= s; return r; }

    blade024  operator ~  () const { return rev(); }

    blade024 rev() const { return blade024(b02.rev(), b24.rev()); }
};

template<typename T>
class multivector {
public:
    blade13<T> b13;
    blade024<T> b024;

    multivector() {}
    multivector(blade13<T> b13, blade024<T> b024) : b13(std::move(b13)), b024(std::move(b024)) {}
    multivector(blade1<T> b1, blade02<T> b02, blade24<T> b24, blade3<T> b3) :
        b13(std::move(b1), std::move(b3)), b024(std::move(b02), std::move(b24)) {}

    multivector& operator += (const multivector& b) { b13 += b.b13; b024 += b.b024; return *this; }
    multivector& operator -= (const multivector& b) { b13 -= b.b13; b024 -= b.b024; return *this; }
    multivector& operator *= (const T& s) { b13 *= s; b024 *= s; return *this; }
    multivector& operator /= (const T& s) { b13 /= s; b024 /= s; return *this; }

    multivector  operator +  ()                     const { return multivector(+b13, +b024); }
    multivector  operator -  ()                     const { return multivector(-b13, -b024); }
    multivector  operator +  (const multivector& b) const { auto r = *this; r += b; return r; }
    multivector  operator -  (const multivector& b) const { auto r = *this; r -= b; return r; }
    multivector  operator *  (const T& s)           const { auto r = *this; r *= s; return r; }
    multivector  operator /  (const T& s)           const { auto r = *this; r /= s; return r; }

    multivector  operator ~  () const { return rev(); }

    multivector rev() const { return multivector(b13.rev(), b024.rev()); }
};

//-------------------------------------------------------------------------------

template<typename B, typename T> struct is_zero_blade_type { static const bool value = false; };
template<typename T> struct is_zero_blade_type<zero<T>, T> { static const bool value = true; };

template<typename B, typename T> struct is_primitive_blade_type { static const bool value = false; };
template<typename T> struct is_primitive_blade_type<blade1<T>, T> { static const bool value = true; };
template<typename T> struct is_primitive_blade_type<blade02<T>, T> { static const bool value = true; };
template<typename T> struct is_primitive_blade_type<blade24<T>, T> { static const bool value = true; };
template<typename T> struct is_primitive_blade_type<blade3<T>, T> { static const bool value = true; };
#define IF_PRIMITIVE_BLADE_TYPE(B, T) std::enable_if_t<is_primitive_blade_type<B, T>::value, bool> = true

template<typename B, typename T> struct is_composite_blade_type { static const bool value = false; };
template<typename T> struct is_composite_blade_type <blade13<T>, T> { static const bool value = true; };
template<typename T> struct is_composite_blade_type <blade024<T>, T> { static const bool value = true; };
template<typename T> struct is_composite_blade_type <multivector<T>, T> { static const bool value = true; };
#define IF_COMPOSITE_BLADE_TYPE(B, T) std::enable_if_t<is_composite_blade_type<B, T>::value, bool> = true

#define IF_NONZERO_BLADE_TYPE(B, T) std::enable_if_t<is_primitive_blade_type<B, T>::value || is_composite_blade_type<B, T>::value, bool> = true
#define IF_BLADE_TYPE(B, T) std::enable_if_t<is_zero_blade_type<B, T>::value || is_primitive_blade_type<B, T>::value || is_composite_blade_type<B, T>::value, bool> = true

//-------------------------------------------------------------------------------

// dual
template<typename T> blade3<T> operator ! (const blade1<T>& b1) { return blade3<T>(b1.e0, b1.v); }
template<typename T> blade24<T> operator ! (const blade02<T>& b02) { return blade24<T>(b02.biE, b02.s); }
template<typename T> blade02<T> operator ! (const blade24<T>& b24) { return blade02<T>(b24.e0123, b24.bie); }
template<typename T> blade1<T> operator ! (const blade3<T>& b3) { return blade1<T>(b3.e123, b3.triP); }
template<typename T> blade13<T> operator ! (const blade13<T>& b) { return blade13<T>(!b.b3, !b.b1); }
template<typename T> blade024<T> operator ! (const blade024<T>& b) { return blade024<T>(!b.b24, !b.b02); }
template<typename T> multivector<T> operator ! (const multivector<T>& m) { return multivector<T>(!m.b13, !m.b024); }

// get
template<typename B> struct get {};
template<typename T> struct get<zero<T>> {
    using type = T;
    static zero<T> b1(zero<T>) { return zero<T>(); }
    static zero<T> b02(zero<T>) { return zero<T>(); }
    static zero<T> b24(zero<T>) { return zero<T>(); }
    static zero<T> b3(zero<T>) { return zero<T>(); }
};
template<typename T> struct get<blade1<T>> {
    using type = T;
    static const blade1<T>& b1(const blade1<T>& b1) { return b1; }
    static zero<T> b02(const blade1<T>& b1) { return zero<T>(); }
    static zero<T> b24(const blade1<T>& b1) { return zero<T>(); }
    static zero<T> b3(const blade1<T>& b1) { return zero<T>(); }
};
template<typename T> struct get<blade3<T>> {
    using type = T;
    static zero<T> b1(const blade3<T>& b3) { return zero<T>(); }
    static zero<T> b02(const blade3<T>& b3) { return zero<T>(); }
    static zero<T> b24(const blade3<T>& b3) { return zero<T>(); }
    static const blade3<T>& b3(const blade3<T>& b3) { return b3; }
};
template<typename T> struct get<blade13<T>> {
    using type = T;
    static const blade1<T>& b1(const blade13<T>& b) { return b.b1; }
    static zero<T> b02(const blade13<T>& b) { return zero<T>(); }
    static zero<T> b24(const blade13<T>& b) { return zero<T>(); }
    static const blade3<T>& b3(const blade13<T>& b) { return b.b3; }
};
template<typename T> struct get<blade02<T>> {
    using type = T;
    static zero<T> b1(const blade02<T>& b02) { return zero<T>(); }
    static const blade02<T>& b02(const blade02<T>& b02) { return b02; }
    static zero<T> b24(const blade02<T>& b02) { return zero<T>(); }
    static zero<T> b3(const blade02<T>& b02) { return zero<T>(); }
};
template<typename T> struct get<blade24<T>> {
    using type = T;
    static zero<T> b1(const blade24<T>& b24) { return zero<T>(); }
    static zero<T> b02(const blade24<T>& b24) { return zero<T>(); }
    static const blade24<T>& b24(const blade24<T>& b24) { return b24; }
    static zero<T> b3(const blade24<T>& b24) { return zero<T>(); }
};
template<typename T> struct get<blade024<T>> {
    using type = T;
    static zero<T> b1(const blade024<T>& b) { return zero<T>(); }
    static const blade02<T>& b02(const blade024<T>& b) { return b.b02; }
    static const blade24<T>& b24(const blade024<T>& b) { return b.b24; }
    static zero<T> b3(const blade024<T>& b) { return zero<T>(); }
};
template<typename T> struct get<multivector<T>> {
    using type = T;
    static const blade1<T>& b1(const multivector<T>& m) { return m.b13.b1; }
    static const blade02<T>& b02(const multivector<T>& m) { return m.b024.b02; }
    static const blade24<T>& b24(const multivector<T>& m) { return m.b024.b24; }
    static const blade3<T>& b3(const multivector<T>& m) { return m.b13.b3; }
};

// combine
template<typename B, typename T, IF_BLADE_TYPE(B, T)> const B& combine_primitive(zero<T>, const B& b) { return b; }
template<typename B, typename T, IF_NONZERO_BLADE_TYPE(B, T)> const B& combine_primitive(const B& b, zero<T>) { return b; }
template<typename T> const blade13<T> combine_primitive(blade1<T> b1, blade3<T> b3) { return blade13<T>(std::move(b1), std::move(b3)); }
template<typename T> const blade024<T> combine_primitive(blade02<T> b02, blade24<T> b24) { return blade024<T>(std::move(b02), std::move(b24)); }

template<typename B, typename T, IF_BLADE_TYPE(B, T)> const B& combine_multivector(zero<T>, const B& b) { return b; }
template<typename B, typename T, IF_NONZERO_BLADE_TYPE(B, T)> const B& combine_multivector(const B& b, zero<T>) { return b; }
template<typename BL, typename BR, typename T = get<BL>::type, IF_NONZERO_BLADE_TYPE(BL, T), IF_NONZERO_BLADE_TYPE(BR, T)>
const multivector<T> combine_multivector(BL bl, BR br) { return multivector<T>(blade13<T>(std::move(bl)), blade024<T>(std::move(br))); }

// add
template<typename B, typename T, IF_BLADE_TYPE(B, T)> const B& operator + (zero<T>, const B& b) { return b; }
template<typename B, typename T, IF_NONZERO_BLADE_TYPE(B, T)> const B& operator + (const B& b, zero<T>) { return b; }

template<typename BL, typename BR, typename T = get<BL>::type, IF_NONZERO_BLADE_TYPE(BL, T), IF_NONZERO_BLADE_TYPE(BR, T)>
auto operator + (const BL& bl, const BR& br) {
    return combine_multivector(
        combine_primitive(
            get<BL>::b1(bl) + get<BR>::b1(br),
            get<BL>::b3(bl) + get<BR>::b3(br)),
        combine_primitive(
            get<BL>::b02(bl) + get<BR>::b02(br),
            get<BL>::b24(bl) + get<BR>::b24(br)));
}

// multiply
template<typename T> blade1<T> operator * (const T& s, const blade1<T>& b) {
    return blade1<T>(s * b.e0, s * b.v);
}
template<typename T> blade02<T> operator * (const T& s, const blade02<T>& b) {
    return blade02<T>(s * b.s, s * b.biE);
}
template<typename T> blade24<T> operator * (const T& s, const blade24<T>& b) {
    return blade24<T>(s * b.bie, s * b.e0123);
}
template<typename T> blade3<T> operator * (const T& s, const blade3<T>& b) {
    return blade3<T>(s * b.e123, s * b.triP);
}

template<typename T> blade024<T> operator * (const blade1<T>& a, const blade1<T>& b) {
    return blade02<T>(a.v & b.v, a.v ^ b.v) + blade24<T>(a.e0 * b.v - a.v * b.e0);
}
template<typename T> blade13<T> operator * (const blade1<T>& a, const blade02<T>& b) {
    return blade1<T>(a.e0 * b.s, a.v * b.s - (a.v ^ b.biE)) + blade3<T>(a.v & b.biE, -a.e0 * b.biE);
}
template<typename T> blade13<T> operator * (const blade1<T>& a, const blade24<T>& b) {
    return blade1<T>(-(a.v & b.bie)) + blade3<T>(a.v * b.e0123 + (a.v ^ b.bie));
}
template<typename T> blade024<T> operator * (const blade1<T>& a, const blade3<T>& b) {
    return blade02<T>(a.v * b.e123) + blade24<T>(-a.v ^ b.triP, a.e0 * b.e123 + (a.v & b.triP));
}
template<typename T> blade13<T> operator * (const blade02<T>& a, const blade1<T>& b) {
    return blade1<T>(a.s * b.e0, a.s * b.v - (a.biE ^ b.v)) + blade3<T>(a.biE & b.v, -a.biE * b.e0);
}
template<typename T> blade02<T> operator * (const blade02<T>& a, const blade02<T>& b) {
    return blade02<T>(a.s * b.s - (a.biE & b.biE), a.biE * b.s + a.s * b.biE - (a.biE ^ b.biE));
}
template<typename T> blade24<T> operator * (const blade02<T>& a, const blade24<T>& b) {
    return blade24<T>(a.s * b.bie - a.biE * b.e0123 - (a.biE ^ b.bie), a.s * b.e0123 + (a.biE & b.bie));
}
template<typename T> blade13<T> operator * (const blade02<T>& a, const blade3<T>& b) {
    return blade1<T>(a.biE & b.triP, -a.biE * b.e123) + blade3<T>(a.s * b.e123, a.s * b.triP - (a.biE ^ b.triP));
}
template<typename T> blade13<T> operator * (const blade24<T>& a, const blade1<T>& b) {
    return blade1<T>(a.bie & b.v) + blade3<T>(-a.e0123 * b.v - (a.bie ^ b.v));
}
template<typename T> blade24<T> operator * (const blade24<T>& a, const blade02<T>& b) {
    return blade24<T>(a.bie * b.s - a.e0123 * b.biE - (a.bie ^ b.biE), a.e0123 * b.s + (a.bie & b.biE));
}
template<typename T> zero<T> operator * (const blade24<T>& a, const blade24<T>& b) {
    return zero<T>();
}
template<typename T> blade13<T> operator * (const blade24<T>& a, const blade3<T>& b) {
    return blade1<T>(-a.e0123 * b.e123) + blade3<T>(-a.bie * b.e123);
}
template<typename T> blade024<T> operator * (const blade3<T>& a, const blade1<T>& b) {
    return blade02<T>(a.e123 * b.v) + blade24<T>(a.triP ^ b.v, -a.e123 * b.e0 - (a.triP & b.v));
}
template<typename T> blade13<T> operator * (const blade3<T>& a, const blade02<T>& b) {
    return blade1<T>(a.triP & b.biE, -a.e123 * b.biE) + blade3<T>(a.e123 * b.s, a.triP * b.s - (a.triP ^ b.biE));
}
template<typename T> blade13<T> operator * (const blade3<T>& a, const blade24<T>& b) {
    return blade1<T>(a.e123 * b.e0123) + blade3<T>(a.e123 * b.bie);
}
template<typename T> blade024<T> operator * (const blade3<T>& a, const blade3<T>& b) {
    return blade02<T>(-a.e123 * b.e123) + blade24<T>(a.triP * b.e123 - a.e123 * b.triP);
}

template<typename B, typename T, IF_BLADE_TYPE(B, T)> zero<T> operator * (zero<T> z, const B& b) { return z; }
template<typename B, typename T, IF_NONZERO_BLADE_TYPE(B, T)> zero<T> operator * (const B& b, zero<T> z) { return z; }
template<typename B, typename T, IF_PRIMITIVE_BLADE_TYPE(B, T)> auto operator * (const blade13<T>& a, const B& b) { return a.b1 * b + a.b3 * b; }
template<typename B, typename T, IF_PRIMITIVE_BLADE_TYPE(B, T)> auto operator * (const blade024<T>& a, const B& b) { return a.b02 * b + a.b24 * b; }
template<typename B, typename T, IF_PRIMITIVE_BLADE_TYPE(B, T)> auto operator * (const multivector<T>& a, const B& b) { return a.b13 * b + a.b024 * b; }
template<typename B, typename T, IF_NONZERO_BLADE_TYPE(B, T)> auto operator * (const B& a, const blade13<T>& b) { return a * b.b1 + a * b.b3; }
template<typename B, typename T, IF_NONZERO_BLADE_TYPE(B, T)> auto operator * (const B& a, const blade024<T>& b) { return a * b.b02 + a * b.b24; }
template<typename B, typename T, IF_NONZERO_BLADE_TYPE(B, T)> auto operator * (const B& a, const multivector<T>& b) { return a * b.b13 + a * b.b024; }

// wedge
//template<typename T> zero operator ^ (T& v, zero z) { return z; }
//template<typename T> zero operator ^ (zero z, T& v) { return z; }
//
//template<typename T> auto operator ^ (const blade1<T>& a, const blade1<T>& b) {
//    return blade02<T>(a.v ^ b.v) + blade24<T>(-a.v * b.e0 + a.e0 * b.v);
//}
//template<typename T> auto operator ^ (const blade1<T>& a, const blade02<T>& b) {
//    return blade1<T>(a.e0 * b.s, a.v * b.s) + blade3<T>(a.v & b.biE, -a.e0 * b.biE);
//}
//template<typename T> auto operator ^ (const blade1<T>& a, const blade24<T>& b) {
//    return blade3<T>(a.v ^ b.bie);
//}
//template<typename T> auto operator ^ (const blade1<T>& a, const blade3<T>& b) {
//    return blade24<T>(a.e0 * b.e123 + (a.v & b.triP));
//}
//template<typename T> auto operator ^ (const blade02<T>& a, const blade1<T>& b) {
//    return blade1<T>(a.s * b.e0, a.s * b.v) + blade3<T>(a.biE & b.v, -a.biE * b.e0);
//}
//template<typename T> auto operator ^ (const blade02<T>& a, const blade02<T>& b) {
//    return blade02<T>(a.s * b.s, a.biE * b.s + a.s * b.biE);
//}
//template<typename T> auto operator ^ (const blade02<T>& a, const blade24<T>& b) {
//    return blade24<T>(a.s * b.bie, a.s * b.e0123 + (a.biE & b.bie));
//}
//template<typename T> auto operator ^ (const blade02<T>& a, const blade3<T>& b) {
//    return blade3<T>(a.s * b.e123, a.s * b.triP);
//}
//template<typename T> auto operator ^ (const blade24<T>& a, const blade1<T>& b) {
//    return blade3<T>(-(a.bie ^ b.v));
//}
//template<typename T> auto operator ^ (const blade24<T>& a, const blade02<T>& b) {
//    return blade24<T>(a.bie * b.s, a.e0123 * b.s + (a.bie & b.biE));
//}
//template<typename T> auto operator ^ (const blade24<T>& a, const blade24<T>& b) {
//    return zero();
//}
//template<typename T> auto operator ^ (const blade24<T>& a, const blade3<T>& b) {
//    return zero();
//}
//template<typename T> auto operator ^ (const blade3<T>& a, const blade1<T>& b) {
//    return blade24<T>(-a.e123 * b.e0 - (a.triP & b.v));
//}
//template<typename T> auto operator ^ (const blade3<T>& a, const blade02<T>& b) {
//    return blade3<T>(a.e123 * b.s, a.triP * b.s);
//}
//template<typename T> auto operator ^ (const blade3<T>& a, const blade24<T>& b) {
//    return zero();
//}
//template<typename T> auto operator ^ (const blade3<T>& a, const blade3<T>& b) {
//    return zero();
//}

// TODO: wedge product overloads
// TODO: dot product overloads
// TODO: join overloads
// TODO: inverse
// TODO: conjugate product overloads (a b ~a)

//-------------------------------------------------------------------------------

// plane with normal `n` and distance `d / norm(n)` from the origin
template<typename T>
blade1<T> plane(vector3d<T> n, T d) {
    return blade1<T>(std::move(d), std::move(n));
}

// line with direction `l` through point `P`
template<typename T>
blade024<T> line(const vector3d<T>& l, const vector3d<T>& P) {
    auto zero = zeroOf(l.z);
    return blade024<T>(blade02<T>(zero, l), blade24<T>(l ^ P, zero));
}

// point with homogenous coordinate `h` at position `{x, y, z} / h`
template<typename T>
blade3<T> point(const vector3d<T>& P, T h) {
    return blade3<T>(std::move(h), -P);
}

//-------------------------------------------------------------------------------

// TODO: projection, rejection, intersection, join
// TODO: distance, angle

} // pga
} // math
} // altruct
