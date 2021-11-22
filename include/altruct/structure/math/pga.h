#pragma once

#include "altruct/algorithm/math/base.h"
#include "altruct/structure/math/vector3d.h"

namespace altruct {
namespace math {

/**
 * Plane-based Geometric Algebra in 3D
 * Clifford Algebra Cl(3, 0, 1) a.k.a Geometric Algebra G(3, 0, 1) in 3D
 *
 *   orthogonal basis:
 *     3 positive vectors: {e1, e2, e3}
 *     0 negative vectors: {}
 *     1 null     vector:  {e0}
 *
 *   multiplication table:
 *          e0    e1   e2   e3
 *     e0    0   e01  e02  e03
 *     e1  -e01   1   e12 -e31
 *     e2  -e02 -e12   1   e23
 *     e3  -e03  e31 -e23   1
 *
 *   elements:   |e0, {e1, e2, e3}|   1   | {e23, e31, e12}| {e01, e02, e03}| e0123  |  e123, {e032, e013, e021}|
 *    square:    | 0   +1  +1  +1 |  +1   |  -1   -1   -1  |   0    0    0  |  0     |   -1     0     0     0   |
 *   ------------+----------------+-------+----------------+----------------+--------+--------------------------+
 *   blades:     |    vector      | salar |   E-bivector   |   e-bivector   | pseudo |        trivector         |
 *    c++ class: |    blade1      | blade0|     blade2E    |     blade2e    | blade4 |          blade3          |
 *   ------------+----------------+-------+----------------+----------------+--------+--------------------------+
 *   primitives: |    plane       |       |          line or screw          |        |          point           |
 *    c++ class: |    blade1      |       |             blade22             |        |          blade3          |
 *   ------------+----------------+-------+----------------+----------------+--------+--------------------------+
 *   operators:  |                |     rotor              |           translator    |                          |
 *    c++ class: |                |    blade02             |             blade24     |                          |
 *   ------------+----------------+-------+----------------+----------------+--------+--------------------------+
 *
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
 */
namespace pga {
template<typename T>
vector3d<T> make_zero_vec(const T& v) { return vector3d<T>(zeroOf(v), zeroOf(v), zeroOf(v)); }

#define PGA_CONSTRUCTORS_1(B, T, p)                                               \
    B() {}                                                                        \
    explicit B(T p) : p(std::move(p)) {}

#define PGA_CONSTRUCTORS_2(B, TL, pl, zl, TR, pr, zr)                             \
    B() {}                                                                        \
    explicit B(const TL& pl) : pl(pl), pr(zr) {}                                  \
    explicit B(const TR& pr) : pl(zl), pr(pr) {}                                  \
    B(TL pl, TR pr) : pl(std::move(pl)), pr(std::move(pr)) {}

#define PGA_CLOSED_OPERATORS_1(B, T, p)                                           \
    B& operator += (const B& rhs) { p += rhs.p; return *this; }                   \
    B& operator -= (const B& rhs) { p -= rhs.p; return *this; }                   \
    B& operator *= (const T& rhs) { p *= rhs; return *this; }                     \
    B& operator /= (const T& rhs) { p /= rhs; return *this; }                     \
                                                                                  \
    B  operator +  ()             const { return B(+p); }                         \
    B  operator -  ()             const { return B(-p); }                         \
    B  operator +  (const B& rhs) const { auto r = *this; r += rhs; return r; }   \
    B  operator -  (const B& rhs) const { auto r = *this; r -= rhs; return r; }   \
    B  operator *  (const T& rhs) const { auto r = *this; r *= rhs; return r; }   \
    B  operator /  (const T& rhs) const { auto r = *this; r /= rhs; return r; }   \
                                                                                  \
    B  operator ~  () const { return rev(); }

#define PGA_CLOSED_OPERATORS_2(B, T, pl, pr)                                      \
    B& operator += (const B& rhs) { pl += rhs.pl; pr += rhs.pr; return *this; }   \
    B& operator -= (const B& rhs) { pl -= rhs.pl; pr -= rhs.pr; return *this; }   \
    B& operator *= (const T& rhs) { pl *= rhs; pr *= rhs; return *this; }         \
    B& operator /= (const T& rhs) { pl /= rhs; pr /= rhs; return *this; }         \
                                                                                  \
    B  operator +  ()             const { return B(+pl, +pr); }                   \
    B  operator -  ()             const { return B(-pl, -pr); }                   \
    B  operator +  (const B& rhs) const { auto r = *this; r += rhs; return r; }   \
    B  operator -  (const B& rhs) const { auto r = *this; r -= rhs; return r; }   \
    B  operator *  (const T& rhs) const { auto r = *this; r *= rhs; return r; }   \
    B  operator /  (const T& rhs) const { auto r = *this; r /= rhs; return r; }   \
                                                                                  \
    B  operator ~  () const { return rev(); }

#define PGA_COMPOSITE_GETTERS(B, pl, pr)                                          \
    const auto& first() const { return pl; }                                      \
    const auto& second() const { return pr; }                                     \
                                                                                  \
    auto rev() const { return B(pl.rev(), pr.rev()); }                            \

template<typename T>
class zero {};

template<typename T>
class blade0 {
public:
    T s;

    PGA_CONSTRUCTORS_1(blade0, T, s);
    PGA_CLOSED_OPERATORS_1(blade0, T, s);

    blade0 rev() const { return blade0(s); } // same
};

template<typename T>
class blade1 {
public:
    T e0; vector3d<T> v;

    PGA_CONSTRUCTORS_2(blade1, T, e0, zeroOf(v.z), vector3d<T>, v, make_zero_vec(e0));
    PGA_CLOSED_OPERATORS_2(blade1, T, e0, v);

    blade1 rev() const { return blade1(e0, v); } // same
};

template<typename T>
class blade2E {
public:
    vector3d<T> biE;

    PGA_CONSTRUCTORS_1(blade2E, vector3d<T>, biE);
    PGA_CLOSED_OPERATORS_1(blade2E, T, biE);

    blade2E rev() const { return blade2E(-biE); } // negate blade2 part
};

template<typename T>
class blade2e {
public:
    vector3d<T> bie;

    PGA_CONSTRUCTORS_1(blade2e, vector3d<T>, bie);
    PGA_CLOSED_OPERATORS_1(blade2e, T, bie);

    blade2e rev() const { return blade2e(-bie); } // negate blade2 part
};

template<typename T>
class blade3 {
public:
    T e123; vector3d<T> triP;

    PGA_CONSTRUCTORS_2(blade3, T, e123, zeroOf(triP.z), vector3d<T>, triP, make_zero_vec(e123));
    PGA_CLOSED_OPERATORS_2(blade3, T, e123, triP);

    blade3 rev() const { return blade3(-e123, -triP); } // negate blade3 part
};

template<typename T>
class blade4 {
public:
    T e0123;

    PGA_CONSTRUCTORS_1(blade4, T, e0123);
    PGA_CLOSED_OPERATORS_1(blade4, T, e0123);

    blade4 rev() const { return blade4(e0123); } // same
};

//-------------------------------------------------------------------------------

template<typename T>
class blade02 {
public:
    blade0<T> b0;
    blade2E<T> b2E;

    PGA_CONSTRUCTORS_2(blade02, blade0<T>, b0, zeroOf(b2E.biE.z), blade2E<T>, b2E, zeroOf(b0.s));
    PGA_CLOSED_OPERATORS_2(blade02, T, b0, b2E);
    PGA_COMPOSITE_GETTERS(blade02, b0, b2E);
};

template<typename T>
class blade22 {
public:
    blade2E<T> b2E;
    blade2e<T> b2e;

    PGA_CONSTRUCTORS_2(blade22, blade2E<T>, b2E, zeroOf(b2e.bie.z), blade2e<T>, b2e, zeroOf(b2E.biE.z));
    PGA_CLOSED_OPERATORS_2(blade22, T, b2E, b2e);
    PGA_COMPOSITE_GETTERS(blade22, b2E, b2e);
};

template<typename T>
class blade24 {
public:
    blade2e<T> b2e;
    blade4<T> b4;

    PGA_CONSTRUCTORS_2(blade24, blade2e<T>, b2e, zeroOf(b4.e0123), blade4<T>, b4, zeroOf(b2e.bie.z));
    PGA_CLOSED_OPERATORS_2(blade24, T, b2e, b4);
    PGA_COMPOSITE_GETTERS(blade24, b2e, b4);
};

template<typename T>
class blade024 {
public:
    blade02<T> b02;
    blade24<T> b24;

    explicit blade024(const blade0<T>& b0) : b02(b0), b24(blade4<T>(zeroOf(b0.s))) {}
    explicit blade024(const blade2E<T>& b2E) : b02(b2E), b24(blade4<T>(zeroOf(b2E.biE.z))) {}
    explicit blade024(const blade2e<T>& b2e) : b02(blade0<T>(zeroOf(b2e.bie.z))), b24(b2e) {}
    explicit blade024(const blade4<T>& b4) : b02(blade0<T>(zeroOf(b4.e0123))), b24(b4) {}
    explicit blade024(const blade22<T>& c) : b02(c.b2E), b24(c.b2e) {}
    PGA_CONSTRUCTORS_2(blade024, blade02<T>, b02, blade0<T>(zeroOf(b24.b4.e0123)), blade24<T>, b24, blade4<T>(zeroOf(b02.b0.s)));
    PGA_CLOSED_OPERATORS_2(blade024, T, b02, b24);
    PGA_COMPOSITE_GETTERS(blade024, b02, b24);
};

template<typename T>
class blade13 {
public:
    blade1<T> b1;
    blade3<T> b3;

    PGA_CONSTRUCTORS_2(blade13, blade1<T>, b1, zeroOf(b3.e123), blade3<T>, b3, zeroOf(b1.e0));
    PGA_CLOSED_OPERATORS_2(blade13, T, b1, b3);
    PGA_COMPOSITE_GETTERS(blade13, b1, b3);
};

template<typename T>
class multivector {
public:
    blade024<T> b024;
    blade13<T> b13;

    multivector(blade0<T> b0, blade1<T> b1, blade2E<T> b2E, blade2e<T> b2e, blade3<T> b3, blade4<T> b4) :
        b13(std::move(b1), std::move(b3)), b024(blade02<T>(std::move(b0), std::move(b2E)), blade24<T>(std::move(b2e), std::move(b4))) {}
    multivector(blade1<T> b1, blade0<T> b0, blade2E<T> b2E, blade2e<T> b2e, blade4<T> b4, blade3<T> b3) :
        b13(std::move(b1), std::move(b3)), b024(blade02<T>(std::move(b0), std::move(b2E)), blade24<T>(std::move(b2e), std::move(b4))) {}
    multivector(blade1<T> b1, blade02<T> b02, blade24<T> b24, blade3<T> b3) :
        b13(std::move(b1), std::move(b3)), b024(std::move(b02), std::move(b24)) {}
    PGA_CONSTRUCTORS_2(multivector, blade024<T>, b024, blade02<T>(blade0<T>(zeroOf(b13.b1.e0))), blade13<T>, b13, blade1<T>(zeroOf(b024.b02.b0.s)));
    PGA_CLOSED_OPERATORS_2(multivector, T, b024, b13);
    PGA_COMPOSITE_GETTERS(multivector, b024, b13);
};

//-------------------------------------------------------------------------------

template<typename B, typename T> struct is_zero_blade_type { static const bool value = false; };
template<typename T> struct is_zero_blade_type<zero<T>, T> { static const bool value = true; };

template<typename B, typename T> struct is_primitive_blade_type { static const bool value = false; };
template<typename T> struct is_primitive_blade_type<blade0<T>, T> { static const bool value = true; };
template<typename T> struct is_primitive_blade_type<blade1<T>, T> { static const bool value = true; };
template<typename T> struct is_primitive_blade_type<blade2E<T>, T> { static const bool value = true; };
template<typename T> struct is_primitive_blade_type<blade2e<T>, T> { static const bool value = true; };
template<typename T> struct is_primitive_blade_type<blade3<T>, T> { static const bool value = true; };
template<typename T> struct is_primitive_blade_type<blade4<T>, T> { static const bool value = true; };
#define IF_PRIMITIVE_BLADE_TYPE(B, T) std::enable_if_t<is_primitive_blade_type<B, T>::value, bool> = true

template<typename B, typename T> struct is_composite_blade_type { static const bool value = false; };
template<typename T> struct is_composite_blade_type <blade02<T>, T> { static const bool value = true; };
template<typename T> struct is_composite_blade_type <blade22<T>, T> { static const bool value = true; };
template<typename T> struct is_composite_blade_type <blade24<T>, T> { static const bool value = true; };
template<typename T> struct is_composite_blade_type <blade024<T>, T> { static const bool value = true; };
template<typename T> struct is_composite_blade_type <blade13<T>, T> { static const bool value = true; };
template<typename T> struct is_composite_blade_type <multivector<T>, T> { static const bool value = true; };
#define IF_COMPOSITE_BLADE_TYPE(B, T) std::enable_if_t<is_composite_blade_type<B, T>::value, bool> = true

#define IF_NONZERO_BLADE_TYPE(B, T) std::enable_if_t<is_primitive_blade_type<B, T>::value || is_composite_blade_type<B, T>::value, bool> = true
#define IF_BLADE_TYPE(B, T) std::enable_if_t<is_zero_blade_type<B, T>::value || is_primitive_blade_type<B, T>::value || is_composite_blade_type<B, T>::value, bool> = true

//-------------------------------------------------------------------------------

// dual
template<typename T> blade4<T> operator ! (const blade0<T>& b0) { return blade4<T>(b0.s); }
template<typename T> blade3<T> operator ! (const blade1<T>& b1) { return blade3<T>(b1.e0, b1.v); }
template<typename T> blade2e<T> operator ! (const blade2E<T>& b2E) { return blade2e<T>(b2E.biE); }
template<typename T> blade2E<T> operator ! (const blade2e<T>& b2e) { return blade2E<T>(b2e.bie); }
template<typename T> blade1<T> operator ! (const blade3<T>& b3) { return blade1<T>(b3.e123, b3.triP); }
template<typename T> blade0<T> operator ! (const blade4<T>& b4) { return blade0<T>(b4.e0123); }
template<typename T> blade24<T> operator ! (const blade02<T>& c) { return blade24<T>(!c.b2E, !c.b0); }
template<typename T> blade22<T> operator ! (const blade22<T>& c) { return blade22<T>(!c.b2e, !c.b2E); }
template<typename T> blade02<T> operator ! (const blade24<T>& c) { return blade02<T>(!c.b4, !c.b2e); }
template<typename T> blade024<T> operator ! (const blade024<T>& c) { return blade024<T>(!c.b24, !c.b02); }
template<typename T> blade13<T> operator ! (const blade13<T>& c) { return blade13<T>(!c.b3, !c.b1); }
template<typename T> multivector<T> operator ! (const multivector<T>& c) { return multivector<T>(!c.b024, !c.b13); }

// get
template<typename B> struct get {};
template<typename T> struct get<zero<T>> {
    using type = T;
    static zero<T> b0(zero<T>) { return zero<T>(); }
    static zero<T> b1(zero<T>) { return zero<T>(); }
    static zero<T> b2E(zero<T>) { return zero<T>(); }
    static zero<T> b2e(zero<T>) { return zero<T>(); }
    static zero<T> b3(zero<T>) { return zero<T>(); }
    static zero<T> b4(zero<T>) { return zero<T>(); }
};
template<typename T> struct get<blade0<T>> {
    using type = T;
    static const blade0<T>& b0(const blade0<T>& b0) { return b0; }
    static zero<T> b1(const blade0<T>& b0) { return zero<T>(); }
    static zero<T> b2E(const blade0<T>& b0) { return zero<T>(); }
    static zero<T> b2e(const blade0<T>& b0) { return zero<T>(); }
    static zero<T> b3(const blade0<T>& b0) { return zero<T>(); }
    static zero<T> b4(const blade0<T>& b0) { return zero<T>(); }
};
template<typename T> struct get<blade1<T>> {
    using type = T;
    static zero<T> b0(const blade1<T>& b1) { return zero<T>(); }
    static const blade1<T>& b1(const blade1<T>& b1) { return b1; }
    static zero<T> b2E(const blade1<T>& b1) { return zero<T>(); }
    static zero<T> b2e(const blade1<T>& b1) { return zero<T>(); }
    static zero<T> b3(const blade1<T>& b1) { return zero<T>(); }
    static zero<T> b4(const blade1<T>& b1) { return zero<T>(); }
};
template<typename T> struct get<blade2E<T>> {
    using type = T;
    static zero<T> b0(const blade2E<T>& b2E) { return zero<T>(); }
    static zero<T> b1(const blade2E<T>& b2E) { return zero<T>(); }
    static const blade2E<T>& b2E(const blade2E<T>& b2E) { return b2E; }
    static zero<T> b2e(const blade2E<T>& b2E) { return zero<T>(); }
    static zero<T> b3(const blade2E<T>& b2E) { return zero<T>(); }
    static zero<T> b4(const blade2E<T>& b2E) { return zero<T>(); }
};
template<typename T> struct get<blade2e<T>> {
    using type = T;
    static zero<T> b0(const blade2e<T>& b2e) { return zero<T>(); }
    static zero<T> b1(const blade2e<T>& b2e) { return zero<T>(); }
    static zero<T> b2E(const blade2e<T>& b2e) { return zero<T>(); }
    static const blade2e<T>& b2e(const blade2e<T>& b2e) { return b2e; }
    static zero<T> b3(const blade2e<T>& b2e) { return zero<T>(); }
    static zero<T> b4(const blade2e<T>& b2e) { return zero<T>(); }
};
template<typename T> struct get<blade3<T>> {
    using type = T;
    static zero<T> b0(const blade3<T>& b3) { return zero<T>(); }
    static zero<T> b1(const blade3<T>& b3) { return zero<T>(); }
    static zero<T> b2E(const blade3<T>& b3) { return zero<T>(); }
    static zero<T> b2e(const blade3<T>& b3) { return zero<T>(); }
    static const blade3<T>& b3(const blade3<T>& b3) { return b3; }
    static zero<T> b4(const blade3<T>& b3) { return zero<T>(); }
};
template<typename T> struct get<blade4<T>> {
    using type = T;
    static zero<T> b0(const blade4<T>& b4) { return zero<T>(); }
    static zero<T> b1(const blade4<T>& b4) { return zero<T>(); }
    static zero<T> b2E(const blade4<T>& b4) { return zero<T>(); }
    static zero<T> b2e(const blade4<T>& b4) { return zero<T>(); }
    static zero<T> b3(const blade4<T>& b4) { return zero<T>(); }
    static const blade4<T>& b4(const blade4<T>& b4) { return b4; }
};
template<typename T> struct get<blade02<T>> {
    using type = T;
    static const blade0<T>& b0(const blade02<T>& c) { return c.b0; }
    static zero<T> b1(const blade02<T>& c) { return zero<T>(); }
    static const blade2E<T>& b2E(const blade02<T>& c) { return c.b2E; }
    static zero<T> b2e(const blade02<T>& c) { return zero<T>(); }
    static zero<T> b3(const blade02<T>& c) { return zero<T>(); }
    static zero<T> b4(const blade02<T>& c) { return zero<T>(); }
};
template<typename T> struct get<blade22<T>> {
    using type = T;
    static zero<T> b0(const blade22<T>& c) { return zero<T>(); }
    static zero<T> b1(const blade22<T>& c) { return zero<T>(); }
    static const blade2E<T>& b2E(const blade22<T>& c) { return c.b2E; }
    static const blade2e<T>& b2e(const blade22<T>& c) { return c.b2e; }
    static zero<T> b3(const blade22<T>& c) { return zero<T>(); }
    static zero<T> b4(const blade22<T>& c) { return zero<T>(); }
}; 
template<typename T> struct get<blade24<T>> {
    using type = T;
    static zero<T> b0(const blade24<T>& c) { return zero<T>(); }
    static zero<T> b1(const blade24<T>& c) { return zero<T>(); }
    static zero<T> b2E(const blade24<T>& c) { return zero<T>(); }
    static const blade2e<T>& b2e(const blade24<T>& c) { return c.b2e; }
    static zero<T> b3(const blade24<T>& c) { return zero<T>(); }
    static const blade4<T>& b4(const blade24<T>& c) { return c.b4; }
};
template<typename T> struct get<blade024<T>> {
    using type = T;
    static const blade0<T>& b0(const blade024<T>& c) { return c.b02.b0; }
    static zero<T> b1(const blade024<T>& c) { return zero<T>(); }
    static const blade2E<T>& b2E(const blade024<T>& c) { return c.b02.b2E; }
    static const blade2e<T>& b2e(const blade024<T>& c) { return c.b24.b2e; }
    static zero<T> b3(const blade024<T>& c) { return zero<T>(); }
    static const blade4<T>& b4(const blade024<T>& c) { return c.b24.b4; }
};
template<typename T> struct get<blade13<T>> {
    using type = T;
    static zero<T> b0(const blade13<T>& c) { return zero<T>(); }
    static const blade1<T>& b1(const blade13<T>& c) { return c.b1; }
    static zero<T> b2E(const blade13<T>& c) { return zero<T>(); }
    static zero<T> b2e(const blade13<T>& c) { return zero<T>(); }
    static const blade3<T>& b3(const blade13<T>& c) { return c.b3; }
    static zero<T> b4(const blade13<T>& c) { return zero<T>(); }
};
template<typename T> struct get<multivector<T>> {
    using type = T;
    static const blade0<T>& b0(const multivector<T>& c) { return c.b024.b02.b0; }
    static const blade1<T>& b1(const multivector<T>& c) { return c.b13.b1; }
    static const blade2E<T>& b2E(const multivector<T>& c) { return c.b024.b02.b2E; }
    static const blade2e<T>& b2e(const multivector<T>& c) { return c.b024.b24.b2e; }
    static const blade3<T>& b3(const multivector<T>& c) { return c.b13.b3; }
    static const blade4<T>& b4(const multivector<T>& c) { return c.b024.b24.b4; }
};

// combine
template<typename T> zero<T> combine024(zero<T>, zero<T>, zero<T>, zero<T>) { return zero<T>(); }
template<typename T> const blade0<T>& combine024(const blade0<T>& b0, zero<T>, zero<T>, zero<T>) { return b0; }
template<typename T> const blade2E<T>& combine024(zero<T>, const blade2E<T>& b2E, zero<T>, zero<T>) { return b2E; }
template<typename T> const blade2e<T>& combine024(zero<T>, zero<T>, const blade2e<T>& b2e, zero<T>) { return b2e; }
template<typename T> const blade4<T>& combine024(zero<T>, zero<T>, zero<T>, const blade4<T>& b4) { return b4; }
template<typename T> blade02<T> combine024(const blade0<T>& b0, const blade2E<T>& b2E, zero<T>, zero<T>) { return blade02<T>(b0, b2E); }
template<typename T> blade22<T> combine024(zero<T>, const blade2E<T>& b2E, const blade2e<T>& b2e, zero<T>) { return blade22<T>(b2E, b2e); }
template<typename T> blade24<T> combine024(zero<T>, zero<T>, const blade2e<T>& b2e, const blade4<T>& b4) { return blade24<T>(b2e, b4); }
template<typename T> blade024<T> combine024(zero<T>, const blade2E<T>& b2E, zero<T>, const blade4<T>& b4) { return blade024<T>(blade02<T>(b2E), blade24<T>(b4)); }
template<typename T> blade024<T> combine024(const blade0<T>& b0, zero<T>, zero<T>, const blade4<T>& b4) { return blade024<T>(blade02<T>(b0), blade24<T>(b4)); }
template<typename T> blade024<T> combine024(const blade0<T>& b0, zero<T>, const blade2e<T>& b2e, zero<T>) { return blade024<T>(blade02<T>(b0), blade24<T>(b2e)); }
template<typename T> blade024<T> combine024(zero<T>, const blade2E<T>& b2E, const blade2e<T>& b2e, const blade4<T>& b4) { return blade024<T>(blade02<T>(b2E), blade24<T>(b2e, b4)); }
template<typename T> blade024<T> combine024(const blade0<T>& b0, zero<T>, const blade2e<T>& b2e, const blade4<T>& b4) { return blade024<T>(blade02<T>(b0), blade24<T>(b2e, b4)); }
template<typename T> blade024<T> combine024(const blade0<T>& b0, const blade2E<T>& b2E, zero<T>, const blade4<T>& b4) { return blade024<T>(blade02<T>(b0, b2E), blade24<T>(b4)); }
template<typename T> blade024<T> combine024(const blade0<T>& b0, const blade2E<T>& b2E, const blade2e<T>& b2e, zero<T>) { return blade024<T>(blade02<T>(b0, b2E), blade24<T>(b2e)); }
template<typename T> blade024<T> combine024(const blade0<T>& b0, const blade2E<T>& b2E, const blade2e<T>& b2e, const blade4<T>& b4) { return blade024<T>(blade02<T>(b0, b2E), blade24<T>(b2e, b4)); }

template<typename T> zero<T> combine13(zero<T>, zero<T>) { return zero<T>(); }
template<typename T> const blade1<T>& combine13(const blade1<T>& b1, zero<T>) { return b1; }
template<typename T> const blade3<T>& combine13(zero<T>, const blade3<T>& b3) { return b3; }
template<typename T> blade13<T> combine13(blade1<T> b1, blade3<T> b3) { return blade13<T>(std::move(b1), std::move(b3)); }

template<typename B, typename T, IF_BLADE_TYPE(B, T)> const B& combine_multivector(zero<T>, const B& b) { return b; }
template<typename B, typename T, IF_NONZERO_BLADE_TYPE(B, T)> const B& combine_multivector(const B& b, zero<T>) { return b; }
template<typename BL, typename BR, typename T = get<BL>::type, IF_NONZERO_BLADE_TYPE(BL, T), IF_NONZERO_BLADE_TYPE(BR, T)>
multivector<T> combine_multivector(BL bl, BR br) { return multivector<T>(blade024<T>(std::move(bl)), blade13<T>(std::move(br))); }

// add
template<typename B, typename T, IF_BLADE_TYPE(B, T)> const B& operator + (zero<T>, const B& b) { return b; }
template<typename B, typename T, IF_NONZERO_BLADE_TYPE(B, T)> const B& operator + (const B& b, zero<T>) { return b; }

template<typename BL, typename BR, typename T = get<BL>::type, IF_NONZERO_BLADE_TYPE(BL, T), IF_NONZERO_BLADE_TYPE(BR, T)>
auto operator + (const BL& bl, const BR& br) {
    return combine_multivector(
        combine024(
            get<BL>::b0(bl) + get<BR>::b0(br),
            get<BL>::b2E(bl) + get<BR>::b2E(br),
            get<BL>::b2e(bl) + get<BR>::b2e(br),
            get<BL>::b4(bl) + get<BR>::b4(br)),
        combine13(
            get<BL>::b1(bl) + get<BR>::b1(br),
            get<BL>::b3(bl) + get<BR>::b3(br)));
}

// multiply
template<typename B, typename T, IF_NONZERO_BLADE_TYPE(B, T)> B operator * (const T& s, const B& b) { return b * s; }

template<typename T> auto operator * (const blade0<T>& a, const blade0<T>& b) { return blade0<T>(a.s * b.s); }
template<typename T> auto operator * (const blade0<T>& a, const blade1<T>& b) { return blade1<T>(a.s * b.e0, a.s * b.v); }
template<typename T> auto operator * (const blade0<T>& a, const blade2E<T>& b) { return blade2E<T>(a.s * b.biE); }
template<typename T> auto operator * (const blade0<T>& a, const blade2e<T>& b) { return blade2e<T>(a.s * b.bie); }
template<typename T> auto operator * (const blade0<T>& a, const blade3<T>& b) { return blade3<T>(a.s * b.e123, a.s * b.triP); }
template<typename T> auto operator * (const blade0<T>& a, const blade4<T>& b) { return blade4<T>(a.s * b.e0123); }
template<typename T> auto operator * (const blade1<T>& a, const blade0<T>& b) { return blade1<T>(a.e0 * b.s, a.v * b.s); }
template<typename T> auto operator * (const blade1<T>& a, const blade1<T>& b) { return blade0<T>(a.v & b.v) + blade2E<T>(a.v ^ b.v) + blade2e<T>(a.e0 * b.v - a.v * b.e0); }
template<typename T> auto operator * (const blade1<T>& a, const blade2E<T>& b) { return blade1<T>(-(a.v ^ b.biE)) + blade3<T>(a.v & b.biE, -a.e0 * b.biE); }
template<typename T> auto operator * (const blade1<T>& a, const blade2e<T>& b) { return blade1<T>(-(a.v & b.bie)) + blade3<T>(a.v ^ b.bie); }
template<typename T> auto operator * (const blade1<T>& a, const blade3<T>& b) { return blade2E<T>(a.v * b.e123) + blade2e<T>(-(a.v ^ b.triP)) + blade4<T>(a.e0 * b.e123 + (a.v & b.triP)); }
template<typename T> auto operator * (const blade1<T>& a, const blade4<T>& b) { return blade3<T>(a.v * b.e0123); }
template<typename T> auto operator * (const blade2E<T>& a, const blade0<T>& b) { return blade2E<T>(a.biE * b.s); }
template<typename T> auto operator * (const blade2E<T>& a, const blade1<T>& b) { return blade1<T>(-(a.biE ^ b.v)) + blade3<T>(a.biE & b.v, a.biE * -b.e0); }
template<typename T> auto operator * (const blade2E<T>& a, const blade2E<T>& b) { return blade0<T>(-(a.biE & b.biE)) + blade2E<T>(-(a.biE ^ b.biE)); }
template<typename T> auto operator * (const blade2E<T>& a, const blade2e<T>& b) { return blade2e<T>(-(a.biE ^ b.bie)) + blade4<T>(a.biE & b.bie); }
template<typename T> auto operator * (const blade2E<T>& a, const blade3<T>& b) { return blade1<T>(a.biE & b.triP, a.biE * -b.e123) + blade3<T>(-(a.biE ^ b.triP)); }
template<typename T> auto operator * (const blade2E<T>& a, const blade4<T>& b) { return blade2e<T>(a.biE * -b.e0123); }
template<typename T> auto operator * (const blade2e<T>& a, const blade0<T>& b) { return blade2e<T>(a.bie * b.s); }
template<typename T> auto operator * (const blade2e<T>& a, const blade1<T>& b) { return blade1<T>(a.bie & b.v) + blade3<T>(-(a.bie ^ b.v)); }
template<typename T> auto operator * (const blade2e<T>& a, const blade2E<T>& b) { return blade2e<T>(-(a.bie ^ b.biE)) + blade4<T>(a.bie & b.biE); }
template<typename T> auto operator * (const blade2e<T>& a, const blade2e<T>& b) { return zero<T>(); }
template<typename T> auto operator * (const blade2e<T>& a, const blade3<T>& b) { return blade3<T>(a.bie * -b.e123); }
template<typename T> auto operator * (const blade2e<T>& a, const blade4<T>& b) { return zero<T>(); }
template<typename T> auto operator * (const blade3<T>& a, const blade0<T>& b) { return blade3<T>(a.e123 * b.s, a.triP * b.s); }
template<typename T> auto operator * (const blade3<T>& a, const blade1<T>& b) { return blade2E<T>(a.e123 * b.v) + blade2e<T>(a.triP ^ b.v) + blade4<T>(-a.e123 * b.e0 - (a.triP & b.v)); }
template<typename T> auto operator * (const blade3<T>& a, const blade2E<T>& b) { return blade1<T>(a.triP & b.biE, -a.e123 * b.biE) + blade3<T>(-(a.triP ^ b.biE)); }
template<typename T> auto operator * (const blade3<T>& a, const blade2e<T>& b) { return blade3<T>(a.e123 * b.bie); }
template<typename T> auto operator * (const blade3<T>& a, const blade3<T>& b) { return blade0<T>(-a.e123 * b.e123) + blade2e<T>(a.triP * b.e123 - a.e123 * b.triP); }
template<typename T> auto operator * (const blade3<T>& a, const blade4<T>& b) { return blade1<T>(a.e123 * b.e0123); }
template<typename T> auto operator * (const blade4<T>& a, const blade0<T>& b) { return blade4<T>(a.e0123 * b.s); }
template<typename T> auto operator * (const blade4<T>& a, const blade1<T>& b) { return blade3<T>(-a.e0123 * b.v); }
template<typename T> auto operator * (const blade4<T>& a, const blade2E<T>& b) { return blade2e<T>(-a.e0123 * b.biE); }
template<typename T> auto operator * (const blade4<T>& a, const blade2e<T>& b) { return zero<T>(); }
template<typename T> auto operator * (const blade4<T>& a, const blade3<T>& b) { return blade1<T>(-a.e0123 * b.e123); }
template<typename T> auto operator * (const blade4<T>& a, const blade4<T>& b) { return zero<T>(); }

template<typename B, typename T, IF_BLADE_TYPE(B, T)> zero<T> operator * (zero<T> z, const B& b) { return z; }
template<typename B, typename T, IF_NONZERO_BLADE_TYPE(B, T)> zero<T> operator * (const B& b, zero<T> z) { return z; }
template<typename BL, typename BR, typename T = get<BL>::type, IF_COMPOSITE_BLADE_TYPE(BL, T), IF_PRIMITIVE_BLADE_TYPE(BR, T)>
auto operator * (const BL& lhs, const BR& rhs) { return (lhs.first() * rhs) + (lhs.second() * rhs); }
template<typename BL, typename BR, typename T = get<BL>::type, IF_NONZERO_BLADE_TYPE(BL, T), IF_COMPOSITE_BLADE_TYPE(BR, T)>
auto operator * (const BL& lhs, const BR& rhs) { return (lhs * rhs.first()) + (lhs * rhs.second()); }

// wedge
template<typename T> auto operator ^ (const blade0<T>& a, const blade0<T>& b) { return blade0<T>(a.s * b.s); }
template<typename T> auto operator ^ (const blade0<T>& a, const blade1<T>& b) { return blade1<T>(a.s * b.e0, a.s * b.v); }
template<typename T> auto operator ^ (const blade0<T>& a, const blade2E<T>& b) { return blade2E<T>(a.s * b.biE); }
template<typename T> auto operator ^ (const blade0<T>& a, const blade2e<T>& b) { return blade2e<T>(a.s * b.bie); }
template<typename T> auto operator ^ (const blade0<T>& a, const blade3<T>& b) { return blade3<T>(a.s * b.e123, a.s * b.triP); }
template<typename T> auto operator ^ (const blade0<T>& a, const blade4<T>& b) { return blade4<T>(a.s * b.e0123); }
template<typename T> auto operator ^ (const blade1<T>& a, const blade0<T>& b) { return blade1<T>(a.e0 * b.s, a.v * b.s); }
template<typename T> auto operator ^ (const blade1<T>& a, const blade1<T>& b) { return blade2E<T>(a.v ^ b.v) + blade2e<T>(a.e0 * b.v - a.v * b.e0); }
template<typename T> auto operator ^ (const blade1<T>& a, const blade2E<T>& b) { return blade3<T>(a.v & b.biE, -a.e0 * b.biE); }
template<typename T> auto operator ^ (const blade1<T>& a, const blade2e<T>& b) { return blade3<T>(a.v ^ b.bie); }
template<typename T> auto operator ^ (const blade1<T>& a, const blade3<T>& b) { return blade4<T>((a.e0 * b.e123 + (a.v & b.triP))); }
template<typename T> auto operator ^ (const blade1<T>& a, const blade4<T>& b) { return zero<T>(); }
template<typename T> auto operator ^ (const blade2E<T>& a, const blade0<T>& b) { return blade2E<T>(a.biE * b.s); }
template<typename T> auto operator ^ (const blade2E<T>& a, const blade1<T>& b) { return blade3<T>(a.biE & b.v, a.biE * -b.e0); }
template<typename T> auto operator ^ (const blade2E<T>& a, const blade2E<T>& b) { return zero<T>(); }
template<typename T> auto operator ^ (const blade2E<T>& a, const blade2e<T>& b) { return blade4<T>(a.biE & b.bie); }
template<typename T> auto operator ^ (const blade2E<T>& a, const blade3<T>& b) { return zero<T>(); }
template<typename T> auto operator ^ (const blade2E<T>& a, const blade4<T>& b) { return zero<T>(); }
template<typename T> auto operator ^ (const blade2e<T>& a, const blade0<T>& b) { return blade2e<T>(a.bie * b.s); }
template<typename T> auto operator ^ (const blade2e<T>& a, const blade1<T>& b) { return blade3<T>(-(a.bie ^ b.v)); }
template<typename T> auto operator ^ (const blade2e<T>& a, const blade2E<T>& b) { return blade4<T>(a.bie & b.biE); }
template<typename T> auto operator ^ (const blade2e<T>& a, const blade2e<T>& b) { return zero<T>(); }
template<typename T> auto operator ^ (const blade2e<T>& a, const blade3<T>& b) { return zero<T>(); }
template<typename T> auto operator ^ (const blade2e<T>& a, const blade4<T>& b) { return zero<T>(); }
template<typename T> auto operator ^ (const blade3<T>& a, const blade0<T>& b) { return blade3<T>(a.e123 * b.s, a.triP * b.s); }
template<typename T> auto operator ^ (const blade3<T>& a, const blade1<T>& b) { return blade4<T>(-a.e123 * b.e0 - (a.triP & b.v)); }
template<typename T> auto operator ^ (const blade3<T>& a, const blade2E<T>& b) { return zero<T>(); }
template<typename T> auto operator ^ (const blade3<T>& a, const blade2e<T>& b) { return zero<T>(); }
template<typename T> auto operator ^ (const blade3<T>& a, const blade3<T>& b) { return zero<T>(); }
template<typename T> auto operator ^ (const blade3<T>& a, const blade4<T>& b) { return zero<T>(); }
template<typename T> auto operator ^ (const blade4<T>& a, const blade0<T>& b) { return blade4<T>(a.e0123 * b.s); }
template<typename T> auto operator ^ (const blade4<T>& a, const blade1<T>& b) { return zero<T>(); }
template<typename T> auto operator ^ (const blade4<T>& a, const blade2E<T>& b) { return zero<T>(); }
template<typename T> auto operator ^ (const blade4<T>& a, const blade2e<T>& b) { return zero<T>(); }
template<typename T> auto operator ^ (const blade4<T>& a, const blade3<T>& b) { return zero<T>(); }
template<typename T> auto operator ^ (const blade4<T>& a, const blade4<T>& b) { return zero<T>(); }

template<typename B, typename T, IF_BLADE_TYPE(B, T)> zero<T> operator ^ (zero<T> z, const B& b) { return z; }
template<typename B, typename T, IF_NONZERO_BLADE_TYPE(B, T)> zero<T> operator ^ (const B& b, zero<T> z) { return z; }
template<typename BL, typename BR, typename T = get<BL>::type, IF_COMPOSITE_BLADE_TYPE(BL, T), IF_PRIMITIVE_BLADE_TYPE(BR, T)>
auto operator ^ (const BL& lhs, const BR& rhs) { return (lhs.first() ^ rhs) + (lhs.second() ^ rhs); }
template<typename BL, typename BR, typename T = get<BL>::type, IF_NONZERO_BLADE_TYPE(BL, T), IF_COMPOSITE_BLADE_TYPE(BR, T)>
auto operator ^ (const BL& lhs, const BR& rhs) { return (lhs ^ rhs.first()) + (lhs ^ rhs.second()); }

// dot
template<typename T> auto operator & (const blade0<T>& a, const blade0<T>& b) { return blade0<T>(a.s * b.s); }
template<typename T> auto operator & (const blade0<T>& a, const blade1<T>& b) { return blade1<T>(a.s * b.e0, a.s * b.v); }
template<typename T> auto operator & (const blade0<T>& a, const blade2E<T>& b) { return blade2E<T>(a.s * b.biE); }
template<typename T> auto operator & (const blade0<T>& a, const blade2e<T>& b) { return blade2e<T>(a.s * b.bie); }
template<typename T> auto operator & (const blade0<T>& a, const blade3<T>& b) { return blade3<T>(a.s * b.e123, a.s * b.triP); }
template<typename T> auto operator & (const blade0<T>& a, const blade4<T>& b) { return blade4<T>(a.s * b.e0123); }
template<typename T> auto operator & (const blade1<T>& a, const blade0<T>& b) { return blade1<T>(a.e0 * b.s, a.v * b.s); }
template<typename T> auto operator & (const blade1<T>& a, const blade1<T>& b) { return blade0<T>(a.v & b.v); }
template<typename T> auto operator & (const blade1<T>& a, const blade2E<T>& b) { return blade1<T>(-(a.v ^ b.biE)); }
template<typename T> auto operator & (const blade1<T>& a, const blade2e<T>& b) { return blade1<T>(-(a.v & b.bie)); }
template<typename T> auto operator & (const blade1<T>& a, const blade3<T>& b) { return blade2E<T>(a.v * b.e123) + blade2e<T>(-(a.v ^ b.triP)); }
template<typename T> auto operator & (const blade1<T>& a, const blade4<T>& b) { return blade3<T>(a.v * b.e0123); }
template<typename T> auto operator & (const blade2E<T>& a, const blade0<T>& b) { return blade2E<T>(a.biE * b.s); }
template<typename T> auto operator & (const blade2E<T>& a, const blade1<T>& b) { return blade1<T>(-(a.biE ^ b.v)); }
template<typename T> auto operator & (const blade2E<T>& a, const blade2E<T>& b) { return blade0<T>(-(a.biE & b.biE)); }
template<typename T> auto operator & (const blade2E<T>& a, const blade2e<T>& b) { return zero<T>(); }
template<typename T> auto operator & (const blade2E<T>& a, const blade3<T>& b) { return blade1<T>(a.biE & b.triP, a.biE * -b.e123); }
template<typename T> auto operator & (const blade2E<T>& a, const blade4<T>& b) { return blade2e<T>(a.biE * -b.e0123); }
template<typename T> auto operator & (const blade2e<T>& a, const blade0<T>& b) { return blade2e<T>(a.bie * b.s); }
template<typename T> auto operator & (const blade2e<T>& a, const blade1<T>& b) { return blade1<T>(a.bie & b.v); }
template<typename T> auto operator & (const blade2e<T>& a, const blade2E<T>& b) { return zero<T>(); }
template<typename T> auto operator & (const blade2e<T>& a, const blade2e<T>& b) { return zero<T>(); }
template<typename T> auto operator & (const blade2e<T>& a, const blade3<T>& b) { return zero<T>(); }
template<typename T> auto operator & (const blade2e<T>& a, const blade4<T>& b) { return zero<T>(); }
template<typename T> auto operator & (const blade3<T>& a, const blade0<T>& b) { return blade3<T>(a.e123 * b.s, a.triP * b.s); }
template<typename T> auto operator & (const blade3<T>& a, const blade1<T>& b) { return blade2E<T>(a.e123 * b.v) + blade2e<T>(a.triP ^ b.v); }
template<typename T> auto operator & (const blade3<T>& a, const blade2E<T>& b) { return blade1<T>(a.triP & b.biE, -a.e123 * b.biE); }
template<typename T> auto operator & (const blade3<T>& a, const blade2e<T>& b) { return zero<T>(); }
template<typename T> auto operator & (const blade3<T>& a, const blade3<T>& b) { return blade0<T>(-a.e123 * b.e123); }
template<typename T> auto operator & (const blade3<T>& a, const blade4<T>& b) { return blade1<T>(a.e123 * b.e0123); }
template<typename T> auto operator & (const blade4<T>& a, const blade0<T>& b) { return blade4<T>(a.e0123 * b.s); }
template<typename T> auto operator & (const blade4<T>& a, const blade1<T>& b) { return blade3<T>(-a.e0123 * b.v); }
template<typename T> auto operator & (const blade4<T>& a, const blade2E<T>& b) { return blade2e<T>(-a.e0123 * b.biE); }
template<typename T> auto operator & (const blade4<T>& a, const blade2e<T>& b) { return zero<T>(); }
template<typename T> auto operator & (const blade4<T>& a, const blade3<T>& b) { return blade1<T>(-a.e0123 * b.e123); }
template<typename T> auto operator & (const blade4<T>& a, const blade4<T>& b) { return zero<T>(); }

template<typename B, typename T, IF_BLADE_TYPE(B, T)> zero<T> operator & (zero<T> z, const B& b) { return z; }
template<typename B, typename T, IF_NONZERO_BLADE_TYPE(B, T)> zero<T> operator & (const B& b, zero<T> z) { return z; }
template<typename BL, typename BR, typename T = get<BL>::type, IF_COMPOSITE_BLADE_TYPE(BL, T), IF_PRIMITIVE_BLADE_TYPE(BR, T)>
auto operator & (const BL& lhs, const BR& rhs) { return (lhs.first() & rhs) + (lhs.second() & rhs); }
template<typename BL, typename BR, typename T = get<BL>::type, IF_NONZERO_BLADE_TYPE(BL, T), IF_COMPOSITE_BLADE_TYPE(BR, T)>
auto operator & (const BL& lhs, const BR& rhs) { return (lhs & rhs.first()) + (lhs & rhs.second()); }

// join
template<typename T> auto operator | (const blade0<T>& a, const blade0<T>& b) { return zero<T>(); }
template<typename T> auto operator | (const blade0<T>& a, const blade1<T>& b) { return zero<T>(); }
template<typename T> auto operator | (const blade0<T>& a, const blade2E<T>& b) { return zero<T>(); }
template<typename T> auto operator | (const blade0<T>& a, const blade2e<T>& b) { return zero<T>(); }
template<typename T> auto operator | (const blade0<T>& a, const blade3<T>& b) { return zero<T>(); }
template<typename T> auto operator | (const blade0<T>& a, const blade4<T>& b) { return blade0<T>(a.s * b.e0123); }
template<typename T> auto operator | (const blade1<T>& a, const blade0<T>& b) { return zero<T>(); }
template<typename T> auto operator | (const blade1<T>& a, const blade1<T>& b) { return zero<T>(); }
template<typename T> auto operator | (const blade1<T>& a, const blade2E<T>& b) { return zero<T>(); }
template<typename T> auto operator | (const blade1<T>& a, const blade2e<T>& b) { return zero<T>(); }
template<typename T> auto operator | (const blade1<T>& a, const blade3<T>& b) { return blade0<T>(-a.e0 * b.e123 - (a.v & b.triP)); }
template<typename T> auto operator | (const blade1<T>& a, const blade4<T>& b) { return blade1<T>(a.e0 * b.e0123, a.v * b.e0123); }
template<typename T> auto operator | (const blade2E<T>& a, const blade0<T>& b) { return zero<T>(); }
template<typename T> auto operator | (const blade2E<T>& a, const blade1<T>& b) { return zero<T>(); }
template<typename T> auto operator | (const blade2E<T>& a, const blade2E<T>& b) { return zero<T>(); }
template<typename T> auto operator | (const blade2E<T>& a, const blade2e<T>& b) { return blade0<T>(a.biE & b.bie); }
template<typename T> auto operator | (const blade2E<T>& a, const blade3<T>& b) { return blade1<T>(-(a.biE ^ b.triP)); }
template<typename T> auto operator | (const blade2E<T>& a, const blade4<T>& b) { return blade2E<T>(a.biE * b.e0123); }
template<typename T> auto operator | (const blade2e<T>& a, const blade0<T>& b) { return zero<T>(); }
template<typename T> auto operator | (const blade2e<T>& a, const blade1<T>& b) { return zero<T>(); }
template<typename T> auto operator | (const blade2e<T>& a, const blade2E<T>& b) { return blade0<T>(a.bie & b.biE); }
template<typename T> auto operator | (const blade2e<T>& a, const blade2e<T>& b) { return zero<T>(); }
template<typename T> auto operator | (const blade2e<T>& a, const blade3<T>& b) { return blade1<T>(a.bie & b.triP, a.bie * -b.e123); }
template<typename T> auto operator | (const blade2e<T>& a, const blade4<T>& b) { return blade2e<T>(a.bie * b.e0123); }
template<typename T> auto operator | (const blade3<T>& a, const blade0<T>& b) { return zero<T>(); }
template<typename T> auto operator | (const blade3<T>& a, const blade1<T>& b) { return blade0<T>(a.e123 * b.e0 + (a.triP & b.v)); }
template<typename T> auto operator | (const blade3<T>& a, const blade2E<T>& b) { return blade1<T>(a.triP ^ b.biE); }
template<typename T> auto operator | (const blade3<T>& a, const blade2e<T>& b) { return blade1<T>(a.triP & b.bie, -a.e123 * b.bie); }
template<typename T> auto operator | (const blade3<T>& a, const blade3<T>& b) { return blade2E<T>(a.triP * -b.e123 + a.e123 * b.triP) + blade2e<T>((a.triP ^ b.triP)); }
template<typename T> auto operator | (const blade3<T>& a, const blade4<T>& b) { return blade3<T>(a.e123 * b.e0123, a.triP * b.e0123); }
template<typename T> auto operator | (const blade4<T>& a, const blade0<T>& b) { return blade0<T>(a.e0123 * b.s); }
template<typename T> auto operator | (const blade4<T>& a, const blade1<T>& b) { return blade1<T>(a.e0123 * b.e0, a.e0123 * b.v); }
template<typename T> auto operator | (const blade4<T>& a, const blade2E<T>& b) { return blade2E<T>(a.e0123 * b.biE); }
template<typename T> auto operator | (const blade4<T>& a, const blade2e<T>& b) { return blade2e<T>(a.e0123 * b.bie); }
template<typename T> auto operator | (const blade4<T>& a, const blade3<T>& b) { return blade3<T>(a.e0123 * b.e123, a.e0123 * b.triP); }
template<typename T> auto operator | (const blade4<T>& a, const blade4<T>& b) { return blade4<T>(a.e0123 * b.e0123); }

template<typename B, typename T, IF_BLADE_TYPE(B, T)> zero<T> operator | (zero<T> z, const B& b) { return z; }
template<typename B, typename T, IF_NONZERO_BLADE_TYPE(B, T)> zero<T> operator | (const B& b, zero<T> z) { return z; }
template<typename BL, typename BR, typename T = get<BL>::type, IF_COMPOSITE_BLADE_TYPE(BL, T), IF_PRIMITIVE_BLADE_TYPE(BR, T)>
auto operator | (const BL& lhs, const BR& rhs) { return (lhs.first() | rhs) + (lhs.second() | rhs); }
template<typename BL, typename BR, typename T = get<BL>::type, IF_NONZERO_BLADE_TYPE(BL, T), IF_COMPOSITE_BLADE_TYPE(BR, T)>
auto operator | (const BL& lhs, const BR& rhs) { return (lhs | rhs.first()) + (lhs | rhs.second()); }


// TODO: inverse
// TODO: conjugate product overloads (a b ~a)

//-------------------------------------------------------------------------------

// plane with normal `n` and distance `d / norm(n)` from the origin
template<typename T>
blade1<T> plane(vector3d<T> n, T d) {
    return blade1<T>(std::move(d), std::move(n));
}

// line with direction `l` through point `P / norm(l)`
template<typename T>
blade22<T> line(const vector3d<T>& l, const vector3d<T>& P) {
    return blade22<T>(blade2E<T>(l), blade2e<T>(l ^ P));
}

// point with homogenous coordinate `h` at position `{x, y, z} / h`
template<typename T>
blade3<T> point(const vector3d<T>& P, T h) {
    return blade3<T>(std::move(h), -P);
}
// point with at position `{x, y, z}`
template<typename T>
blade3<T> point(const vector3d<T>& P) {
    return blade3<T>(identityOf(P.z), -P);
}

//-------------------------------------------------------------------------------

// TODO: projection, rejection, intersection, join
// TODO: distance, angle

} // pga
} // math
} // altruct
