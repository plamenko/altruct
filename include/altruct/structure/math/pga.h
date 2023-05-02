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
 *    square:    | 0   +1  +1  +1 |  +1   |  -1   -1   -1  |   0    0    0  |   0    |   -1     0     0     0   |
 *   ------------+----------------+-------+----------------+----------------+--------+--------------------------+
 *   blades:     |    vector      | salar |   E-bivector   |   e-bivector   | pseudo |        trivector         |
 *   c++ class:  |    blade1      | blade0|     blade2E    |     blade2e    | blade4 |          blade3          |
 *   ------------+----------------+-------+----------------+----------------+--------+--------------------------+
 *   plane       | #---#---#---#  |       |                |                |        |                          |
 *   line        |                |       |   #----#----#  |   #----#----#  |        |                          |
 *   point       |                |       |                |                |        |    #-----#-----#-----#   |
 *   ------------+----------------+-------+----------------+----------------+--------+--------------------------+
 *   rotor       |                |   #   |   #----#----#  |                |        |                          |
 *   translator  |                |   #   |                |   #----#----#  |        |                          |
 *   motor       |                |   #   |   #----#----#  |   #----#----#  |   #    |                          |
 *
 *   plane:      blade1
 *   line:       blade22 = blade2E + blade2e
 *   point:      blade3
 *
 *   rotor:      blade02E = blade0 + blade2E
 *   translator: blade02e = blade0 + balde2e
 *   motor:      blade024 = blade0 + blade2E + blade2e + blade4
 *
 *
 *   reverse:                  ~a = a.rev(); change sign of blade2 and blade3
 *   dual:                     !a = {e123, triP, e0123, bie, biE, s,     e0,   v   }
 *                                 !{e0,   v,    s,     biE, bie, e0123, e123, triP}
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
    auto norm2() const { return pl.norm2() + pr.norm2(); }                        \
    auto ninf2() const { return pl.ninf2() + pr.ninf2(); }                        \
    auto diff2() const { return pl.diff2() + pr.diff2(); }

template<typename T>
class zero {};

template<typename T>
class blade0 {
public:
    T s;

    PGA_CONSTRUCTORS_1(blade0, T, s);
    PGA_CLOSED_OPERATORS_1(blade0, T, s);

    blade0 rev() const { return blade0(s); } // same
    T norm2() const { return sqT(s); }
    T ninf2() const { return zeroOf(s); }
    T diff2() const { return norm2(); }
    auto inv() const { return rev() / norm2(); }
};

template<typename T>
class blade1 {
public:
    T e0; vector3d<T> v;

    PGA_CONSTRUCTORS_2(blade1, T, e0, zeroOf(v.z), vector3d<T>, v, make_zero_vec(e0));
    PGA_CLOSED_OPERATORS_2(blade1, T, e0, v);

    blade1 rev() const { return blade1(e0, v); } // same
    T norm2() const { return v.abs2(); }
    T ninf2() const { return sqT(e0); }
    T diff2() const { return norm2(); }
    auto inv() const { return rev() / norm2(); }
};

template<typename T>
class blade2E {
public:
    vector3d<T> biE;

    PGA_CONSTRUCTORS_1(blade2E, vector3d<T>, biE);
    PGA_CLOSED_OPERATORS_1(blade2E, T, biE);

    blade2E rev() const { return blade2E(-biE); } // negate blade2 part
    T norm2() const { return biE.abs2(); }
    T ninf2() const { return zeroOf(biE.z); }
    T diff2() const { return -norm2(); }
    auto inv() const { return rev() / norm2(); }
};

template<typename T>
class blade2e {
public:
    vector3d<T> bie;

    PGA_CONSTRUCTORS_1(blade2e, vector3d<T>, bie);
    PGA_CLOSED_OPERATORS_1(blade2e, T, bie);

    blade2e rev() const { return blade2e(-bie); } // negate blade2 part
    T norm2() const { return zeroOf(bie.z); }
    T ninf2() const { return bie.abs2(); }
    T diff2() const { return norm2(); }
    auto inv() const { return rev() / norm2(); }
};

template<typename T>
class blade3 {
public:
    T e123; vector3d<T> triP;

    PGA_CONSTRUCTORS_2(blade3, T, e123, zeroOf(triP.z), vector3d<T>, triP, make_zero_vec(e123));
    PGA_CLOSED_OPERATORS_2(blade3, T, e123, triP);

    blade3 rev() const { return blade3(-e123, -triP); } // negate blade3 part
    T norm2() const { return sqT(e123); }
    T ninf2() const { return triP.abs2(); }
    T diff2() const { return -norm2(); }
    auto inv() const { return rev() / norm2(); }
};

template<typename T>
class blade4 {
public:
    T e0123;

    PGA_CONSTRUCTORS_1(blade4, T, e0123);
    PGA_CLOSED_OPERATORS_1(blade4, T, e0123);

    blade4 rev() const { return blade4(e0123); } // same
    T norm2() const { return zeroOf(e0123); }
    T ninf2() const { return sqT(e0123); }
    T diff2() const { return norm2(); }
    auto inv() const { return rev() / norm2(); }
};

//-------------------------------------------------------------------------------

template<typename T>
class blade02E {
public:
    blade0<T> b0;
    blade2E<T> b2E;

    PGA_CONSTRUCTORS_2(blade02E, blade0<T>, b0, zeroOf(b2E.biE.z), blade2E<T>, b2E, zeroOf(b0.s));
    PGA_CLOSED_OPERATORS_2(blade02E, T, b0, b2E);
    PGA_COMPOSITE_GETTERS(blade02E, b0, b2E);
    auto inv() const { return rev() / norm2(); }
};

template<typename T>
class blade02e {
public:
    blade0<T> b0;
    blade2e<T> b2e;

    PGA_CONSTRUCTORS_2(blade02e, blade0<T>, b0, zeroOf(b2e.bie.z), blade2e<T>, b2e, zeroOf(b0.s));
    PGA_CLOSED_OPERATORS_2(blade02e, T, b0, b2e);
    PGA_COMPOSITE_GETTERS(blade02e, b0, b2e);
    auto inv() const { return rev() / norm2(); }
};

template<typename T>
class blade22 {
public:
    blade2E<T> b2E;
    blade2e<T> b2e;

    PGA_CONSTRUCTORS_2(blade22, blade2E<T>, b2E, zeroOf(b2e.bie.z), blade2e<T>, b2e, zeroOf(b2E.biE.z));
    PGA_CLOSED_OPERATORS_2(blade22, T, b2E, b2e);
    PGA_COMPOSITE_GETTERS(blade22, b2E, b2e);
    auto inv() const { return rev() / norm2(); } // only works when b2E and b2e are perpendicular! use blade024::inv() otherwise
};

template<typename T>
class blade2E4 {
public:
    blade2E<T> b2E;
    blade4<T> b4;

    PGA_CONSTRUCTORS_2(blade2E4, blade2E<T>, b2E, zeroOf(b4.e0123), blade4<T>, b4, zeroOf(b2E.biE.z));
    PGA_CLOSED_OPERATORS_2(blade2E4, T, b2E, b4);
    PGA_COMPOSITE_GETTERS(blade2E4, b2E, b4);
    auto inv() const { return rev() / norm2(); }
};

template<typename T>
class blade2e4 {
public:
    blade2e<T> b2e;
    blade4<T> b4;

    PGA_CONSTRUCTORS_2(blade2e4, blade2e<T>, b2e, zeroOf(b4.e0123), blade4<T>, b4, zeroOf(b2e.bie.z));
    PGA_CLOSED_OPERATORS_2(blade2e4, T, b2e, b4);
    PGA_COMPOSITE_GETTERS(blade2e4, b2e, b4);
    auto inv() const { return rev() / norm2(); }
};

template<typename T>
class blade024 {
public:
    blade02E<T> b02;
    blade2e4<T> b24;

    explicit blade024(const blade0<T>& b0) : b02(b0), b24(blade4<T>(zeroOf(b0.s))) {}
    explicit blade024(const blade2E<T>& b2E) : b02(b2E), b24(blade4<T>(zeroOf(b2E.biE.z))) {}
    explicit blade024(const blade2e<T>& b2e) : b02(blade0<T>(zeroOf(b2e.bie.z))), b24(b2e) {}
    explicit blade024(const blade4<T>& b4) : b02(blade0<T>(zeroOf(b4.e0123))), b24(b4) {}
    explicit blade024(const blade02e<T>& c) : b02(c.b0), b24(c.b2e) {}
    explicit blade024(const blade22<T>& c) : b02(c.b2E), b24(c.b2e) {}
    explicit blade024(const blade2E4<T>& c) : b02(c.b2E), b24(c.b4) {}
    PGA_CONSTRUCTORS_2(blade024, blade02E<T>, b02, blade0<T>(zeroOf(b24.b4.e0123)), blade2e4<T>, b24, blade4<T>(zeroOf(b02.b0.s)));
    PGA_CLOSED_OPERATORS_2(blade024, T, b02, b24);
    PGA_COMPOSITE_GETTERS(blade024, b02, b24);
    auto inv() const {
        T n2 = norm2(); T t = ((b02.b2E.biE & b24.b2e.bie) - (b02.b0.s * b24.b4.e0123)) * castOf(b02.b0.s, 2) / n2;
        return blade024(b02.rev(), b24.rev() + blade2e4<T>(blade2e<T>(b02.b2E.biE), blade4<T>(b02.b0.s)) * t) / n2;
    }
};

template<typename T>
class blade13 {
public:
    blade1<T> b1;
    blade3<T> b3;

    PGA_CONSTRUCTORS_2(blade13, blade1<T>, b1, zeroOf(b3.e123), blade3<T>, b3, zeroOf(b1.e0));
    PGA_CLOSED_OPERATORS_2(blade13, T, b1, b3);
    PGA_COMPOSITE_GETTERS(blade13, b1, b3);
    auto inv() const {
        T n2 = norm2(); T t = ((b1.v & b3.triP) + (b1.e0 * b3.e123)) * castOf(b1.e0, 2) / n2;
        return blade13(blade1<T>(b1.e0 - b3.e123 * t, b1.v), blade3<T>(-b3.e123, b1.v * t - b3.triP)) / n2;
    }
};

template<typename T>
class multivector {
public:
    blade024<T> b024;
    blade13<T> b13;

    multivector(blade0<T> b0, blade1<T> b1, blade2E<T> b2E, blade2e<T> b2e, blade3<T> b3, blade4<T> b4) :
        b13(std::move(b1), std::move(b3)), b024(blade02E<T>(std::move(b0), std::move(b2E)), blade2e4<T>(std::move(b2e), std::move(b4))) {}
    multivector(blade1<T> b1, blade0<T> b0, blade2E<T> b2E, blade2e<T> b2e, blade4<T> b4, blade3<T> b3) :
        b13(std::move(b1), std::move(b3)), b024(blade02E<T>(std::move(b0), std::move(b2E)), blade2e4<T>(std::move(b2e), std::move(b4))) {}
    multivector(blade1<T> b1, blade02E<T> b02, blade2e4<T> b24, blade3<T> b3) :
        b13(std::move(b1), std::move(b3)), b024(std::move(b02), std::move(b24)) {}
    PGA_CONSTRUCTORS_2(multivector, blade024<T>, b024, blade02E<T>(blade0<T>(zeroOf(b13.b1.e0))), blade13<T>, b13, blade1<T>(zeroOf(b024.b02.b0.s)));
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
template<typename T> struct is_composite_blade_type <blade02E<T>, T> { static const bool value = true; };
template<typename T> struct is_composite_blade_type <blade02e<T>, T> { static const bool value = true; };
template<typename T> struct is_composite_blade_type <blade22<T>, T> { static const bool value = true; };
template<typename T> struct is_composite_blade_type <blade2E4<T>, T> { static const bool value = true; };
template<typename T> struct is_composite_blade_type <blade2e4<T>, T> { static const bool value = true; };
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
template<typename T> blade2e4<T> operator ! (const blade02E<T>& c) { return blade2e4<T>(!c.b2E, !c.b0); }
template<typename T> blade2E4<T> operator ! (const blade02e<T>& c) { return blade2E4<T>(!c.b2e, !c.b0); }
template<typename T> blade22<T> operator ! (const blade22<T>& c) { return blade22<T>(!c.b2e, !c.b2E); }
template<typename T> blade02e<T> operator ! (const blade2E4<T>& c) { return blade02e<T>(!c.b4, !c.b2E); }
template<typename T> blade02E<T> operator ! (const blade2e4<T>& c) { return blade02E<T>(!c.b4, !c.b2e); }
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
template<typename T> struct get<blade02E<T>> {
    using type = T;
    static const blade0<T>& b0(const blade02E<T>& c) { return c.b0; }
    static zero<T> b1(const blade02E<T>& c) { return zero<T>(); }
    static const blade2E<T>& b2E(const blade02E<T>& c) { return c.b2E; }
    static zero<T> b2e(const blade02E<T>& c) { return zero<T>(); }
    static zero<T> b3(const blade02E<T>& c) { return zero<T>(); }
    static zero<T> b4(const blade02E<T>& c) { return zero<T>(); }
};
template<typename T> struct get<blade02e<T>> {
    using type = T;
    static const blade0<T>& b0(const blade02e<T>& c) { return c.b0; }
    static zero<T> b1(const blade02e<T>& c) { return zero<T>(); }
    static zero<T> b2E(const blade02e<T>& c) { return zero<T>(); }
    static const blade2e<T>& b2e(const blade02e<T>& c) { return c.b2e; }
    static zero<T> b3(const blade02e<T>& c) { return zero<T>(); }
    static zero<T> b4(const blade02e<T>& c) { return zero<T>(); }
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
template<typename T> struct get<blade2E4<T>> {
    using type = T;
    static zero<T> b0(const blade2E4<T>& c) { return zero<T>(); }
    static zero<T> b1(const blade2E4<T>& c) { return zero<T>(); }
    static const blade2E<T>& b2E(const blade2E4<T>& c) { return c.b2E; }
    static zero<T> b2e(const blade2E4<T>& c) { return zero<T>(); }
    static zero<T> b3(const blade2E4<T>& c) { return zero<T>(); }
    static const blade4<T>& b4(const blade2E4<T>& c) { return c.b4; }
};
template<typename T> struct get<blade2e4<T>> {
    using type = T;
    static zero<T> b0(const blade2e4<T>& c) { return zero<T>(); }
    static zero<T> b1(const blade2e4<T>& c) { return zero<T>(); }
    static zero<T> b2E(const blade2e4<T>& c) { return zero<T>(); }
    static const blade2e<T>& b2e(const blade2e4<T>& c) { return c.b2e; }
    static zero<T> b3(const blade2e4<T>& c) { return zero<T>(); }
    static const blade4<T>& b4(const blade2e4<T>& c) { return c.b4; }
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
template<typename T> blade02E<T> combine024(const blade0<T>& b0, const blade2E<T>& b2E, zero<T>, zero<T>) { return blade02E<T>(b0, b2E); }
template<typename T> blade02e<T> combine024(const blade0<T>& b0, zero<T>, const blade2e<T>& b2e, zero<T>) { return blade02e<T>(b0, b2e); }
template<typename T> blade22<T> combine024(zero<T>, const blade2E<T>& b2E, const blade2e<T>& b2e, zero<T>) { return blade22<T>(b2E, b2e); }
template<typename T> blade2E4<T> combine024(zero<T>, const blade2E<T>& b2E, zero<T>, const blade4<T>& b4) { return blade2E4<T>(b2E, b4); }
template<typename T> blade2e4<T> combine024(zero<T>, zero<T>, const blade2e<T>& b2e, const blade4<T>& b4) { return blade2e4<T>(b2e, b4); }
template<typename T> blade024<T> combine024(const blade0<T>& b0, zero<T>, zero<T>, const blade4<T>& b4) { return blade024<T>(blade02E<T>(b0), blade2e4<T>(b4)); }
template<typename T> blade024<T> combine024(zero<T>, const blade2E<T>& b2E, const blade2e<T>& b2e, const blade4<T>& b4) { return blade024<T>(blade02E<T>(b2E), blade2e4<T>(b2e, b4)); }
template<typename T> blade024<T> combine024(const blade0<T>& b0, zero<T>, const blade2e<T>& b2e, const blade4<T>& b4) { return blade024<T>(blade02E<T>(b0), blade2e4<T>(b2e, b4)); }
template<typename T> blade024<T> combine024(const blade0<T>& b0, const blade2E<T>& b2E, zero<T>, const blade4<T>& b4) { return blade024<T>(blade02E<T>(b0, b2E), blade2e4<T>(b4)); }
template<typename T> blade024<T> combine024(const blade0<T>& b0, const blade2E<T>& b2E, const blade2e<T>& b2e, zero<T>) { return blade024<T>(blade02E<T>(b0, b2E), blade2e4<T>(b2e)); }
template<typename T> blade024<T> combine024(const blade0<T>& b0, const blade2E<T>& b2E, const blade2e<T>& b2e, const blade4<T>& b4) { return blade024<T>(blade02E<T>(b0, b2E), blade2e4<T>(b2e, b4)); }

template<typename T> zero<T> combine13(zero<T>, zero<T>) { return zero<T>(); }
template<typename T> const blade1<T>& combine13(const blade1<T>& b1, zero<T>) { return b1; }
template<typename T> const blade3<T>& combine13(zero<T>, const blade3<T>& b3) { return b3; }
template<typename T> blade13<T> combine13(blade1<T> b1, blade3<T> b3) { return blade13<T>(std::move(b1), std::move(b3)); }

template<typename B, typename T, IF_BLADE_TYPE(B, T)> const B& combine_multivector(zero<T>, const B& b) { return b; }
template<typename B, typename T, IF_NONZERO_BLADE_TYPE(B, T)> const B& combine_multivector(const B& b, zero<T>) { return b; }
template<typename BL, typename BR, typename T = typename get<BL>::type, IF_NONZERO_BLADE_TYPE(BL, T), IF_NONZERO_BLADE_TYPE(BR, T)>
multivector<T> combine_multivector(BL bl, BR br) { return multivector<T>(blade024<T>(std::move(bl)), blade13<T>(std::move(br))); }

// add
template<typename B, typename T, IF_BLADE_TYPE(B, T)> const B& operator + (const B& b, zero<T>) { return b; }
template<typename B, typename T, IF_NONZERO_BLADE_TYPE(B, T)> const B& operator + (zero<T>, const B& b) { return b; }

template<typename BL, typename BR, typename T = typename get<BL>::type, IF_NONZERO_BLADE_TYPE(BL, T), IF_NONZERO_BLADE_TYPE(BR, T)>
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

// subtract
template<typename B, typename T, IF_BLADE_TYPE(B, T)> const B& operator - (const B& b, zero<T>) { return b; }
template<typename B, typename T, IF_NONZERO_BLADE_TYPE(B, T)> const B& operator - (zero<T>, const B& b) { return -b; }

template<typename BL, typename BR, typename T = typename get<BL>::type, IF_NONZERO_BLADE_TYPE(BL, T), IF_NONZERO_BLADE_TYPE(BR, T)>
auto operator - (const BL& bl, const BR& br) {
    return combine_multivector(
        combine024(
            get<BL>::b0(bl) - get<BR>::b0(br),
            get<BL>::b2E(bl) - get<BR>::b2E(br),
            get<BL>::b2e(bl) - get<BR>::b2e(br),
            get<BL>::b4(bl) - get<BR>::b4(br)),
        combine13(
            get<BL>::b1(bl) - get<BR>::b1(br),
            get<BL>::b3(bl) - get<BR>::b3(br)));
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
template<typename BL, typename BR, typename T = typename get<BL>::type, IF_COMPOSITE_BLADE_TYPE(BL, T), IF_PRIMITIVE_BLADE_TYPE(BR, T)>
auto operator * (const BL& lhs, const BR& rhs) { return (lhs.first() * rhs) + (lhs.second() * rhs); }
template<typename BL, typename BR, typename T = typename get<BL>::type, IF_NONZERO_BLADE_TYPE(BL, T), IF_COMPOSITE_BLADE_TYPE(BR, T)>
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
template<typename BL, typename BR, typename T = typename get<BL>::type, IF_COMPOSITE_BLADE_TYPE(BL, T), IF_PRIMITIVE_BLADE_TYPE(BR, T)>
auto operator ^ (const BL& lhs, const BR& rhs) { return (lhs.first() ^ rhs) + (lhs.second() ^ rhs); }
template<typename BL, typename BR, typename T = typename get<BL>::type, IF_NONZERO_BLADE_TYPE(BL, T), IF_COMPOSITE_BLADE_TYPE(BR, T)>
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
template<typename BL, typename BR, typename T = typename get<BL>::type, IF_COMPOSITE_BLADE_TYPE(BL, T), IF_PRIMITIVE_BLADE_TYPE(BR, T)>
auto operator & (const BL& lhs, const BR& rhs) { return (lhs.first() & rhs) + (lhs.second() & rhs); }
template<typename BL, typename BR, typename T = typename get<BL>::type, IF_NONZERO_BLADE_TYPE(BL, T), IF_COMPOSITE_BLADE_TYPE(BR, T)>
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
template<typename BL, typename BR, typename T = typename get<BL>::type, IF_COMPOSITE_BLADE_TYPE(BL, T), IF_PRIMITIVE_BLADE_TYPE(BR, T)>
auto operator | (const BL& lhs, const BR& rhs) { return (lhs.first() | rhs) + (lhs.second() | rhs); }
template<typename BL, typename BR, typename T = typename get<BL>::type, IF_NONZERO_BLADE_TYPE(BL, T), IF_COMPOSITE_BLADE_TYPE(BR, T)>
auto operator | (const BL& lhs, const BR& rhs) { return (lhs | rhs.first()) + (lhs | rhs.second()); }

//-------------------------------------------------------------------------------

// sandwich product: a % b = (-1)^m b a ~b for orthogonal b; a % b = (-1)^m b a b^-1 for orthonormal b; m = number of reflections
//   blade0    0 reflections = nop
// > blade1    1 reflection
//   blade2E   ?
//   blade2e   ?
//   blade3    ?
//   blade4    ?
//   blade02E  2 reflections (non-parallel planes through origin) = rotation at origin
//   blade02e  2 reflections (parallel planes) = translation
//   blade22   2 reflections (perpendicular planes) = 180deg rotation around intersection
//   blade2E4  ?
//   blade2e4  ?
//   blade024  4 reflections = generic roto-translation
//   blade13   3 reflections

#define PGA_DEF_024(b024, c2)                                               \
    const auto& b024##s = b024.b02.b0.s, b024##0123 = b024.b24.b4.e0123;    \
    const auto& b024##biE = b024.b02.b2E.biE, b024##bie = b024.b24.b2e.bie; \
    const auto c2 = castOf(b024##s, 2);

#define PGA_DEF_13(b13, c2)                                                 \
    const auto& b13##0 = b13.b1.e0, b13##123 = b13.b3.e123;                 \
    const auto& b13##v = b13.b1.v, b13##triP = b13.b3.triP;                 \
    const auto c2 = castOf(b13##0, 2);

#define PGA_DOT_CROSS(a, b, a_o_b, a_x_b)                                   \
    const auto a_x_b = (a ^ b); const auto a_o_b = (a & b);

template<typename T> auto operator % (const blade0<T>& a, const blade0<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade0<T>& a, const blade1<T>& b) { return a * -b.norm2(); }
template<typename T> auto operator % (const blade0<T>& a, const blade2E<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade0<T>& a, const blade2e<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade0<T>& a, const blade3<T>& b) { return a * -b.norm2(); }
template<typename T> auto operator % (const blade0<T>& a, const blade4<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade0<T>& a, const blade02E<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade0<T>& a, const blade02e<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade0<T>& a, const blade22<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade0<T>& a, const blade2E4<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade0<T>& a, const blade2e4<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade0<T>& a, const blade024<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade0<T>& a, const blade13<T>& b) { return a * -b.norm2(); }

template<typename T> auto operator % (const blade1<T>& a, const blade0<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade1<T>& a, const blade1<T>& b) { return a * b.norm2() - b * (a.v & b.v) * castOf(a.e0, 2); }
template<typename T> auto operator % (const blade1<T>& a, const blade2E<T>& b) { return blade1<T>(a.e0 * b.norm2(), a.v * b.diff2() + b.biE * ((a.v & b.biE) * castOf(a.e0, 2))); }
template<typename T> auto operator % (const blade1<T>& a, const blade2e<T>& b) { return blade1<T>(zeroOf(a.e0)); }
template<typename T> auto operator % (const blade1<T>& a, const blade3<T>& b) { return blade1<T>(a.e0 * b.norm2() + b.e123 * (a.v & b.triP) * castOf(a.e0, 2), a.v * b.diff2()); }
template<typename T> auto operator % (const blade1<T>& a, const blade4<T>& b) { return blade1<T>(zeroOf(a.e0)); }
template<typename T> auto operator % (const blade1<T>& a, const blade02E<T>& b) {
    const auto c2 = castOf(a.e0, 2); PGA_DOT_CROSS(a.v, b.b2E.biE, av_o_bbiE, av_x_bbiE);
    return blade1<T>(a.e0 * b.norm2(), a.v * b.diff2() + (b.b2E.biE * av_o_bbiE + b.b0.s * av_x_bbiE) * c2);
}
template<typename T> auto operator % (const blade1<T>& a, const blade02e<T>& b) {
    return blade1<T>(a.e0 * b.norm2() + b.b0.s * (a.v & b.b2e.bie) * castOf(a.e0, 2), a.v * b.diff2());
}
template<typename T> auto operator % (const blade1<T>& a, const blade22<T>& b) {
    const auto c2 = castOf(a.e0, 2); PGA_DOT_CROSS(a.v, b.b2E.biE, av_o_bbiE, av_x_bbiE);
    return blade1<T>(a.e0 * b.norm2() + (b.b2e.bie & av_x_bbiE) * c2, a.v * b.diff2() + b.b2E.biE * (av_o_bbiE * c2));
}
template<typename T> auto operator % (const blade1<T>& a, const blade2E4<T>& b) {
    const auto av_o_bbiE_c2 = (a.v & b.b2E.biE) * castOf(a.e0, 2);
    return blade1<T>(a.e0 * b.norm2() + b.b4.e0123 * av_o_bbiE_c2, a.v * b.diff2() + b.b2E.biE * av_o_bbiE_c2);
}
template<typename T> auto operator % (const blade1<T>& a, const blade2e4<T>& b) {
    return blade1<T>(zeroOf(a.e0));
}
template<typename T> auto operator % (const blade1<T>& a, const blade024<T>& b) {
    PGA_DEF_024(b, c2); PGA_DOT_CROSS(a.v, bbiE, av_o_bbiE, av_x_bbiE);
    return blade1<T>(a.e0 * b.norm2() + (bs * (a.v & bbie) + b0123 * av_o_bbiE + (bbie & av_x_bbiE)) * c2, a.v * b.diff2() + (bbiE * av_o_bbiE + bs * av_x_bbiE) * c2);
}
template<typename T> auto operator % (const blade1<T>& a, const blade13<T>& b) {
    PGA_DEF_13(b, c2); PGA_DOT_CROSS(a.v, bv, av_o_bv, av_x_bv);
    return blade1<T>(a.e0 * b.norm2() + (-b0 * av_o_bv + b123 * (a.v & btriP) - (btriP & av_x_bv)) * c2, a.v * b.diff2() + (b123 * av_x_bv - bv * av_o_bv) * c2);
}

template<typename T> auto operator % (const blade2E<T>& a, const blade0<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade2E<T>& a, const blade1<T>& b) { return blade2E<T>(a.biE * b.diff2() - b.v * ((a.biE & b.v) * castOf(a.biE.z, 2))) + blade2e<T>(b.e0 * castOf(a.biE.z, 2) * (a.biE ^ b.v)); }
template<typename T> auto operator % (const blade2E<T>& a, const blade2E<T>& b) { return blade2E<T>(a.biE * b.diff2() + b.biE * ((a.biE & b.biE) * castOf(a.biE.z, 2))); }
template<typename T> auto operator % (const blade2E<T>& a, const blade2e<T>& b) { return blade2E<T>(zeroOf(a.biE.z)); }
template<typename T> auto operator % (const blade2E<T>& a, const blade3<T>& b) { return blade2E<T>(a.biE * b.diff2()) + blade2e<T>(b.e123 * castOf(a.biE.z, 2) * (a.biE ^ b.triP)); }
template<typename T> auto operator % (const blade2E<T>& a, const blade4<T>& b) { return blade2E<T>(zeroOf(a.biE.z)); }
template<typename T> auto operator % (const blade2E<T>& a, const blade02E<T>& b) {
    const auto c2 = castOf(a.biE.z, 2); PGA_DOT_CROSS(a.biE, b.b2E.biE, abiE_o_bbiE, abiE_x_bbiE);
    return blade2E<T>(a.biE * b.diff2() + (b.b2E.biE * abiE_o_bbiE + b.b0.s * abiE_x_bbiE) * c2);
}
template<typename T> auto operator % (const blade2E<T>& a, const blade02e<T>& b) {
    return blade2E<T>(a.biE * b.diff2()) + blade2e<T>(b.b0.s * castOf(a.biE.z, 2) * (a.biE ^ b.b2e.bie));
}
template<typename T> auto operator % (const blade2E<T>& a, const blade22<T>& b) {
    const auto abiE_o_bbiE_c2 = (a.biE & b.b2E.biE) * castOf(a.biE.z, 2);
    const auto abiE_o_bbie_c2 = (a.biE & b.b2e.bie) * castOf(a.biE.z, 2);
    return blade2E<T>(a.biE * b.diff2() + b.b2E.biE * abiE_o_bbiE_c2) + blade2e<T>(b.b2E.biE * abiE_o_bbie_c2 + b.b2e.bie * abiE_o_bbiE_c2);
}
template<typename T> auto operator % (const blade2E<T>& a, const blade2E4<T>& b) {
    const auto c2 = castOf(a.biE.z, 2); PGA_DOT_CROSS(a.biE, b.b2E.biE, abiE_o_bbiE, abiE_x_bbiE);
    return blade2E<T>(a.biE * b.diff2() + b.b2E.biE * (abiE_o_bbiE * c2)) + blade2e<T>(-c2 * b.b4.e0123 * abiE_x_bbiE);
}
template<typename T> auto operator % (const blade2E<T>& a, const blade2e4<T>& b) {
    return blade2E<T>(zeroOf(a.biE.z));
}
template<typename T> auto operator % (const blade2E<T>& a, const blade024<T>& b) {
    PGA_DEF_024(b, c2); PGA_DOT_CROSS(a.biE, bbiE, abiE_o_bbiE, abiE_x_bbiE); PGA_DOT_CROSS(a.biE, bbie, abiE_o_bbie, abiE_x_bbie);
    return blade2E<T>(a.biE * b.diff2() + (bbiE * abiE_o_bbiE + bs * abiE_x_bbiE) * c2) +
        blade2e<T>((a.biE * (b0123 * bs * -c2) + bbiE * abiE_o_bbie + bs * abiE_x_bbie + bbie * abiE_o_bbiE - b0123 * abiE_x_bbiE) * c2);
}
template<typename T> auto operator % (const blade2E<T>& a, const blade13<T>& b) {
    PGA_DEF_13(b, c2); PGA_DOT_CROSS(a.biE, bv, abiE_o_bv, abiE_x_bv); PGA_DOT_CROSS(a.biE, btriP, abiE_o_btriP, abiE_x_btriP);
    return blade2E<T>(a.biE * b.diff2() + (b123 * abiE_x_bv - bv * abiE_o_bv) * c2) +
        blade2e<T>((a.biE * (b0 * b123 * -c2) - btriP * abiE_o_bv + b0 * abiE_x_bv - bv * abiE_o_btriP + b123 * abiE_x_btriP) * c2);
}

template<typename T> auto operator % (const blade2e<T>& a, const blade0<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade2e<T>& a, const blade1<T>& b) { return blade2e<T>(a.bie * -b.diff2() + b.v * (a.bie & b.v) * castOf(a.bie.z, 2)); }
template<typename T> auto operator % (const blade2e<T>& a, const blade2E<T>& b) { return blade2e<T>(a.bie * b.diff2() + b.biE * ((a.bie & b.biE) * castOf(a.bie.z, 2))); }
template<typename T> auto operator % (const blade2e<T>& a, const blade2e<T>& b) { return blade2e<T>(zeroOf(a.bie.z)); }
template<typename T> auto operator % (const blade2e<T>& a, const blade3<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade2e<T>& a, const blade4<T>& b) { return blade2e<T>(zeroOf(a.bie.z)); }
template<typename T> auto operator % (const blade2e<T>& a, const blade02E<T>& b) {
    const auto c2 = castOf(a.bie.z, 2); PGA_DOT_CROSS(a.bie, b.b2E.biE, abie_o_bbiE, abie_x_bbiE);
    return blade2e<T>(a.bie * b.diff2() + (b.b0.s * abie_x_bbiE + b.b2E.biE * abie_o_bbiE) * c2);
}
template<typename T> auto operator % (const blade2e<T>& a, const blade02e<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade2e<T>& a, const blade22<T>& b) { return a % b.b2E; }
template<typename T> auto operator % (const blade2e<T>& a, const blade2E4<T>& b) { return a % b.b2E; }
template<typename T> auto operator % (const blade2e<T>& a, const blade2e4<T>& b) { return blade2e<T>(zeroOf(a.bie.z)); }
template<typename T> auto operator % (const blade2e<T>& a, const blade024<T>& b) {
    PGA_DEF_024(b, c2); PGA_DOT_CROSS(a.bie, bbiE, abie_o_bbiE, abie_x_bbiE);
    return blade2e<T>(a.bie * b.diff2() + (bs * abie_x_bbiE + bbiE * abie_o_bbiE) * c2);
}
template<typename T> auto operator % (const blade2e<T>& a, const blade13<T>& b) {
    PGA_DEF_13(b, c2); PGA_DOT_CROSS(a.bie, bv, abie_o_bv, abie_x_bv);
    return blade2e<T>(a.bie * -b.diff2() + (bv * abie_o_bv - b123 * abie_x_bv) * c2);
}

template<typename T> auto operator % (const blade3<T>& a, const blade0<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade3<T>& a, const blade1<T>& b) { return blade3<T>(a.e123 * -b.norm2(), a.triP * -b.diff2() + b.v * (((a.triP & b.v) + a.e123 * b.e0) * castOf(a.e123, 2))); }
template<typename T> auto operator % (const blade3<T>& a, const blade2E<T>& b) { return blade3<T>(a.e123 * b.norm2(), a.triP * b.diff2() + b.biE * ((a.triP & b.biE) * castOf(a.e123, 2))); }
template<typename T> auto operator % (const blade3<T>& a, const blade2e<T>& b) { return blade3<T>(zeroOf(a.e123)); }
template<typename T> auto operator % (const blade3<T>& a, const blade3<T>& b) { return blade3<T>(a.e123 * -b.norm2(), a.triP * -b.diff2() - b.triP * (a.e123 * b.e123 * castOf(a.e123, 2))); }
template<typename T> auto operator % (const blade3<T>& a, const blade4<T>& b) { return blade3<T>(zeroOf(a.e123)); }
template<typename T> auto operator % (const blade3<T>& a, const blade02E<T>& b) {
    const auto c2 = castOf(a.e123, 2); const auto& bs = b.b0.s; const auto& bbiE = b.b2E.biE;
    return blade3<T>(a.e123 * b.norm2(), a.triP * b.diff2() + (bbiE * (a.triP & bbiE) + bs * (a.triP ^ bbiE)) * c2);
}
template<typename T> auto operator % (const blade3<T>& a, const blade02e<T>& b) {
    const auto c2 = castOf(a.e123, 2); const auto& bs = b.b0.s; const auto& bbie = b.b2e.bie;
    return blade3<T>(a.e123 * b.norm2(), a.triP * b.diff2() - bbie * (a.e123 * bs * c2));
}
template<typename T> auto operator % (const blade3<T>& a, const blade22<T>& b) {
    const auto c2 = castOf(a.e123, 2); const auto& bbiE = b.b2E.biE; const auto& bbie = b.b2e.bie;
    return blade3<T>(a.e123 * b.norm2(), a.triP * b.diff2() + (bbiE * (a.triP & bbiE) + a.e123 * (bbiE ^ bbie)) * c2);
}
template<typename T> auto operator % (const blade3<T>& a, const blade2E4<T>& b) {
    const auto c2 = castOf(a.e123, 2); const auto& bbiE = b.b2E.biE; const auto& b0123 = b.b4.e0123;
    return blade3<T>(a.e123 * b.norm2(), a.triP * b.diff2() + bbiE * (((a.triP & bbiE) - a.e123 * b0123) * c2));
}
template<typename T> auto operator % (const blade3<T>& a, const blade2e4<T>& b) {
    return blade3<T>(zeroOf(a.e123));
}
template<typename T> auto operator % (const blade3<T>& a, const blade024<T>& b) {
    PGA_DEF_024(b, c2); PGA_DOT_CROSS(a.triP, bbiE, atriP_o_bbiE, atriP_x_bbiE);
    return blade3<T>(a.e123 * b.norm2(), a.triP * b.diff2() + (bbiE * (atriP_o_bbiE - a.e123 * b0123) - bbie * (a.e123 * bs) + bs * atriP_x_bbiE + a.e123 * (bbiE ^ bbie)) * c2);
}
template<typename T> auto operator % (const blade3<T>& a, const blade13<T>& b) {
    PGA_DEF_13(b, c2); PGA_DOT_CROSS(a.triP, bv, atriP_o_bv, atriP_x_bv);
    return blade3<T>(a.e123 * -b.norm2(), a.triP * -b.diff2() + (bv * (atriP_o_bv + a.e123 * b0) - btriP * (a.e123 * b123) - b123 * atriP_x_bv - a.e123 * (bv ^ btriP)) * c2);
}

template<typename T> auto operator % (const blade4<T>& a, const blade0<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade4<T>& a, const blade1<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade4<T>& a, const blade2E<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade4<T>& a, const blade2e<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade4<T>& a, const blade3<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade4<T>& a, const blade4<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade4<T>& a, const blade02E<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade4<T>& a, const blade02e<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade4<T>& a, const blade22<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade4<T>& a, const blade2E4<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade4<T>& a, const blade2e4<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade4<T>& a, const blade024<T>& b) { return a * b.norm2(); }
template<typename T> auto operator % (const blade4<T>& a, const blade13<T>& b) { return a * b.norm2(); }

// composite a can be distribued, but composite b can not! needs explicit overloads for all types of b
template<typename BL, typename BR, typename T = typename get<BL>::type, IF_COMPOSITE_BLADE_TYPE(BL, T), IF_NONZERO_BLADE_TYPE(BR, T)>
auto operator % (const BL& lhs, const BR& rhs) { return (lhs.first() % rhs) + (lhs.second() % rhs); }

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

// translator by vector `t`
template<typename T>
blade02e<T> translator(const vector3d<T>& t) {
    auto id = identityOf(t.z);
    auto ht = t / castOf(t.z, 2);
    return blade02e<T>(blade0<T>(id), blade2e<T>(ht));
}

// rotor along axis `n` by angle `a` given `cos(a/2)` and `sin(a/2)`
template<typename T>
blade02E<T> rotor(const vector3d<T>& n, T cos_a2, T sin_a2) {
    return blade02E<T>(blade0<T>(cos_a2), blade2E<T>(n.unit() * -sin_a2));
}
// rotor along axis `n` by angle `a` in radians
template<typename T>
blade02E<T> rotor(const vector3d<T>& n, T a) {
    auto ha = a / castOf(a, 2);
    return rotor(n, cos(ha), sin(ha));
}

//-------------------------------------------------------------------------------

// TODO: projection, rejection
// TODO: distance, angle

} // pga
} // math
} // altruct
