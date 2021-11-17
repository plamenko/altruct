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
 *   reverse:                  a.rev(); change sign of blade2 and blade3
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
class blade1 {
public:
    T e0; vector3d<T> v;

    blade1() {}
    blade1(T e0) : e0(e0), v({ zeroOf(e0), zeroOf(e0) , zeroOf(e0) }) {}
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

    blade1 rev() const { return blade1(e0, v); } // same
};

template<typename T>
class blade02 {
public:
    T s; vector3d<T> biE;

    blade02() {}
    blade02(T s, vector3d<T> biE) : s(std::move(s)), biE(std::move(biE)) {}

    blade02 rev() const { return blade02(s, -biE); } // negate blade2 part
};

template<typename T>
class blade24 {
public:
    vector3d<T> bie; T e0123;

    blade24() {}
    blade24(vector3d<T> bie, T e0123) : bie(std::move(bie)), e0123(std::move(e0123)) {}

    blade24 rev() const { return blade24(-bie, e0123); } // negate blade2 part
};

template<typename T>
class blade3 {
public:
    T e123; vector3d<T> triP;

    blade3() {}
    blade3(T e123, vector3d<T> triP) : e123(std::move(e123)), triP(std::move(triP)) {}

    blade3 rev() const { return blade3(-e123, -triP); } // negate blade3 part
};

template<typename T>
class blade13 {
public:
    blade1<T> b1;
    blade3<T> b3;

    blade13() {}
    blade13(blade1<T> b1, blade3<T> b3) : b1(std::move(b1)), b3(std::move(b3)) {}

    blade13 rev() const { return blade13(b1.rev(), b3.rev()); }
};

template<typename T>
class blade024 {
public:
    blade02<T> b02;
    blade24<T> b24;

    blade024() {}
    blade024(blade02<T> b02, blade24<T> b24) : b02(std::move(b02)), b24(std::move(b24)) {}

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

    multivector rev() const { return multivector(b13.rev(), b024.rev()); }
};

//-------------------------------------------------------------------------------

// dual
template<typename T> blade3<T> operator ! (const blade1<T>& b1) { return blade3<T>(b1.e0, b1.v); }
template<typename T> blade24<T> operator ! (const blade02<T>& b02) { return blade24<T>(b02.biE, b02.s); }
template<typename T> blade02<T> operator ! (const blade24<T>& b24) { return blade02<T>(b24.e0123, b24.bie); }
template<typename T> blade1<T> operator ! (const blade3<T>& b3) { return blade1<T>(b3.e123, b3.triP); }
template<typename T> blade13<T> operator ! (const blade13<T>& b) { return blade13<T>(!b.b3, !b.b1); }
template<typename T> blade024<T> operator ! (const blade024<T>& b) { return blade024<T>(!b.b24, !b.b02); }
template<typename T> multivector<T> operator ! (const multivector<T>& m) { return multivector<T>(!m.b13, !m.b024); }

// multiply
template<typename T> blade024<T> operator * (const blade1<T>& a, const blade1<T>& b) {
    return blade024<T>(blade02<T>(a.v & b.v, a.v ^ b.v), blade24<T>(a.e0 * b.v - b.e0 * a.v, zeroOf(a.e0)));
}
template<typename T> blade13<T> operator * (const blade1<T>& a, const blade02<T>& b) {
    return blade13<T>(blade1<T>(a.e0 * b.s, a.v * b.s - (a.v ^ b.biE)), blade3<T>(a.v & b.biE, -a.e0 * b.biE));
}
template<typename T> blade13<T> operator * (const blade1<T>& a, const blade24<T>& b) {
    return blade13<T>(blade1<T>(-a.v & b.bie), blade3<T>(zeroOf(b.e0123), a.v * b.e0123 + (a.v ^ b.bie)));
}
template<typename T> blade024<T> operator * (const blade1<T>& a, const blade3<T>& b) {
    return blade024<T>(blade02<T>(zeroOf(a.e0), a.v * b.e123), blade24<T>(-a.v ^ b.triP, a.e0 * b.e123 + (a.v & b.triP)));
}
// TODO: multiply overloads

// TODO: wedge overloads
// TODO: dot overloads
// TODO: join overloads

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
    return blade024<T>(blade02<T>(zero, l), blade24<T>(cross(l, P), zero));
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
