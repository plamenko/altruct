#pragma once

#include "altruct/algorithm/math/base.h"

namespace altruct {
namespace math {

/**
 * Clifford Algebra Cl(3, 0) a.k.a Geometric Algebra G(3, 0) in 3D
 *   orthogonal basis: {e1, e2, e3}
 *   elements: {1, e1, e2, e3, e2e3, e3e1, e1e2, e1e2e3}
 *   even part (rotor): {1, e2e3, e3e1, e1e2} 
 *   odd part (vector): {e1, e2, e3, e1e2e3}
 *   representation: s + x e1 + y e2 + z e3 + yz e2e3 + zx e3e1 + xy e1e2 + w e1e2e3
 **/
namespace cl3 {

//------------------------------------------------------------------------------------------------/
// Even subalgebra of G(3,0), isomorphic to quaternion
template<typename T>
class rotor {
public:
    T s;          // scalar
    T yz, zx, xy; // bivectors

    rotor() : s(1), yz(0), zx(0), xy(0) {}
    rotor(T s) : s(std::move(s)), yz(zeroOf(s)), zx(zeroOf(s)), xy(zeroOf(s)) {}
    rotor(T s, T yz, T zx, T xy) : s(std::move(s)), yz(std::move(yz)), zx(std::move(zx)), xy(std::move(xy)) {}

    bool operator == (const rotor& r) const { return (s == r.s) && (yz == r.yz) && (zx == r.zx) && (xy == r.xy); }
    bool operator != (const rotor& r) const { return !(*this == r); }
    bool operator <  (const rotor& r) const { 
        if (s != r.s) return s < r.s;
        if (yz != r.yz) return yz < r.yz;
        if (zx != r.zx) return zx < r.zx;
                        return xy < r.xy;
    }
    bool operator >  (const rotor& r) const { return (r < *this); }
    bool operator <= (const rotor& r) const { return !(r < *this); }
    bool operator >= (const rotor& r) const { return !(*this < r); }

    rotor& operator += (const rotor& r) { s += r.s; yz += r.yz; zx += r.zx; xy += r.xy; return *this; }
    rotor& operator -= (const rotor& r) { s -= r.s; yz -= r.yz; zx -= r.zx; xy -= r.xy; return *this; }
    rotor& operator *= (const rotor& r) { return *this = *this * r; }
    rotor& operator /= (const rotor& r) { return *this = *this / r; }

    rotor& operator *= (const T& t) { s *= t; yz *= t; zx *= t; xy *= t; return *this; }
    rotor& operator /= (const T& t) { s /= t; yz /= t; zx /= t; xy /= t; return *this; }

    rotor  operator +  (const rotor& r) const { return rotor(s + r.s, yz + r.yz, zx + r.zx, xy + r.xy); }
    rotor  operator -  (const rotor& r) const { return rotor(s - r.s, yz - r.yz, zx - r.zx, xy - r.xy); }
    rotor  operator -  ()               const { return rotor(-s, -yz, -zx, -xy); }
    rotor  operator *  (const rotor& r) const {
        return rotor(s * r.s - yz * r.yz - zx * r.zx - xy * r.xy,
            s * r.yz + yz * r.s - zx * r.xy + xy * r.zx,
            s * r.zx + yz * r.xy + zx * r.s - xy * r.yz,
            s * r.xy - yz * r.zx + zx * r.yz + xy * r.s);
    }
    rotor  operator /  (const rotor& r) const { return (*this) * r.inv(); }

    rotor  operator *  (const T& t) const { return rotor(s * t, yz * t, zx * t, xy * t); }
    rotor  operator /  (const T& t) const { return rotor(s / t, yz / t, zx / t, xy / t); }

    rotor unit() const { return *this / abs1(); }
    rotor inv()  const { return rev() / abs2(); }
    rotor rev() const { return rotor(s, -yz, -zx, -xy); }
    rotor conj() const { return rotor(s, -yz, -zx, -xy); }
    T     abs2() const { return s * s + xy * xy + yz * yz + zx * zx; }
    T     abs1() const { return sqrtT(abs2()); }
};

template<typename T>
rotor<T> operator * (const T& t, const rotor<T> r) { return r * t; }
template<typename T>
rotor<T> operator / (const T& t, const rotor<T> r) { return t * r.inv(); }

template<typename T, typename I>
struct castT<rotor<T>, I> {
    using rotor_t = ::altruct::math::cl3::rotor<T>;
    static rotor_t of(const I& s) {
        return rotor_t(castOf<T>(s));
    }
    static rotor_t of(const rotor_t& ref, const I& s) {
        return rotor_t(castOf(ref.s, s));
    }
};
template<typename T>
struct castT<rotor<T>, rotor<T>> : nopCastT<::altruct::math::cl3::rotor<T>> {};

template<typename T>
struct identityT<rotor<T>> {
    using rotor_t = ::altruct::math::cl3::rotor<T>;
    static rotor_t of(const rotor_t& r) {
        return rotor_t(identityOf(r.s), zeroOf(r.yz), zeroOf(r.zx), zeroOf(r.xy));
    }
};

template<typename T>
struct zeroT<rotor<T>> {
    using rotor_t = ::altruct::math::cl3::rotor<T>;
    static rotor_t of(const rotor_t& r) {
        return rotor_t(zeroOf(r.s), zeroOf(r.yz), zeroOf(r.zx), zeroOf(r.xy));
    }
};

template<typename T>
struct conjugateT<rotor<T>> {
    using rotor_t = ::altruct::math::cl3::rotor<T>;
    static rotor_t of(const rotor_t& r) {
        return r.conj();
    }
};

//------------------------------------------------------------------------------------------------/
// Odd part of G(3,0)
template<typename T>
class vector {
public:
    T x, y, z; // vectors
    T w;       // trivector (volume)

    vector() : x(0), y(0), z(0), w(0) {}
    vector(T x, T y, T z) : x(std::move(x)), y(std::move(y)), z(std::move(z)), w(zeroOf(z)) {}
    vector(T x, T y, T z, T w) : x(std::move(x)), y(std::move(y)), z(std::move(z)), w(std::move(w)) {}

    bool operator == (const vector& v) const { return (x == v.x) && (y == v.y) && (z == v.z) && (w == v.w); }
    bool operator != (const vector& v) const { return !(*this == v); }
    bool operator <  (const vector& v) const {
        if (x != v.x) return x < v.x;
        if (y != v.y) return y < v.y;
        if (z != v.z) return z < v.z;
        return w < v.w;
    }
    bool operator >  (const vector& v) const { return (v < *this); }
    bool operator <= (const vector& v) const { return !(v < *this); }
    bool operator >= (const vector& v) const { return !(*this < v); }

    vector& operator += (const vector& v) { x += v.x; y += v.y; z += v.z; w += v.w; return *this; }
    vector& operator -= (const vector& v) { x -= v.x; y -= v.y; z -= v.z; w -= v.w; return *this; }

    vector& operator *= (const rotor<T>& r) { return *this = *this * r; }
    vector& operator /= (const rotor<T>& r) { return *this = *this / r; }

    vector& operator *= (const T& t) { x *= t; y *= t; z *= t; w *= t; return *this; }
    vector& operator /= (const T& t) { x /= t; y /= t; z /= t; w /= t; return *this; }

    vector   operator +  (const vector& v) const { return vector(x + v.x, y + v.y, z + v.z, w + v.w); }
    vector   operator -  (const vector& v) const { return vector(x - v.x, y - v.y, z - v.z, w - v.w); }
    vector   operator -  ()                const { return vector(-x, -y, -z, -w); }

    vector   operator *  (const rotor<T>& r)  const {
        return vector(x * r.s - y * r.xy + z * r.zx - w * r.yz,
            x * r.xy + y * r.s - z * r.yz - w * r.zx,
            -x * r.zx + y * r.yz + z * r.s - w * r.xy,
            x * r.yz + y * r.zx + z * r.xy + w * r.s);
    }
    vector   operator /  (const rotor<T>& r) const { return *this * r.inv(); }

    rotor<T> operator *  (const vector& v)  const {
        return rotor<T>(x * v.x + y * v.y + z * v.z - w * v.w,
            x * v.w + y * v.z - z * v.y + w * v.x,
            -x * v.z + y * v.w + z * v.x + w * v.y,
            x * v.y - y * v.x + z * v.w + w * v.z);
    }
    rotor<T> operator /  (const vector& v) const { return *this * v.inv(); }

    vector   operator *  (const T& t) const { return vector(x * t, y * t, z * t, w * t); }
    vector   operator /  (const T& t) const { return vector(x / t, y / t, z / t, w / t); }

    vector unit() const { return *this / abs1(); }
    vector inv()  const { return rev() / abs2(); }
    vector rev() const { return vector(x, y, z, -w); }
    vector conj() const { return vector(-x, -y, -z, w); }
    T      abs2() const { return x * x + y * y + z * z + w * w; }
    T      abs1() const { return sqrtT(abs2()); }

    vector reflect(const vector<T>& v) const { return -v * (*this) * v.rev(); } // v must be normalized
    vector rotate(const rotor<T>& r)  const { return  r * (*this) * r.rev(); } // r must be normalized
};

template<typename T>
vector<T> operator * (const T& t, const vector<T> v) { return v * t; }
template<typename T>
vector<T> operator / (const T& t, const vector<T> v) { return t * v.inv(); }

template<typename T>
vector<T> operator * (const rotor<T> r, const vector<T>& v) {
    return vector<T>(r.s * v.x - r.yz * v.w - r.zx * v.z + r.xy * v.y,
        r.s * v.y + r.yz * v.z - r.zx * v.w - r.xy * v.x,
        r.s * v.z - r.yz * v.y + r.zx * v.x - r.xy * v.w,
        r.s * v.w + r.yz * v.x + r.zx * v.y + r.xy * v.z);
}
template<typename T>
vector<T> operator / (const rotor<T> r, const vector<T>& v) { return r * v.inv(); }

template<typename T>
vector<T> make_reflector(const vector<T>& over) { return over.unit(); }

template<typename T>
rotor<T> make_rotor(const vector<T>& from, const vector<T>& to) {
    T id = identityOf(from.w);
    return (to.unit() * from.unit() + id).unit();
}

template<typename T, typename I>
struct castT<vector<T>, I> {
    using vector_t = ::altruct::math::cl3::vector<T>;
    static vector_t of(const I& w) {
        T s = castOf<T>(w);
        return vector_t(zeroOf(s), zeroOf(s), zeroOf(s), s);
    }
    static vector_t of(const vector_t& ref, const I& w) {
        T s = castOf(ref.w, w);
        return vector_t(zeroOf(s), zeroOf(s), zeroOf(s), s);
    }
};
template<typename T>
struct castT<vector<T>, vector<T>> : nopCastT<::altruct::math::cl3::vector<T>> {};

template<typename T>
struct identityT<vector<T>> {
    using vector_t = ::altruct::math::cl3::vector<T>;
    static vector_t of(const vector_t& v) {
        // there is no identity because vector * vector = rotor != vector; return zero
        return zeroOf(v);
    }
};

template<typename T>
struct zeroT<vector<T>> {
    using vector_t = ::altruct::math::cl3::vector<T>;
    static vector_t of(const vector_t& v) {
        return vector_t(zeroOf(v.x), zeroOf(v.y), zeroOf(v.z), zeroOf(v.w));
    }
};

template<typename T>
struct conjugateT<vector<T>> {
    using vector_t = ::altruct::math::cl3::vector<T>;
    static vector_t of(const vector_t& v) {
        return v.conj();
    }
};

//------------------------------------------------------------------------------------------------/
// Algebra of G(3,0)
template<typename T>
class multivector {
public:
    rotor<T> r; // even part
    vector<T> v;  // odd part

    multivector() {}
    multivector(rotor<T> r) : r(std::move(r)) {}
    multivector(vector<T> v) : r(zeroOf(v.w)), v(std::move(v)) {}
    multivector(rotor<T> r, vector<T> v) : r(std::move(r)), v(std::move(v)) {}

    bool operator == (const multivector& m) const { return (r == m.r) && (v == m.v); }
    bool operator != (const multivector& m) const { return !(*this == m); }
    bool operator <  (const multivector& m) const { return (r != m.r) ? (r < m.r) : (v < m.v); }
    bool operator >  (const multivector& m) const { return (m < *this); }
    bool operator <= (const multivector& m) const { return !(m < *this); }
    bool operator >= (const multivector& m) const { return !(*this < m); }

    multivector& operator += (const rotor<T>& r2) { r += r2; return *this; }
    multivector& operator -= (const rotor<T>& r2) { r -= r2; return *this; }

    multivector& operator += (const vector<T>& v2) { v += v2; return *this; }
    multivector& operator -= (const vector<T>& v2) { v -= v2; return *this; }

    multivector& operator += (const multivector& m) { r += m.r; v += m.v; return *this; }
    multivector& operator -= (const multivector& m) { r -= m.r; v -= m.v; return *this; }

    multivector& operator *= (const rotor<T>& r2) { r *= r2; v *= r2; return *this; }
    multivector& operator /= (const rotor<T>& r2) { r /= r2; v /= r2; return *this; } // return *this *= r2.inv();

    multivector& operator *= (const vector<T>& v2) { return *this = *this * v2; }
    multivector& operator /= (const vector<T>& v2) { return *this = *this / v2; }

    multivector& operator *= (const multivector& m) { return *this = *this * m; }
    multivector& operator /= (const multivector& m) { return *this = *this / m; }

    multivector& operator *= (const T& t) { r *= t; v *= t; return *this; }
    multivector& operator /= (const T& t) { r /= t; v /= t; return *this; }

    multivector  operator +  (const rotor<T>& r2) const { return multivector(r + r2, v); }
    multivector  operator -  (const rotor<T>& r2) const { return multivector(r - r2, v); }

    multivector  operator +  (const vector<T>& v2) const { return multivector(r, v + v2); }
    multivector  operator -  (const vector<T>& v2) const { return multivector(r, v - v2); }

    multivector  operator +  (const multivector& m) const { return multivector(r + m.r, v + m.v); }
    multivector  operator -  (const multivector& m) const { return multivector(r - m.r, v - m.v); }

    multivector  operator -  ()                     const { return multivector(-r, -v); }

    multivector  operator *  (const rotor<T>& r2)  const { multivector res = *this; res *= r2; return res; }
    multivector  operator /  (const rotor<T>& r2)  const { multivector res = *this; res /= r2; return res; }

    multivector  operator *  (const vector<T>& v2)  const { return multivector(v * v2, r * v2); }
    multivector  operator /  (const vector<T>& v2)  const { return multivector(v / v2, r / v2); } // return *this * v2.inv();

    multivector  operator *  (const multivector& m)  const { return multivector(r * m.r + v * m.v, r * m.v + v * m.r); }
    multivector  operator /  (const multivector& m)  const { return *this * m.inv(); }

    multivector  operator *  (const T& t) const { return multivector(r * t, v * t); }
    multivector  operator /  (const T& t) const { return multivector(r / t, v / t); }

    multivector inv() const {
        if (r.s == 0 && r.yz == 0 && r.zx == 0 && r.xy == 0) return v.inv();
        // {(r - v r^-1 v)^-1, (v - r v^-1 r)^-1}
        auto t = -v * r.inv();
        auto ri = (r + t * v).inv();
        auto vi = ri * t;
        return multivector(ri, vi);
    }
    multivector rev() const {
        return multivector(r.rev(), v.rev());
    }
    
    multivector conj() const {
        return multivector(r.conj(), v.conj());
    }
};

template<typename T>
multivector<T> operator + (const rotor<T>& r, const vector<T> v) { return multivector<T>(r, v); }
template<typename T>
multivector<T> operator + (const vector<T>& v, const rotor<T> r) { return multivector<T>(r, v); }
template<typename T>
multivector<T> operator - (const rotor<T>& r, const vector<T> v) { return multivector<T>(r, -v); }
template<typename T>
multivector<T> operator - (const vector<T>& v, const rotor<T> r) { return multivector<T>(-r, v); }

template<typename T>
multivector<T> operator + (const rotor<T>& r, const multivector<T> m) { return multivector<T>(r + m.r, m.v); }
template<typename T>
multivector<T> operator + (const vector<T>& v, const multivector<T> m) { return multivector<T>(m.r, v + m.v); }
template<typename T>
multivector<T> operator - (const rotor<T>& r, const multivector<T> m) { return multivector<T>(r - m.r, -m.v); }
template<typename T>
multivector<T> operator - (const vector<T>& v, const multivector<T> m) { return multivector<T>(-m.r, v - m.v); }

template<typename T>
multivector<T> operator * (const T& t, const multivector<T> m) { return m * t; }
template<typename T>
multivector<T> operator / (const T& t, const multivector<T> m) { return t * m.inv(); }
template<typename T>
multivector<T> operator * (const rotor<T>& r, const multivector<T> m) { return multivector<T>(r * m.r, r * m.v); }
template<typename T>
multivector<T> operator / (const rotor<T>& r, const multivector<T> m) { return r * m.inv(); }
template<typename T>
multivector<T> operator * (const vector<T>& v, const multivector<T> m) { return multivector<T>(v * m.v, v * m.r); }
template<typename T>
multivector<T> operator / (const vector<T>& v, const multivector<T> m) { return v * m.inv(); }

template<typename T, typename I>
struct castT<multivector<T>, I> {
    using multivector_t = ::altruct::math::cl3::multivector<T>;
    static multivector_t of(const I& s) {
        return multivector_t(castOf<T>(s));
    }
    static multivector_t of(const multivector_t& ref, const I& s) {
        return multivector_t(castOf(ref.r, s));
    }
};
template<typename T>
struct castT<multivector<T>, multivector<T>> : nopCastT<::altruct::math::cl3::multivector<T>> {};

template<typename T>
struct identityT<multivector<T>> {
    using multivector_t = ::altruct::math::cl3::multivector<T>;
    static multivector_t of(const multivector_t& m) {
        return multivector_t(identityOf(m.r));
    }
};

template<typename T>
struct zeroT<multivector<T>> {
    using multivector_t = ::altruct::math::cl3::multivector<T>;
    static multivector_t of(const multivector_t& m) {
        return multivector_t(zeroOf(m.r));
    }
};

template<typename T>
struct conjugateT<multivector<T>> {
    using multivector_t = ::altruct::math::cl3::multivector<T>;
    static multivector_t of(const multivector_t& m) {
        return m.conj();
    }
};

} // cl3
} // math
} // altruct
