#pragma once

#include "altruct/algorithm/math/base.h"

namespace altruct {
namespace math {

template <typename T>
class fraction {
public:
    T p, q;
    fraction(const T &p = 0) : p(p), q(identityT<T>::of(p)) { }
    fraction(const T &p, const T &q) : p(p), q(q) { reduce(); }

    void reduce() {
        T e0 = zeroT<T>::of(p), e1 = identityT<T>::of(p);
        T g = absT(gcd(p, q));
        if (g != e1) p /= g, q /= g;
        if (q < e0) p = -p, q = -q;
    }

    // WARNING: comparison may overflow!
    T cmp(const fraction& f) const { return p * f.q - f.p * q; }

    bool operator == (const fraction &f) const { return (p == f.p && q == f.q); }
    bool operator != (const fraction &f) const { return (p != f.p || q != f.q); }
    bool operator <  (const fraction &f) const { return cmp(f) <  0; }
    bool operator <= (const fraction &f) const { return cmp(f) <= 0; }
    bool operator >  (const fraction &f) const { return cmp(f) >  0; }
    bool operator >= (const fraction &f) const { return cmp(f) >= 0; }

    fraction& operator += (const fraction &f) { p = p * f.q + f.p * q; q *= f.q; reduce(); return *this; }
    fraction& operator -= (const fraction &f) { p = p * f.q - f.p * q; q *= f.q; reduce(); return *this; }
    fraction& operator *= (const fraction &f) { p = p * f.p;           q *= f.q; reduce(); return *this; }
    fraction& operator /= (const fraction &f) { T f_p = f.p; p *= f.q; q *= f_p; reduce(); return *this; }
    fraction& operator %= (const fraction &f) { p = zeroT<T>::of(p), q = identityT<T>::of(q); return *this; }

    fraction  operator +  (const fraction &f) const { return fraction(p * f.q + f.p * q, q * f.q); }
    fraction  operator -  (const fraction &f) const { return fraction(p * f.q - f.p * q, q * f.q); }
    fraction  operator -  (                 ) const { return fraction(-p, q); }
    fraction  operator *  (const fraction &f) const { return fraction(p * f.p,           q * f.q); }
    fraction  operator /  (const fraction &f) const { return fraction(p * f.q,           q * f.p); }
    fraction  operator %  (const fraction &f) const { return fraction(zeroT<T>::of(p)); }

    fraction reciprocal() const { return fraction(q, p); }
};

template<typename T>
struct identityT<fraction<T>> {
    static fraction<T> of(const fraction<T>& f) {
        return fraction<T>(identityT<T>::of(f.p));
    }
};

template<typename T>
struct zeroT<fraction<T>> {
    static fraction<T> of(const fraction<T>& f) {
        return fraction<T>(zeroT<T>::of(f.p));
    }
};

} // math
} // altruct
