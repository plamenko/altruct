#pragma once

#include "altruct/algorithm/math/base.h"

namespace altruct {
namespace math {

/**
 * Moebius Transformation {a, b, c, d, s}
 *
 * For s == +1:
 *          a z + b
 *   f(z) = -------
 *          c z + d
 *
 * For s == -1:
 *           a z* + b
 *   f'(z) = -------
 *           c z* + d
 *   where z* is conjugate of z, and
 *   f' is anti-moebius transformation (not to be confused with inverse)
 *
 * @param C - the underlying type (e.g. with_infinity<complex<double>> or with_infinity<complex<mod>>)
 */
template<typename C>
class moebius_tr {
public:
    C a, b, c, d; int s;

    // construct from int, but only if C is not integral to avoid constructor clashing
    template <typename I = C, typename = std::enable_if_t<!std::is_integral<I>::value>>
    moebius_tr(int a) : a(a), b(0), c(0), d(1), s(+1) {}
    moebius_tr(const C& a = 0) : a(a), b(zeroOf(a)), c(zeroOf(a)), d(identityOf(a)), s(+1) {}
    moebius_tr(const C& a, const C& b, const C& c, const C& d, int s = +1) : a(a), b(b), c(c), d(d), s(s) {}
    moebius_tr(const moebius_tr& rhs) : a(rhs.a), b(rhs.b), c(rhs.c), d(rhs.d), s(rhs.s) {}

    bool operator == (const moebius_tr& rhs) const { return (a == rhs.a && b == rhs.b && c == rhs.c && d == rhs.d && s == rhs.s); }
    bool operator != (const moebius_tr& rhs) const { return !(*this == rhs); }
    bool operator <  (const moebius_tr& rhs) const {
        if (a != rhs.a) return a < rhs.a;
        if (b != rhs.b) return b < rhs.b;
        if (c != rhs.c) return c < rhs.c;
        if (d != rhs.d) return d < rhs.d;
                        return s < rhs.s;
    }
    bool operator >  (const moebius_tr& rhs) const { return (rhs < *this); }
    bool operator <= (const moebius_tr& rhs) const { return !(rhs < *this); }
    bool operator >= (const moebius_tr& rhs) const { return !(*this < rhs); }

    moebius_tr  operator *  (const moebius_tr& rhs) const { moebius_tr t(*this); t *= rhs; return t; }
    moebius_tr  operator /  (const moebius_tr& rhs) const { moebius_tr t(*this); t /= rhs; return t; }  
    moebius_tr& operator *= (const moebius_tr& rhs) {
        return *this = moebius_tr(
            a * adapt(rhs.a) + b * adapt(rhs.c),
            a * adapt(rhs.b) + b * adapt(rhs.d),
            c * adapt(rhs.a) + d * adapt(rhs.c),
            c * adapt(rhs.b) + d * adapt(rhs.d),
            s * rhs.s);
    }
    moebius_tr& operator /= (const moebius_tr& rhs) { return *this *= rhs.inverse(); }

    moebius_tr& normalize() {
        if (a != zeroOf(a)) {
            auto ai = invert(a);
            a = identityOf(a), b *= ai, c *= ai, d *= ai;
        } else if (b != zeroOf(b)) {
            auto bi = invert(b);
            b = identityOf(b), c *= bi, d *= bi;
        } else {
            c = zeroOf(c), d = identityOf(d);
        }
        return *this;
    }

    moebius_tr inverse() const {
        C w = invert(a * d - b * c);
        return moebius_tr(adapt(d) * w, -adapt(b) * w, -adapt(c) * w, adapt(a) * w, s);
    }

    C operator()(C z) const {
        if (infinityT<C>::is(z)) return a / c;
        return (a * adapt(z) + b) / (c * adapt(z) + d);
    }

    static C invert(C z) { return identityOf(z) / z; }
    static C conj(C z) { return conjugateT<C>::of(z); }
    C adapt(C z) const { return (s < 0) ? conj(z) : z; }

    static moebius_tr scaling(C z) { return moebius_tr(z, zeroOf(z), zeroOf(z), identityOf(z), +1); }
    static moebius_tr translation(C z) { return moebius_tr(identityOf(z), z, zeroOf(z), identityOf(z), +1); }
    static moebius_tr rotation(C z) { return moebius_tr(z, zeroOf(z), zeroOf(z), identityOf(z), +1); }
    static moebius_tr flip_x(C id = 1) { return moebius_tr(-identityOf(id), zeroOf(id), zeroOf(id), identityOf(id), -1); }
    static moebius_tr flip_y(C id = 1) { return moebius_tr(identityOf(id), zeroOf(id), zeroOf(id), identityOf(id), -1); }
    static moebius_tr inversion(C id = 1) { return moebius_tr(zeroOf(id), identityOf(id), identityOf(id), zeroOf(id), -1); }
};

template<typename C, typename I>
struct castT<moebius_tr<C>, I> {
    typedef moebius_tr<C> moeb_tr;
    static moeb_tr of(const I& a) {
        return moeb_tr(castOf<C>(a));
    }
    static moeb_tr of(const moeb_tr& ref, const I& a) {
        return moeb_tr(castOf(ref.a, a));
    }
};
template<typename C>
struct castT<moebius_tr<C>, moebius_tr<C>> : nopCastT<moebius_tr<C>>{};

template<typename C>
struct identityT<moebius_tr<C>> {
    typedef moebius_tr<C> moeb_tr;
    static moeb_tr of(const moeb_tr& t) {
        return moeb_tr(identityOf(t.a));
    }
};

template<typename C>
struct zeroT<moebius_tr<C>> {
    typedef moebius_tr<C> moeb_tr;
    static moeb_tr of(const moeb_tr& t) {
        return moeb_tr(zeroOf(t.a));
    }
};

} // math
} // altruct
