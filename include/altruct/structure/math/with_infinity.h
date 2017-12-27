#pragma once

#include "altruct/algorithm/math/base.h"

namespace altruct {
namespace math {

/**
 * Extends the underlying structure with the point at infinity.
 *
 * The following holds:
 *        +---+---+---+        +---+---+---+        +---+---+---+        +---+---+---+
 *      + | z | 0 |inf|      - | z | 0 |inf|      * | z | 0 |inf|      / | z | 0 |inf|
 *    +---+---+---+---+    +---+---+---+---+    +---+---+---+---+    +---+---+---+---+
 *    | z | z | z |inf|    | z | z | z |inf|    | z | z | 0 |inf|    | z | z |inf| 0 |
 *    +---+---+---+---+    +---+---+---+---+    +---+---+---+---+    +---+---+---+---+
 *    | 0 | z | 0 |inf|    | 0 | z | 0 |inf|    | 0 | 0 | 0 | ? |    | 0 | 0 | ? | 0 |
 *    +---+---+---+---+    +---+---+---+---+    +---+---+---+---+    +---+---+---+---+
 *    |inf|inf|inf|inf|    |inf|inf|inf| ? |    |inf|inf| ? |inf|    |inf|inf|inf| ? |
 *    +---+---+---+---+    +---+---+---+---+    +---+---+---+---+    +---+---+---+---+
 *
 *  Where:
 *    z   - a number different than 0 and inf
 *    0   - zero
 *    inf - infinity
 *    ?   - undefined
 *
 * @param T - the underlying type
 */
template<typename T>
class with_infinity {
public:
    T v; int is_inf;

    // construct from int, but only if T is not integral to avoid constructor clashing
    template <typename I = T, typename = std::enable_if_t<!std::is_integral<I>::value>>
    with_infinity(int v) : v(v), is_inf(0) {}
    with_infinity(const T& v = 0, int is_inf = 0) : v(v), is_inf(is_inf) {}
    with_infinity(const with_infinity& rhs) : v(rhs.v), is_inf(rhs.is_inf) {}

    bool operator == (const with_infinity &rhs) const { return (v == rhs.v) && (is_inf == rhs.is_inf); }
    bool operator != (const with_infinity &rhs) const { return !(*this == rhs); }
    bool operator <  (const with_infinity &rhs) const { return (is_inf != rhs.is_inf) ? (is_inf < rhs.is_inf) : (v < rhs.v); }
    bool operator >  (const with_infinity& rhs) const { return (rhs < *this); }
    bool operator <= (const with_infinity& rhs) const { return !(rhs < *this); }
    bool operator >= (const with_infinity& rhs) const { return !(*this < rhs); }

    with_infinity  operator +  (const with_infinity &rhs) const { with_infinity t(*this); t += rhs; return t; }
    with_infinity  operator -  (const with_infinity &rhs) const { with_infinity t(*this); t -= rhs; return t; }
    with_infinity  operator -  ()                         const { with_infinity t(-v);              return t; }
    with_infinity  operator *  (const with_infinity &rhs) const { with_infinity t(*this); t *= rhs; return t; }
    with_infinity  operator /  (const with_infinity &rhs) const { with_infinity t(*this); t /= rhs; return t; }
    with_infinity  operator %  (const with_infinity &rhs) const { with_infinity t(*this); t %= rhs; return t; }

    with_infinity& operator += (const with_infinity &rhs) { if (rhs.is_inf) is_inf = 1; else v += rhs.v; return *this; }
    with_infinity& operator -= (const with_infinity &rhs) { if (rhs.is_inf) is_inf = 1; else v -= rhs.v; return *this; }
    with_infinity& operator *= (const with_infinity &rhs) { if (rhs.is_inf) is_inf = 1; else v *= rhs.v; return *this; }
    with_infinity& operator /= (const with_infinity &rhs) { return *this *= rhs.inverse(); }
    with_infinity& operator %= (const with_infinity &rhs) { v %= rhs.v; return *this; }

    with_infinity inverse() const {
        if (is_inf) return with_infinity(zeroOf(v));
        if (v == zeroOf(v)) return with_infinity(v, 1);
        return with_infinity(identityOf(v) / v);
    }
};

template<typename T, typename I>
struct castT<with_infinity<T>, I> {
    typedef with_infinity<T> winf;
    static winf of(const I& x) {
        return winf(castOf<T>(x));
    }
    static winf of(const winf& ref, const I& x) {
        return winf(castOf(ref.v, x));
    }
};
template<typename T>
struct castT<with_infinity<T>, with_infinity<T>> : nopCastT<with_infinity<T>>{};

template<typename T>
struct identityT<with_infinity<T>> {
    typedef with_infinity<T> winf;
    static winf of(const winf& x) {
        return winf(identityOf(x.v));
    }
};

template<typename T>
struct zeroT<with_infinity<T>> {
    typedef with_infinity<T> winf;
    static winf of(const winf& x) {
        return winf(zeroOf(x.v));
    }
};

template<typename T>
struct infinityT<with_infinity<T>> {
    typedef with_infinity<T> winf;
    static bool is(const winf& x) { return x.is_inf != 0; }
};

template<typename T>
struct conjugateT<with_infinity<T>> {
    typedef with_infinity<T> winf;
    static winf of(const winf& x) {
        return winf(conjugateT<T>::of(x.v), x.is_inf);
    }
};

} // math
} // altruct
