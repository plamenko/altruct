#pragma once

#include <type_traits>
#include <vector>

#include "algorithm/math/base.h"

namespace altruct {
namespace math {

/**
 * Based on the:
 *   https://en.wikipedia.org/wiki/Pythagorean_triple
 *   https://en.wikipedia.org/wiki/Eisenstein_triple
 *   https://en.wikipedia.org/wiki/Integer_triangle#Integer_triangles_with_a_120.C2.B0_angle
 */

/**
 * A triple.
 */
template<typename T>
struct triple {
    T a, b, c;
    bool operator < (const triple& rhs) const {
        if (a != rhs.a) return a < rhs.a;
        if (b != rhs.b) return b < rhs.b;
                        return c < rhs.c;
    }
    bool operator == (const triple& rhs) const {
        return a == rhs.a && b == rhs.b && c == rhs.c;
    }
};

/**
 * Generates all the Pythagorean triples up to the specified limit.
 * (Sides of a triangle where one of the angles is 90 degrees.)
 *
 * Pythagorean triple (a, b, c) satisfies:
 *   0 < a < b < c <= c_max
 *   a^2 + b^2 = c^2
 */
template<typename I, typename F>
void pythagorean_triples(I c_max, bool only_primitive, F visitor) {
    // a' = m^2 - n^2, b' = 2mn, c' = m^2 + n^2
    // a = k a', b = k b', c = k c'
    I m_max = sqrtT<I>(c_max - 1);
    for (I m = 1; m <= m_max; m++) {
        I m2 = sqT<I>(m);
        I n_min = (m % 2) + 1; // different parity
        I n_max = min(m - 1, sqrtT<I>(c_max - m2));
        for (I n = n_min; n <= n_max; n += 2) {
            if (gcd(m, n) != 1) continue;
            I n2 = sqT<I>(n), mn = m * n;
            I a = m2 - n2, b = mn * 2, c = m2 + n2;
            if (a > b) swap(a, b);
            I k_max = only_primitive ? 1 : c_max / c;
            I ka = 0, kb = 0, kc = 0;
            for (I k = 1; k <= k_max; k++) {
                ka += a, kb += b, kc += c;
                visitor(ka, kb, kc);
            }
        }
    }
}
template<typename I>
std::vector<triple<I>> pythagorean_triples(I c_max, bool only_primitive) {
    std::vector<triple<I>> vt;
    auto visitor = [&](I a, I b, I c){ vt.push_back({ a, b, c }); };
    pythagorean_triples(c_max, only_primitive, visitor);
    return vt;
}

/**
 * Generates all the Pythagorean triples with one leg fixed.
 *
 * @param f - factorization of the fixed leg
 */
template<typename I, typename P, typename F>
void pythagorean_triples_fixed_leg(I leg, const std::vector<std::pair<P, int>>& f, F visitor) {
    // a^2 = c^2 - b^2 = (c-b)(c+b) = d*e
    I leg2 = sqT<I>(leg);
    auto f2 = f;
    for (auto& t : f2) {
        t.second *= 2;
    }
    vector<I> vd; divisors(vd, f2);
    for (const auto& d : vd) {
        I e = leg2 / d;
        I c2 = e + d, b2 = e - d;
        if (b2 <= 0 || b2 % 2 != 0) continue;
        visitor(b2 / 2, c2 / 2);
    }
}
template<typename I, typename P>
std::vector<triple<I>> pythagorean_triples_fixed_leg(I leg, const std::vector<std::pair<P, int>>& f) {
    std::vector<triple<I>> vt;
    auto visitor = [&](I b, I c){ vt.push_back({ leg, b, c }); };
    pythagorean_triples_fixed_leg(leg, f, visitor);
    return vt;
}

/**
 * Generates all the Eisenstein 60 triples up to the specified limit.
 * (Sides of a triangle where one of the angles is 60 degrees.)
 *
 * Eisenstein triple (a, b, c) satisfies:
 *   0 < a < b < c <= c_max
 *   a^2 - ac + c^2 = b^2  // 60 degrees (oposite to the side b)
 */ 
template<typename I, typename F>
void eisenstein_triples60(I c_max, bool only_primitive, F visitor) {
    // a = (m^2 - n^2), b = (m^2 - mn + n^2), c = (2mn - n^2)
    // With coprime integers m, n with 0 < n < m.
    // The angle of 60 degrees is opposite to the side b.
    // If m+n==0 (mod 3), then gcd(a,b,c)=3, else gcd(a,b,c)=1.
    // Two different pairs (m, n) and (m, m-n) generate the same triple.
    // Unfortunately the two pairs can both be of gcd=3, so we can't simply
    // skip that case. Instead, duplicates can be simply avoided by n going
    // only till m/2. (The only solution for n==m/2 is (3,3,3) for m=2,n=1.)
    I m_max = min(sqrtT<I>(c_max * 4), (c_max * 3 + 1) / 2);
    for (I m = 1; m <= m_max; m++) {
        I m2 = sqT<I>(m);
        for (I n = 1; n <= m / 2; n++) {
            if (gcd(m, n) != 1) continue;
            I n2 = sqT<I>(n), mn = m * n;
            I a = m2 - n2, b = m2 - mn + n2, c = mn * 2 - n2;
            if (a > c) swap(a, c);
            if ((m + n) % 3 == 0) {
                a /= 3, b /= 3, c /= 3;
            }
            if (c > c_max) continue;
            I k_max = only_primitive ? 1 : c_max / c;
            I ka = 0, kb = 0, kc = 0;
            for (I k = 1; k <= k_max; k++) {
                ka += a, kb += b, kc += c;
                visitor(ka, kb, kc);
            }
        }
    }
}
template<typename I>
std::vector<triple<I>> eisenstein_triples60(I c_max, bool only_primitive) {
    std::vector<triple<I>> vt;
    auto visitor = [&](I a, I b, I c){ vt.push_back({ a, b, c }); };
    eisenstein_triples60(c_max, only_primitive, visitor);
    return vt;
}

/**
 * Generates all the Eisenstein 120 triples up to the specified limit.
 * (Sides of a triangle where one of the angles is 120 degrees.)
 *
 * Eisenstein triple (a, b, c) satisfies:
 *   0 < a < b < c <= c_max
 *   a^2 + ab + b^2 = c^2  // 120 degrees (oposite to the side c)
 */
template<typename I, typename F>
void eisenstein_triples120(I c_max, bool only_primitive, F visitor) {
    // a = (m^2 - n^2), b = (2mn + n^2), c = (m^2 + mn + n^2)
    // With coprime integers m, n with 0 < n < m.
    // The angle of 120 degrees is opposite to the side c.
    // If m-n==0 (mod 3), then gcd(a,b,c)=3, else gcd(a,b,c)=1.
    // However, unlike with the 60 degree case, here we can avoid duplicates
    // by simply skipping the gcd=3 case. This is because the biggest side c
    // can only be generated with a single pair (m, n). This means that each
    // primitive triple can be generated in precisely two ways: once directly
    // with gcd=1, and once indirectly with gcd=3 (which we can then skip).
    I m_max = sqrtT<I>(c_max);
    for (I m = 1; m <= m_max; m++) {
        I m2 = sqT<I>(m);
        for (I n = 1; n < m; n++) {
            if ((m - n) % 3 == 0) continue;
            if (gcd(m, n) != 1) continue;
            I n2 = sqT<I>(n), mn = m * n;
            I a = m2 - n2, b = mn * 2 + n2, c = m2 + mn + n2;
            if (a > b) swap(a, b);
            if (c > c_max) continue;
            I k_max = only_primitive ? 1 : c_max / c;
            I ka = 0, kb = 0, kc = 0;
            for (I k = 1; k <= k_max; k++) {
                ka += a, kb += b, kc += c;
                visitor(ka, kb, kc);
            }
        }
    }
}
template<typename I>
std::vector<triple<I>> eisenstein_triples120(I c_max, bool only_primitive) {
    std::vector<triple<I>> vt;
    auto visitor = [&](I a, I b, I c){ vt.push_back({ a, b, c }); };
    eisenstein_triples120(c_max, only_primitive, visitor);
    return vt;
}

} // math
} // altruct
