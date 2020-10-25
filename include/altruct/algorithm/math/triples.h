#pragma once

#include <type_traits>
#include <vector>

#include "altruct/algorithm/math/base.h"
#include "altruct/algorithm/math/primes.h"

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
 * Minimum element of the triple.
 */
template<typename I>
I min_element(const triple<I>& t) {
    return std::min(std::min(t.a, t.b), t.c);
}

/**
 * Maximum element of the triple.
 */
template<typename I>
I max_element(const triple<I>& t) {
    return std::max(std::max(t.a, t.b), t.c);
}

/**
 * Returns whether the triple is primitive.
 */
template<typename I>
bool is_primitive(const triple<I>& t) {
    return gcd(gcd(t.a, t.b), t.c) == 1;
}

/**
 * Returns the next Pythagorean triple given the transformation matrix and the current triple.
 */
template<typename I>
triple<I> next_pythagorean_triple(const I(&m)[3][3], const triple<I>& t) {
    return triple<I>{
        m[0][0] * t.a + m[0][1] * t.b + m[0][2] * t.c,
            m[1][0] * t.a + m[1][1] * t.b + m[1][2] * t.c,
            m[2][0] * t.a + m[2][1] * t.b + m[2][2] * t.c,
    };
}

/**
 * Returns whether the triple is a generator with invariant a^2 + b^2 + k = c^2.
 */
template<typename I>
bool is_pythagorean_triple_generator(const triple<I>& t) {
    static const I mi[3][3][3] = {
        { {-2, -1, 2}, {1, 2, -2}, {-2, -2, 3} },
        { {1, 2, -2}, {-2, -1, 2}, {-2, -2, 3} },
        { {2, 1, -2}, {1, 2, -2}, {-2, -2, 3} },
    };
    return is_primitive(t) &&
        min_element(next_pythagorean_triple(mi[0], t)) < 0 &&
        min_element(next_pythagorean_triple(mi[1], t)) < 0 &&
        min_element(next_pythagorean_triple(mi[2], t)) < 0;
}

/**
 * Generators of the Pythagorean triples with invariant a^2 + b^2 + k = c^2.
 *
 * @param k - some fixed non-negative k used in the invariant above.
 */
template<typename I>
std::vector<triple<I>> pythagorean_triples_generators(I k) {
    vector<triple<int>> vg;
    if (k == 0) {
        vg.push_back({ 3, 4, 5 });
    }
    else if (k == 1) {
        vg.push_back({ 2, 2, 3 });
    }
    else if (k >= 2) {
        int a_max = isqrt(k) * 2;
        for (int a = 0; a <= a_max; a++) {
            int b_max = k / max(a, 1) + a;
            for (int b = a; b <= b_max; b++) {
                int64_t c2 = isq(a) + isq(b) + k; // TODO: isq, isqrt
                int c = isqrt(c2);
                triple<int> t{ a, b, c };
                if (isq(c) == c2 && is_pythagorean_triple_generator(t)) vg.push_back(t);
            }
        }
    }
    return vg;
}

/**
 * Generates all the primitive Pythagorean triples with invariant a^2 + b^2 + k = c^2.
 *
 * @param c_max - generates triples up to c_max
 * @param generators - generators for the fixed k
 */
template<typename I>
std::vector<triple<I>> generate_pythagorean_triples(I n, const std::vector<triple<I>>& generators) {
    static const I m[3][3][3] = {
        { {-2, 1, 2}, {-1, 2, 2}, {-2, 2, 3} },
        { {1, -2, 2}, {2, -1, 2}, {2, -2, 3} },
        { {2, 1, 2}, {1, 2, 2}, {2, 2, 3} },
    };
    std::vector<triple<I>> vt = generators;
    for (size_t i = 0; i < vt.size(); i++) {
        for (int j = 0; j < 3; j++) {
            if (j == 1 && vt[i].a == vt[i].b) continue; // same triple as for j = 0
            if (j == 2 && vt[i].a == 0) continue; // same triple as for j = 0
            auto t = next_pythagorean_triple(m[j], vt[i]);
            if (max_element(t) <= n && is_primitive(t)) vt.push_back(t);
        }
    }
    return vt;
}

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
        I n_max = std::min(m - 1, sqrtT<I>(c_max - m2));
        for (I n = n_min; n <= n_max; n += 2) {
            if (gcd(m, n) != 1) continue;
            I n2 = sqT<I>(n), mn = m * n;
            I a = m2 - n2, b = mn * 2, c = m2 + n2;
            if (a > b) std::swap(a, b);
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
    std::vector<I> vd; divisors(vd, f2);
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
    I m_max = std::min(sqrtT<I>(c_max * 4), (c_max * 3 + 1) / 2);
    for (I m = 1; m <= m_max; m++) {
        I m2 = sqT<I>(m);
        for (I n = 1; n <= m / 2; n++) {
            if (gcd(m, n) != 1) continue;
            I n2 = sqT<I>(n), mn = m * n;
            I a = m2 - n2, b = m2 - mn + n2, c = mn * 2 - n2;
            if (a > c) std::swap(a, c);
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
            if (a > b) std::swap(a, b);
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
