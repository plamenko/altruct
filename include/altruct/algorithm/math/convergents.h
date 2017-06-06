#pragma once

#include "altruct/algorithm/math/base.h"

#include <vector>
#include <utility>
#include <algorithm>

namespace altruct {
namespace math {

/**
 * Convergent of `sqrt(S)` with denominator bigger than `Q`.
 *
 * Complexity: O(log Q)
 */
template<typename I>
std::pair<I, I> sqrt_convergent(I S, I Q) {
    I a0 = sqrtT(S);
    if (a0 * a0 == S) return{ a0, 1 };
    I po = 1, pn = 0;
    I qo = 0, qn = 1;
    I m = 0, d = 1, a = a0, t;
    for (;;) {
        t = pn + a * po, pn = po, po = t;
        t = qn + a * qo, qn = qo, qo = t;
        if (qo > Q) break;
        I m_n = d * a - m;
        I d_n = (S - m_n * m_n) / d;
        I a_n = (a0 + m_n) / d_n;
        m = m_n, d = d_n, a = a_n;
    }
    return{ po, qo }; // qn <= Q < qo
}

/**
 * Continued fraction from a rational number.
 *
 * I.e. `p/q = [a0, a1, ..., an]`
 *
 * Complexity: O(log q)
 */
template<typename I>
std::vector<I> continued_fraction(I p, I q) {
    std::vector<I> va;
    while (q != 0) {
        va.push_back(p / q);
        I r = p % q; p = q; q = r;
    }
    return va;
}

/**
 * Convergents of a continued fraction.
 *
 * @param semi_convergents: number of semi-convergents per convergent.
 *   Pass infinity to get all of them (i.e. all the best approximations),
 *   or zero (default) to get none (returns only the proper convergents).
 */
template<typename I>
std::vector<std::pair<I, I>> convergents(const std::vector<I> &va, I semi_convergents = I(0)) {
    std::vector<std::pair<I, I>>vpq;
    I p2 = 0, p1 = 1, p = 0;
    I q2 = 1, q1 = 0, q = 1;
    for (const auto& a_i : va) {
        // Technically, if a_i is even, a_i / 2 is admissible only if the
        // corresponding semiconvergent is better than the previous convergent.
        // We do not check for it here however.
        I a0 = a_i - std::min<I>(semi_convergents, a_i / 2);
        for (I a = a0; a <= a_i; a++) {
            p = a * p1 + p2;
            q = a * q1 + q2;
            vpq.push_back({ p, q });
        }
        p2 = p1, p1 = p;
        q2 = q1, q1 = q;
    }
    return vpq;
}

/**
 * Finds the lattice point `{x, y}` closest to the line `A x + B y + C`,
 * where `x` is within the interval `[x_min, x_max]`.
 * I.e. `x` within `[x_min, x_max]` for `{x, y}` that minimizes `abs(A x + B y + C)`.
 * For any given `x`, `y = round((A x + C) / -B)`, so this essentially solves:
 * `x` within `[x_min, x_max]` that minimizes `abs(A x + C + B round((A x + C) / -B))`.
 */
template<typename I>
I line_closest_lattice_point(I a, I b, I c, I x_min, I x_max) {
    if (x_min >= x_max) return x_min;
    if (a == 0) return x_min;
    if (b == 0) return boundT<I>(div_round<I>(c, -a), x_min, x_max);
    if (a < 0) a = -a, c = -c;
    if (b < 0) b = -b;
    a %= b;
    if (a == 0) return x_min;
    auto dist = [&](I x, I y) { return absT<I>(a * x + b * y + c); };
    // Reciprocally, for any fixed `round((A x + C) / -B) = y`,
    // the best `x` is `x = round((B y + C) / -A)`. This gives
    // `B y + C + A round((B y + C) / -A)`, which is the same
    // problem but with reduced limits. In particular:
    // {A, B, C} -> {B, A%B, C}, like in Euclidean algorithm,
    // and size of range for `y` is (A%B / B) of that of `x`.
    I y_min = div_ceil<I>(a * (x_max * 2 + 1) + c * 2, -b * 2);
    I y_max = div_floor<I>(a * (x_min * 2 - 1) + c * 2, -b * 2);
    I y = line_closest_lattice_point<I>(b, -a, c, y_min, y_max);
    I x = boundT<I>(div_round<I>(b * y + c, -a), x_min, x_max);
    if (x != x_min && dist(x_min, y_min) < dist(x, y)) x = x_min, y = y_min;
    if (x != x_max && dist(x_max, y_max) < dist(x, y)) x = x_max, y = y_max;
    return x;
}

/**
 * Finds `x` within `[x_min, x_max]` that minimizes `A x + B floor((C x + D) / E)`.
 * Similar to `line_closest_lattice_point` but a slighlty different condition.
 * When plotted points look like a ladder hence the function name.
 */
template<typename I>
I minimize_floor_ladder(I a, I b, I c, I d, I e, I x_min, I x_max) {
    if (x_min >= x_max) return x_min;
    if (e < 0) c = -c, d = -d, e = -e;
    if (c < 0) b = -b, c = -c, d = e - d - 1;
    a += b * (c / e), c %= e;
    I x_ext = (a < 0) ? x_max : x_min;
    if (b == 0 || c == 0) return x_ext;
    I y_min = div_floor<I>(c * x_min + d, e);
    I y_max = div_floor<I>(c * x_max + d, e);
    if (y_min == y_max) return x_ext;
    I f = (a < 0) ? (e - d - 1) : (c - d - 1);
    if (a < 0) {
        I y = minimize_floor_ladder<I>(b, a, e, f, c, y_min, y_max - 1);
        I x = div_floor<I>(e * y + f, c);
        return ((a * x_max + b * y_max) < (a * x + b * y)) ? x_max : x;
    } else {
        I y = minimize_floor_ladder<I>(b, a, e, f, c, y_min + 1, y_max);
        I x = div_floor<I>(e * y + f, c);
        return ((a * x_min + b * y_min) < (a * x + b * y)) ? x_min : x;
    }
}

} // math
} // altruct
