#pragma once

#include "altruct/algorithm/math/base.h"
#include "altruct/structure/math/fraction.h"
#include "altruct/structure/math/modulo.h"

#include <unordered_map>
#include <vector>
#include <algorithm>

namespace altruct {
namespace math {

/**
 * The next element in the Farey sequence of order n.
 *
 * if f_prev is the left neighbour or -inf, f_next is the right neighbour of f
 * if f_prev is the right neighbour or +inf, f_next is the left neighbour of f
 */
template<typename I>
fraction<I> farey_neighbour(const I& n, const fraction<I>& f_prev, const fraction<I>& f) {
    I p = f_prev.p, q = f_prev.q;
    if (f_prev.q == 0) {
        gcd_ex(f.q, f.p, &p, &q);
        (f_prev.p < 0) ? p = -p : q = -q;
    }
    I k = (n + q) / f.q;
    return fraction<I>(k * f.p - p, k * f.q - q);
}

/**
 * Decimal expansion of a fraction p/q in base b.
 * Requirs 0 < p/q < 1 and b >= 2.
 */
template<typename I>
std::pair<std::vector<int>, size_t> repeating_decimal(int b, I p, I q) {
    std::vector<int> digits;
    size_t pos = 0;
    std::unordered_map<I, size_t> positions;
    while (!positions.count(p)) {
        positions[p] = pos;
        I bp = b * p;
        int d = castOf<int, I>(bp / q);
        digits.push_back(d);
        p = bp - d * q;
        if (p == 0) {
            return { digits, 0 };
        }
        pos += 1;
    }
    return { digits, pos - positions[p] };
}

/**
 * Gives n-th digit after decimal point in base b of a fraction p/q.
 * Requires 0 < p/q < 1.
 */
template<typename I>
int rational_digit(I n, int b, I p, I q, I* r = nullptr) {
    using modx = moduloX<I>;
    p = (powT(modx(b, q), n) * modx(p, q)).v;
    if (r) *r = p;
    return castOf<int, I>((p * b) / q);
}

/**
 * Gives digits n to n+len after decimal point in base b of a fraction p/q.
 * Requires 0 < p/q < 1.
 */
template<typename I>
std::vector<int> rational_digits(I n, int len, int b, I p, I q) {
    std::vector<int> digits;
    if (len-- > 0) {
        digits.push_back(rational_digit<I>(n, b, p, q, &p));
    }
    while (len-- > 0) {
        digits.push_back(rational_digit<I>(1, b, p, q, &p));
    }
    return digits;
}

} // math
} // altruct
