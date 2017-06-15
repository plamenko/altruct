#pragma once

#include "base.h"

#include <set>
#include <vector>
#include <algorithm>
#include <unordered_map>

namespace altruct {
namespace math {

/**
 * Finds representation of prime `p` as a sum of two squares `a^2 + b^2 == p`.
 *
 * Complexity: O(sqrt P)
 */
template<typename I>
std::pair<I, I> squares_r_prime(I p, I hint_a = 0, I hint_b = 0) {
    static std::unordered_map<I, std::pair<I, I>> mss;
    if (sqT(hint_a) + sqT(hint_b) == p) return mss[p] = { hint_a, hint_b };
    if (mss.count(p)) return mss[p];
    for (I a = 1; sqT(a) * 2 <= p; a++) {
        I b = sqrtT(p - sqT(a));
        if (sqT(a) + sqT(b) == p) {
            return mss[p] = { a, b };
        }
    }
    return{ 0, 0 };
}

/**
 * Precomputes `squares_r_prime` up to `n`.
 *
 * Complexity: O(n)
 */
template<typename I>
void squares_r_prime_init(I n) {
    for (I a = 1; sqT(a) * 2 <= n; a++) {
        for (I b = a; sqT(a) + sqT(b) <= n; b++) {
            squares_r_prime(sqT(a) + sqT(b), a, b);
        }
    }
}

/**
 * Representations of `n` as a sum of two squares `a^2 + b^2 == n`.
 *
 * @param vf - prime factorization of n
 * @param unique_only - if true, the sign and order won't be taken into account
 *                      i.e. (1, 2) is considered the same as (-2, -1)
 */
template<typename P, typename I = P>
std::vector<std::pair<I, I>> squares_r_list(const std::vector<std::pair<P, int>> &vf, bool unique_only) {
    std::vector<std::pair<I, I>> v{ { 0, 1 } };
    I q = 1;
    for (const auto& f : vf) {
        I p = f.first; int e = f.second;
        if (p % 4 == 1) {
            auto cd = squares_r_prime(p);
            auto c = cd.first, d = cd.second;
            while (e-- > 0) {
                std::set<std::pair<I, I>> s;
                for (auto& t : v) {
                    I a = t.first, b = t.second;
                    I e1 = absT(a * c - b * d), f1 = (a * d + b * c);
                    if (e1 > f1) std::swap(e1, f1); s.insert({ e1, f1 });
                    I e2 = absT(a * d - b * c), f2 = (a * c + b * d);
                    if (e2 > f2) std::swap(e2, f2); s.insert({ e2, f2 });
                }
                v.assign(s.begin(), s.end());
            }
        } else if (p % 4 == 3) {
            if (e % 2 == 1) {
                return{};
            }
            q *= powT(p, e / 2);
        } else if (p == 2) {
            if (e % 2 == 1) {
                for (auto& t : v) {
                    t = { t.second - t.first, t.second + t.first };
                }
            }
            q *= powT(p, e / 2);
        }
    }
    for (auto& t : v) {
        t.first *= q, t.second *= q;
    }
    if (!unique_only) {
        for (int i = (int)v.size() - 1; i >= 0; i--) {
            I a = v[i].first, b = v[i].second;
            if (a != b) v.push_back({ b, a });
        }
        for (int i = (int)v.size() - 1; i >= 0; i--) {
            I a = v[i].first, b = v[i].second;
            if (a != 0) v.push_back({ -a, b });
        }
        for (int i = (int)v.size() - 1; i >= 0; i--) {
            I a = v[i].first, b = v[i].second;
            if (b != 0) v.push_back({ a, -b });
        }
    }
    std::sort(v.begin(), v.end());
    return v;
}

/**
 * Calculates the number of representations of n as a sum of 2 squares from a factorization of `n`.
 *
 * E.g. some unique representations:
 *   4 = 0^2 + 2^2
 *   5 = 1^2 + 2^2
 *   8 = 2^2 + 2^2
 *   25 = 3^2 + 4^2 = 0^2 + 5^2
 *
 * @param unique_only - if true, the sign and order won't be taken into account
 *                      i.e. (1, 2) is considered the same as (-2, -1)
 */
template<typename P, typename I = P>
I squares_r(const std::vector<std::pair<P, int>> &vf, bool unique_only) {
    I B = 1; int s = 1, q = 1;
    for (const auto& f : vf) {
        if (f.first % 4 == 1) {
            B *= f.second + 1;
        } else if (f.first % 4 == 3) {
            if (f.second % 2 == 1) B = 0;
        } else if (f.first == 2) {
            if (f.second % 2 == 1) s = -1;
        }
        if (f.second % 2 == 1) q = 0;
    }
    if (!unique_only) return B * 4;
    if (B % 2 == 1) B -= s;
    return B / 2 + q;
}

} // math
} // altruct
