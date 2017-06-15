#pragma once

#include "altruct/algorithm/math/base.h"

#include <vector>
#include <algorithm>

namespace altruct {
namespace math {

/**
 * Permutation represented in its cycle notation.
 *
 * Having an associative product, a neutral element, and inverses for all its elements,
 * makes the set of all permutations of S into a group, called the symmetric group of S.
 *
 * Note: elements are represented by numbers in the `[0, n)` range.
 *
 * @param I - integral type used for indices and values
 */
template<typename I>
class permutation {
public:
    typedef std::vector<I> line_t;
    typedef std::vector<I> cycle_t;
    typedef std::vector<cycle_t> cycles_t;
    typedef std::pair<I, I> transposition_t;
    typedef std::vector<transposition_t> transpositions_t;

    cycles_t cycles;
    I n;

    // identity
    permutation(I n = 0) : n(n) {}
    // from cycles notation
    permutation(const cycles_t& cycles, I n) : cycles(cycles), n(n) {}
    // from line notation
    permutation(const line_t& line) : cycles(line_to_cycles(line)), n((I)line.size()) {}
    // from transpositions
    permutation(const transpositions_t& transpositions, I n) : permutation(transpositions_to_line(transpositions, n)) {}

    // to cycles notation
    cycles_t to_cycles() const { return cycles; }
    // to cycles including 1-cycles
    cycles_t to_all_cycles() const { return all_cycles(cycles, n); }
    // to line notation
    line_t to_line() const { return cycles_to_line(cycles, n); }
    // to transpositions notation
    transpositions_t to_transpositions() const { return cycles_to_transpositions(cycles); }

    // comparison operators; note: size `n` of the list matters for the comparison
    bool operator == (const permutation& p2) const { return n == p2.n && cycles == p2.cycles; }
    bool operator != (const permutation& p2) const { return !(*this == p2); }
    bool operator <  (const permutation& p2) const { return (n != p2.n) ? (n < p2.n) : (cycles < p2.cycles); }
    bool operator >  (const permutation& p2) const { return (p2 < *this); }
    bool operator <= (const permutation& p2) const { return !(p2 < *this); }
    bool operator >= (const permutation& p2) const { return !(*this < p2); }

    // product is equal to function composition: (p1 * p2)(x) == p1(p2(x))
    permutation operator * (const permutation& p2) const { auto line = p2.to_line(); return permutation(apply_to(line)); }
    permutation& operator *= (const permutation& p2) { return *this = *this * p2; }

    // product by inverse
    permutation operator / (const permutation& p2) const { return *this * p2.inv(); }
    permutation& operator /= (const permutation& p2) { return *this = *this / p2; }

    // applies this permutation onto a line
    line_t& apply_to(line_t& line) const { return apply_cycles_to_line(line, cycles, n); }

    // this permutation applied `T` times
    permutation pow(int64_t T) const {
        if (T < 0) return inv().pow(-T);
        cycles_t vc;
        std::vector<bool> u(n);
        for (const auto& cycle : cycles) {
            I L = (I)cycle.size();
            for (I i = 0; i < L; i++) {
                cycle_t c;
                for (I j = i; !u[cycle[j]]; j = I((j + T) % L)) {
                    c.push_back(cycle[j]);
                    u[cycle[j]] = true;
                }
                if (c.size() > 1) vc.emplace_back(c);
            }
        }
        return permutation(vc, n);
    }

    // inverse
    permutation inv() const {
        cycles_t vc = cycles;
        for (auto& c : vc) {
            std::reverse(c.begin() + 1, c.end());
        }
        return permutation(vc, n);
    }

    // `T-th` root of this permutation, with specified parity
    // @param parity - {-1: optimal, 0: force_even, 1: force_odd}
    permutation root(int64_t T, int parity = -1) const {
        auto vc = to_all_cycles();
        std::vector<cycles_t> len_to_cycles(n + 1);
        for (const auto& c : vc) len_to_cycles[c.size()].push_back(c);
        // calc min number of transpositions, and parity changing length
        I parity_len = -1;
        I transpositions = 0;
        for (I l = 0; l <= n; l++) {
            if (len_to_cycles[l].empty()) continue;
            I m = (I)len_to_cycles[l].size();
            int64_t L = l * gcd_max<int64_t>(l, T);
            int64_t g = gcd(L, T), G = gcd(L * 2, T);
            if (m % g != 0) return permutation(); // root doesn't exist
            if (G == g * 2 && m >= G) parity_len = l;
            transpositions += I((L - 1) * (m / g));
        }
        int wrong_parity = (parity >= 0) && (transpositions % 2 != parity);
        if (wrong_parity && parity_len == -1) return permutation(); // root doesn't exist
        // build result
        permutation result(n);
        for (I l = 0; l <= n; l++) {
            if (len_to_cycles[l].empty()) continue;
            I m = (I)len_to_cycles[l].size();
            int64_t L = l * gcd_max<int64_t>(l, T);
            int64_t g = gcd(L, T), G;
            for (I i = 0; i < m; i += I(G)) {
                G = (wrong_parity && l == parity_len && i == 0) ? g * 2 : g;
                int64_t o = (T / G) % l;
                // combine G l-cycles into one G*l-cycle
                cycle_t c(G * l);
                for (I h = 0; h < l; h++) {
                    I k = I((h * o) % l * G);
                    for (I j = 0; j < G; j++) {
                        c[k + j] = len_to_cycles[l][i + j][h];
                    }
                }
                result.cycles.push_back(c);
            }
        }
        return result;
    }

    // static methods

    static cycles_t all_cycles(const cycles_t& cycles, I n) {
        cycles_t vc = cycles;
        std::vector<bool> u(n);
        for (const auto& cycle : cycles) {
            for (auto e : cycle) {
                u[e] = true;
            }
        }
        for (I i = 0; i < n; i++) {
            if (!u[i]) vc.push_back({ i });
        }
        return vc;
    }

    static line_t identity_line(I n) {
        line_t line(n);
        for (I i = 0; i < n; i++) {
            line[i] = i;
        }
        return line;
    }

    static line_t& expand_line(line_t& line, I n) {
        for (I i = (I)line.size(); i < n; i++) {
            line.push_back(i);
        }
        return line;
    }

    static line_t& apply_cycles_to_line(line_t& line, const cycles_t& cycles, I n = 0) {
        expand_line(line, n);
        for (const auto& c : cycles) {
            for (size_t i = 1; i < c.size(); i++) {
                std::swap(line[c[i - 1]], line[c[i]]);
            }
        }
        return line;
    }

    static line_t cycles_to_line(const cycles_t& cycles, I n) {
        auto line = identity_line(n);
        return apply_cycles_to_line(line, cycles, n);
    }

    static cycles_t line_to_cycles(const line_t& line) {
        I n = (I)line.size();
        cycles_t cycles;
        std::vector<bool> u(n);
        for (I i = 0; i < n; i++) {
            if (u[i]) continue;
            cycle_t c;
            for (I j = i; !u[j]; j = line[j]) {
                c.push_back(j);
                u[j] = true;
            }
            if (c.size() > 1) cycles.emplace_back(c);
        }
        return cycles;
    }

    static line_t& apply_transpositions_to_line(line_t& line, const transpositions_t& transpositions) {
        for (auto t : transpositions) {
            std::swap(line[t.first], line[t.second]);
        }
        return line;
    }

    static line_t transpositions_to_line(const transpositions_t& transpositions, I n) {
        auto line = identity_line(n);
        return apply_transpositions_to_line(line, transpositions);
    }

    static transpositions_t line_to_transpositions(const line_t& line) {
        return cycles_to_transpositions(line_to_cycles(line));
    }

    static transpositions_t cycles_to_transpositions(const cycles_t& cycles) {
        transpositions_t transpositions;
        for (const auto& c : cycles) {
            for (size_t i = 1; i < c.size(); i++) {
                transpositions.push_back({ c[i - 1], c[i] });
            }
        }
        return transpositions;
    }

    static cycles_t transpositions_to_cycles(const transpositions_t& transpositions) {
        return line_to_cycles(transpositions_to_line(transpositions, size(transpositions)));
    }

    static I size(const cycles_t& cycles) {
        I n = 0;
        for (const auto& c : cycles) {
            for (size_t i = 0; i < c.size(); i++) {
                n = std::max(n, c[i] + 1);
            }
        }
        return n;
    }

    static I size(const transpositions_t& transpositions) {
        I n = 0;
        for (auto t : transpositions) {
            n = std::max(n, t.first + 1);
            n = std::max(n, t.second + 1);
        }
        return n;
    }
};

template<typename I>
struct identityT<permutation<I>> {
    static permutation<I> of(const permutation<I>& p) {
        return permutation<I>(p.n);
    }
};

} // math
} // altruct
