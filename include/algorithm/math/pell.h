#pragma once

#include <type_traits>
#include <vector>

#include "algorithm/math/base.h"
#include "algorithm/math/primes.h"
#include "structure/math/quadratic.h"

namespace altruct {
namespace math {

/**
 * Based on the following paper:
 *   John P. Robertson - Solving the generalized Pell equation "x^2 - D y^2 = N" (2004)
 */

/**
 * PQa algorithm for Pell's equations.
 *
 * Computes continued fraction expansion of the quadratic irrational
 * `(P0 + sqrt(D)) / Q0` and some auxiliary variables:
 *
 * (P0 + sqrt(D)) / Q0 = a0 + 1/( a1 + 1/( a2 + 1/( a3 ...
 * (P0 + sqrt(D)) / Q0 = lim[Ai / Bi, i -> inf]
 * gcd(Ai, Bi) = 1
 * Gi = Q0 Ai - P0 Bi
 * Gi^2 - D Bi^2 = (-1)^(i+1) Q[i+1] Q0;  (used in the solutions)
 *
 * @param D - integer such that: D > 0, D not a square;
 * @param Q - Q terms of the expansion;
 *            integer Q[0] must be given initially: Q0 != 0;
 *            further terms get populated by the method;
 * @param P - P terms of the expansion;
 *            integer P[0] must be given initially: P0^2 = D (mod Q0)
 *            further terms get populated by the method;
 * @param a - the main coefficients of the expansion; see above;
 *            initially empty, terms get populated by the method;
 * @param A, B - Ai / Bi are the convergents to the continued fraction;
 *            initially empty, terms get populated by the method;
 * @param G - terms used in the solution of the pell equation;
 *            initially empty, terms get populated by the method;
 * @return  - the length of the period
 */
template<typename I>
int pell_PQa(I D, std::vector<I>& P, std::vector<I>& Q, std::vector<I>& a, std::vector<I>& A, std::vector<I>& B, std::vector<I>& G) {
    I qD = sqrtT<I>(D);
    int i = (int)a.size(); int i0 = i;
    I P0 = P[i];
    I Q0 = Q[i];
    I A2 = (i >= 2) ? A[i - 2] : 0; I A1 = (i >= 2) ? A[i - 1] : 1;
    I B2 = (i >= 2) ? B[i - 2] : 1; I B1 = (i >= 2) ? B[i - 1] : 0;
    I G2 = (i >= 2) ? G[i - 2] : -P0; I G1 = (i >= 2) ? G[i - 1] : +Q0;
    for (;;) {
        I a0 = div_floor<I>(P[i] + (Q[i] < 0 ? qD + 1 : qD), Q[i]); a.push_back(a0);
        I A0 = a0 * A1 + A2; A.push_back(A0); A2 = A1; A1 = A0;
        I B0 = a0 * B1 + B2; B.push_back(B0); B2 = B1; B1 = B0;
        I G0 = a0 * G1 + G2; G.push_back(G0); G2 = G1; G1 = G0;
        i++;
        P0 = a0 * Q0 - P0; P.push_back(P0);
        Q0 = (D - P0 * P0) / Q0; Q.push_back(Q0);
        if (i0 == 0 && Q0 > 0 && P0 > 0 && sqT<I>(P0) < D && sqT<I>(P0 - Q0) < D && D < sqT<I>(P0 + Q0)) i0 = i;
        if (i0 != 0 && i > i0 && P[i] == P[i0] && Q[i] == Q[i0]) return i - i0;
    }
}

/**
 * Solves the ordinary Pell's equation.
 *
 * Computes the minimal positive solution `(x0, y0)` that satisfies `x^2 - D y^2 = +/-1`.
 * All solutions to the `x^2 - D y^2 = +1` equation are given as `(x0 + y0 sqrt(d))^n`.
 * All solutions to the `x^2 - D y^2 = -1` equation are given as `(x0 + y0 sqrt(d))^n`, for odd n.
 *
 * @param D - integer such that: D > 0, D not a square;
 * @param N - integer that can only be +/-1
 * @param x0, y0 - the minimal positive solution
 * @return the length of the period if the solution exists, 0 otherwise;
 */
template<typename I>
int pell1(I D, I N, I& x0, I& y0) {
    std::vector<I> P{ 0 }, Q{ 1 }, a, A, B, G;
    int l = pell_PQa<I>(D, P, Q, a, A, B, G);
    if (l % 2 == 1) {
        // Length of the period is odd.
        // All solutions to the `x^2 - D y^2 = -1` equation are
        // given as `(G[kl-1], B[kl-1])` for positive odd `k`.
        // The minimal solution is `(G[l-1], B[l-1])`.
        if (N == -1) {
            x0 = G[l - 1], y0 = B[l - 1];
            return l;
        }
        // All solutions to the `x^2 - D y^2 = +1` equation are
        // given as `(G[kl-1], B[kl-1])` for positive even `k`.
        // The minimal solution is `(G[2l-1], B[2l-1])`.
        if (N == +1) {
            // calculate another period
            pell_PQa<I>(D, P, Q, a, A, B, G);
            x0 = G[2 * l - 1], y0 = B[2 * l - 1];
            return l;
        }
    } else {
        // Length of the period is odd.
        // The `x^2 - D y^2 = -1` equation has no solutions.
        // All solutions to the `x^2 - D y^2 = +1` equation are
        // given as `(G[kl-1], B[kl-1])` for positive `k`.
        // The minimal solution is `(G[l-1], B[l-1])`.
        if (N == +1) {
            x0 = G[l - 1], y0 = B[l - 1];
            return l;
        }
    }
    return 0;
}

/**
 * Solves the generalized Pell's equation when `1 < N^2 < D`.
 *
 * Computes the minimal positive solutions `(xc0, yc0)` that satisfy `x^2 - D y^2 = N`,
 * for each equivalence class.
 * For each equivalence class, all solutions to the `x^2 - D y^2 = N` equation
 * are given as `(xc0 + yc0 sqrt(D)) * (x0 + y0 sqrt(d))^n`,
 * where `(x0, y0)` is the minimal solution to the `x^2 - D y^2 = 1` equation.
 *
 * @param D - integer such that: D > 0, D not a square;
 * @param N - integer such that: 1 < N^2 < D
 * @param x0, y0 - the minimal positive solutions for each equivalence class
 */
template<typename I>
void pellS(I D, I N, std::vector<I>& xc0, std::vector<I>& yc0) {
    std::vector<I> P{ 0 }, Q{ 1 }, a, A, B, G;
    int l = pell_PQa<I>(D, P, Q, a, A, B, G);
    if (l % 2 == 1 || Q[l] != 1) {
        pell_PQa<I>(D, P, Q, a, A, B, G);
        l *= 2;
    }
    for (int i = 0; i < l; i++) {
        // m = Gi^2 - D Bi^2 = (-1)^(i+1)*Q[i+1]
        // if m == N / f^2 for some f > 0,
        // add (f*Gi,f*Bi) to the list of solutions
        if (N % Q[i + 1] != 0) continue;
        I f2 = N / Q[i + 1];
        if (i % 2 == 0) f2 = -f2;
        if (f2 < 0) continue;
        I f = sqrtT<I>(f2);
        if (sqT<I>(f) != f2) continue;
        xc0.push_back(f * G[i]);
        yc0.push_back(f * B[i]);
    }
}

/**
 * Solves the generalized Pell's equation.
 *
 * Computes the minimal positive solutions `(xc0, yc0)` that satisfy `x^2 - D y^2 = N`,
 * for each equivalence class.
 * For each equivalence class, all solutions to the `x^2 - D y^2 = N` equation
 * are given as `(xc0 + yc0 sqrt(D)) * (x0 + y0 sqrt(d))^n`,
 * where `(x0, y0)` is the minimal solution to the `x^2 - D y^2 = 1` equation.
 *
 * @param D - integer such that: D > 0, D not a square;
 * @param N - integer such that: N != 0
 * @param fN - prime factorization of N
 * @param x0, y0 - the minimal positive solutions for each equivalence class
 */
template<typename I, typename PF>
void pell(I D, I N, const std::vector<std::pair<PF, int>>& fN, std::vector<I>& xc0, std::vector<I>& yc0) {
    std::vector<std::pair<PF, int>> fN2;
    for (const auto& f : fN) {
        if (f.second >= 2) fN2.push_back({ f.first, f.second / 2 });
    }
    std::vector<I> vf; divisors(vf, fN2);
    for (const auto& f : vf) {
        I m = N / sqT<I>(f);
        I mm = absT<I>(m);
        int ms = (m < 0) ? -1 : +1;
        I mh = mm / 2;
        // find all z such that z^2 = D (mod |m|)
        // TODO: this is fine for small N, but for larger N
        // we should use modular square root algorithms.
        std::vector<I> vz;
        for (I z = 0; z < mm; z += 1) {
            if ((z * z - D) % mm == 0) vz.push_back(z);
        }
        for (const auto& z : vz) {
            I zh = (z <= mh) ? z : z - mm;
            std::vector<I> P{ zh }, Q{ mm }, a, A, B, G;
            int l = pell_PQa<I>(D, P, Q, a, A, B, G);
            int i = 1; while (i < (int)Q.size() && Q[i] != +1 && Q[i] != -1) i++;
            if (i >= (int)Q.size()) continue;
            // check for solution by testing Q[i] == (-1)^i sign(m)
            if (Q[i] != ((i % 2 == 1) ? -ms : +ms)) {
                // G[i-1]^2 - D B[i-1]^2 = (-1)^i Q[i] |m| = -m;
                // do another period to get +m
                pell_PQa<I>(D, P, Q, a, A, B, G);
                i += l;
            }
            if (Q[i] != ((i % 2 == 1) ? -ms : +ms)) {
                // G[i-1]^2 - D B[i-1]^2 = (-1)^i Q[i] |m| = -m;
                // -m again, +m will never happen
            } else {
                // G[i-1]^2 - D B[i-1]^2 = (-1)^i Q[i] |m| = +m;
                xc0.push_back(f * G[i - 1]);
                yc0.push_back(f * B[i - 1]);
            }
        }
    }
}

/**
 * Lists all the primitive solutions to the ordinary Pell's equation, up to the specified limit.
 *
 * @param D - integer such that: D > 0, D not a square;
 * @param N - can only be +/-1
 * @param x_max - only solutions with x <= x_max; 0 means unbounded
 * @param y_max - only solutions with y <= y_max; 0 means unbounded
 * @param count - only the first `count` solutions; 0 means unbounded
 * @return the solutions
 */
template<typename I>
std::vector<quadraticX<I>> pell1(I D, I N, I x_max, I y_max, int count) {
    std::vector<quadraticX<I>> vs;
    quadraticX<I> s0(0, 0, D);
    if (pell1<I>(D, N, s0.a, s0.b)) {
        auto s = s0;
        if (N == -1) s0 *= s0;
        while (true) {
            auto t = s; s *= s0;
            if (t.a < 0) t.a = -t.a;
            if (t.b < 0) t.b = -t.b;
            if (x_max > 0 && t.a > x_max) break;
            if (y_max > 0 && t.b > y_max) break;
            if (count > 0 && (int)vs.size() >= count) break;
            vs.push_back(t);
        }
    }
    return vs;
}

/**
 * Lists all the primitive solutions to the generalized Pell's equation, up to the specified limit.
 *
 * @param D - integer such that: D > 0, D not a square;
 * @param N - integer such that: N != 0
 * @param x_max - only solutions with x <= x_max; 0 means unbounded
 * @param y_max - only solutions with y <= y_max; 0 means unbounded
 * @param count - only the first `count` solutions; 0 means unbounded
 * @return the solutions
 */
template<typename I, typename P, typename F>
void pell(I D, I N, const std::vector<std::pair<P, int>>& fN, I x_max, I y_max, int count, F visitor) {
    int size = 0;
    quadraticX<I> s0(0, 0, D);
    if (pell1<I>(D, 1, s0.a, s0.b)) {
        std::vector<I> vxc0, vyc0;
        pell<I>(D, N, fN, vxc0, vyc0);
        std::vector<quadraticX<I>> vsc(vxc0.size());
        for (int c = 0; c < (int)vsc.size(); c++) {
            vsc[c] = { vxc0[c], vyc0[c], D };
        }
        while (true) {
            auto size0 = size;
            for (auto& s : vsc) {
                auto t = s; s *= s0;
                if (t.a < 0) t.a = -t.a;
                if (t.b < 0) t.b = -t.b;
                if (x_max > 0 && t.a > x_max) continue;
                if (y_max > 0 && t.b > y_max) continue;
                if (count > 0 && size >= count) continue;
                visitor(t); size++;
            }
            if (size == size0) break;
        }
    }
}
template<typename I, typename P>
std::vector<quadraticX<I>> pell(I D, I N, const std::vector<std::pair<P, int>>& fN, I x_max, I y_max, int count) {
    std::vector<quadraticX<I>> vs;
    pell<I>(D, N, fN, x_max, y_max, count, [&](const quadraticX<I>& s) { vs.push_back(s); });
    return vs;
}

} // math
} // altruct
