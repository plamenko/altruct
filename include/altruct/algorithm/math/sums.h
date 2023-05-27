#pragma once

#include "altruct/algorithm/math/base.h"
#include "altruct/algorithm/math/recurrence.h"
#include "altruct/structure/math/polynom.h"

namespace altruct {
namespace math {

/**
 * Calculates `Sum[Floor[(a * k + b) / q], {k, 0, n - 1}]` in `O(log min(q, n))`.
 *
 * Note: `q` must not be zero.
 */
template<typename T, typename I>
T sum_ratio(I a, I b, I q, I n, T zero = T(0)) {
    auto TT = [&](I n) { return castOf(zero, n); };
    if (n < 1) return zero;
    if (q < 0) return sum_ratio<T, I>(-a, -b, -q, n, zero);
    if (a < 0) return -sum_ratio<T, I>(-a, q - 1 - b, q, n, zero);
    if (b < 0) { I m = -b / q + 1; return sum_ratio<T, I>(a, b + q * m, q, n, zero) - TT(n) * TT(m); }
    T s = zero;
    I i = 0;
    for (i = 0; n > 0; i++) {
        I n1 = n - 1;
        T t = (n % 2 == 1) ? TT(n) * TT(n1 / 2) : TT(n / 2) * TT(n1);
        s += TT(b / q) * TT(n) + TT(a / q) * t;
        b %= q, a %= q; if (a == 0) break;
        n = (a * n1 + b) / q;
        b = (q - 1) - b, std::swap(a, q);
        s += TT(n) * TT(n1);
        s = -s;
    }
    return (i % 2 == 1) ? -s : s;
}

/**
 * Calculates `Sum[f[k], {k, a, b}]` in `O(b - a)`.
 *
 * @param f - `f[n]`
 */
template<typename T, typename I, typename F>
T sum(F f, I a, I b, T zero = T(0)) {
    T r = zero;
    for (I k = b; k >= a; k--) {
        r += f(k);
    }
    return r;
}

/**
 * Calculates sum of powers: `Sum[k^p, {k, 1, n}]` in `O(p)`.
 *
 * Note: for powers bigger than 3, T must be a field.
 * Example fields that are suitable: double, fraction, modulo.
 *
 * @param B - vector of Bernoulli numbers of the second kind up to `p`.
 */
template<typename T, typename I>
T sum_pow(int p, I N, const std::vector<T>& B) {
    T e0 = zeroOf(B[0]), e1 = identityOf(B[0]);
    T n = castOf(e1, N);
    if (p == 0) return n;
    if (p == 1) return n * (n + e1) / 2;
    if (p == 2) return n * (n + e1) * (n * 2 + e1) / 6;
    if (p == 3) return sqT<T>(n * (n + e1) / 2);
    //Faulhaber's formula
    T r = e0;
    T n_k = e1;
    T bin = e1;
    for (int k = 0; k <= p; k++) {
        n_k *= n;
        bin *= p - k + 1;
        bin /= k + 1;
        r += bin * B[p - k] * n_k;
    }
    return r / (p + 1);
}

/**
 * Calculates `Sum[k^p, {k, 1, n}]` in `O(p^2)`.
 *
 * Note: for powers bigger than 3, T must be a field.
 * Example fields that are suitable: double, fraction, modulo.
 */
template<typename T, typename I>
T sum_pow(int p, I n, T id = T(1)) {
    static std::vector<T> B;
    if (B.size() < p + 1) {
        int sz = (int)(B.size() + B.size() / 2);
        B = bernoulli_b<T>(std::max(p, sz), id);
    }
    return sum_pow(p, n, B);
}

/**
 * Calculates `Sum[k^m x^k, {k, 1, n}]` in `O(m^2)`.
 *
 * Note: `x != 1` must hold. For `x == 1` use `sum_pow`.
 */
template<typename T, typename I>
T sum_powx(int m, T x, I n) {
    T T0 = zeroT<T>::of(x), T1 = identityT<T>::of(x);
    T Tn = castOf(T1, n);
    polynom<T> p{ T0, T1 }, q{ T1 }, z{ T0, -T1, T1 };
    for (int k = 1; k <= m; k++) {
        auto Tk = castOf(T1, k);
        p = z * p.derivative() - (polynom<T>{Tn, Tk - Tn}) * p;
        q = z * q.derivative() - (polynom<T>{T0, Tk - T0}) * q;
    }
    return (powT(x, n) * p(x) - q(x)) / (powT(x - T1, m + 1)) - (m == 0 ? T1 : T0);
}

/**
 * Calculates `Sum[f[n/k], {k, 1, n}]` in `O(sqrt(n))`.
 *
 * @param f - `f[n]`
 */
template<typename T, typename I, typename F>
T sum_sqrt(F f, I n, T zero = T(0)) {
    if (n < 1) return zero;
    I q = sqrtT(n);
    T r = zero;
    for (I k = 1; k <= n / q; k++) {
        r += f(n / k);
    }
    for (I m = 1; m < q; m++) {
        r += f(m) * castOf(zero, (n / m) - (n / (m + 1)));
    }
    return r;
}

/**
 * Calculates `Sum[f[k] * g[n/k], {k, 1, n}]` in `O(sqrt(n))`.
 *
 * @param sf - `sf[n] = Sum[f[k], {k, 1, n}]`
 */
template<typename T, typename I, typename F1, typename F2>
T sum_sqrt2m(F1 sf, F2 g, I n, T zero = T(0)) {
    if (n < 1) return zero;
    I q = sqrtT(n);
    T r = zero;
    T sf0 = sf(n);
    for (I k = 1; k <= n / q; k++) {
        r += (sf(k) - sf(k - 1)) * g(n / k); // f(k) * g(n / k)
    }
    for (I m = 1; m < q; m++) {
        T sf1 = sf(n / (m + 1));
        r += (sf0 - sf1) * g(m);
        sf0 = sf1;
    }
    return r;
}

/**
 * Calculates `Sum[f[k] * g[n/k], {k, 1, n}]` in `O(sqrt(n))`.
 *
 * @param sf - `sf[n] = Sum[f[k], {k, 1, n}]`
 */
template<typename T, typename I, typename F1, typename F2, typename F3>
T sum_sqrt2m(F1 f, F2 sf, F3 g, I n, T zero = T(0)) {
    if (n < 1) return zero;
    I q = sqrtT(n);
    T r = zero;
    T sf0 = sf(n);
    for (I k = 1; k <= n / q; k++) {
        r += f(k) * g(n / k);
    }
    for (I m = 1; m < q; m++) {
        T sf1 = sf(n / (m + 1));
        r += (sf0 - sf1) * g(m);
        sf0 = sf1;
    }
    return r;
}

/**
 * Calculates `Sum[f[k, n/k], {k, 1, n}]` in `O(sqrt(n))`.
 *
 * @param sf - `sf[n, m] = Sum[f[k, m], {k, 1, n}]`
 */
template<typename T, typename I, typename F>
T sum_sqrt2(F sf, I n, T zero = T(0)) {
    if (n < 1) return zero;
    I q = sqrtT(n);
    T r = zero;
    T sf0 = sf(n, 1);
    for (I k = 1; k <= n / q; k++) {
        r += sf(k, n / k) - sf(k - 1, n / k); // f(k, n / k)
    }
    for (I m = 1; m < q; m++) {
        r += sf(n / m, m) - sf(n / (m + 1), m);
    }
    return r;
}

} // math
} // altruct
