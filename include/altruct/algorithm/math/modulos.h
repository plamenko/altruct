#pragma once

#include "altruct/structure/math/modulo.h"
#include "altruct/structure/math/quadratic.h"
#include "altruct/structure/math/prime_holder.h"

#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>

namespace altruct {
namespace math {

/**
 * Chinese Remainder
 *
 * Calculates `a` and `n` so that:
 *   n = lcm(n1, n2)
 *   a % n1 == a1
 *   a % n2 == a2
 *   0 <= a < n
 * `n1` and `n2` don't have to be coprime.
 * Correctly handles 64bit result for long long type.
 */
template<typename T>
void chinese_remainder(T* a, T* n, T a1, T n1, T a2, T n2) {
    T e0 = zeroT<T>::of(*a);
    T ni1, ni2; T g = modulo_gcd_ex(n1, n2, ni1, ni2);
    if ((a2 - a1) % g != e0) { *n = e0, *a = e0; return; }
    T t1 = modulo_mul(a1, ni2, n1);
    T t2 = modulo_mul(a2, ni1, n2);
    n1 /= g; n2 /= g; *n = n1 * n2 * g;
    *a = modulo_add(modulo_mul(t1, n2, *n), modulo_mul(t2, n1, *n), *n);
}
template<typename T>
void chinese_remainder(T* ar, T* nr, T a, T n) {
    chinese_remainder(ar, nr, *ar, *nr, a, n);
}
template<typename T>
T chinese_remainder(T a1, T n1, T a2, T n2) {
    T a, n; chinese_remainder(&a, &n, a1, n1, a2, n2);
    return a;
}

/**
 * Calculates coefficients of the mixed radix equation using Garner Algorithm.
 *
 * Let `u` be given as a set of congruence equations: `u == a[i] (mod p[i])`.
 * Let `q[i]` be the product of the first `i` elements of `p`. `q[0] = 1`.
 * Then `a` can be represented in the mixed radix form: `u = Sum[x[i] * q[i]]`.
 * Note: `u` is unique modulo `q[n]` as per Chinese Remainder Theorem.
 *
 * @param vap - vector of (remainder, modulus) pairs, moduli should be pairwise relatively prime
 * @return vx - vector of (coeeficient, modulus) pairs, moduli are same as in `vap`
 */
template<typename V>
V garner(const V& vap) {
    V vx(vap.size());
    for (int i = 0; i < vx.size(); i++) {
        auto y = vap[i];
        for (int j = 0; j < i; j++) {
            y -= vx[j];
            y /= vap[j].M();
        }
        vx[i] = y;
    }
    return vx;
}

/**
 * Jacobi symbol
 *
 * For prime `m`, this is equivalent to Legendre symbol:
 *    0 - if `n` is `0` mod `m`
 *   +1 - if `n` is a quadratic residue mod `m`
 *   -1 - if `n` is a quadratic nonresidue mod `m`
 *
 * For composite `m`, result of `+1` only means that `n` is
 * a quadratic nonresidue for an even number (zero or more)
 * of prime factors of `m`. In order for `n` to be a quadratic
 * residue, it has to be a residue for each prime factor of `m`.
 */
template <typename I>
int jacobi(I n, I m) {
    int e, j = 1;
    for (;;) {
        if (m == 1) return j;
        n %= m;
        if (n == 0) return 0;
        for (e = 0; (n % 2) == 0; e++) n /= 2;
        if ((e % 2) == 1 && (m % 8) == 3) j = -j;
        if ((e % 2) == 1 && (m % 8) == 5) j = -j;
        if ((n % 4) == 3 && (m % 4) == 3) j = -j;
        std::swap(n, m);
    }
}

/**
 * Square root of `y.v` modulo prime `y.M`
 *
 * @param M - the modulo<I, ...> type
 * @param y - integer such that: y = x^2 (mod p), (y, p) = 1
 */
template <typename M>
M sqrt_cipolla(const M& y) {
    if (y.M() == 2) return y;
    M e0 = zeroOf(y), e1 = identityOf(y);
    // find a quadratic nonresidue `d` modulo `p`
    M a = e0, d = e0;
    do {
        a += e1, d = a * a - y;
    } while (powT(d, (y.M() - 1) / 2) == 1); // jacobi(d, p) == 1
    // r = (a + sqrt(d)) ^ ((p + 1) / 2)
    return powT(quadraticX<M>(a, e1, d), (y.M() + 1) / 2).a;
}

/**
 * Square root of `y` modulo prime `p`
 *
 * @param y - integer such that: y = x^2 (mod p), (y, p) = 1
 */
template <typename I>
I sqrt_cipolla(const I& y, const I& p) {
    return sqrt_cipolla(moduloX<I>(y, p)).v;
}

/**
 * Square roots of `y` modulo `2^k`.
 * All solutions are given by {x1, x2, -x2, -x1}.
 *
 * @param y - integer such that: y = x^2 (mod 2^k), (y, 2) = 1
 * @return {x1, x2}, or {0, 0} if there is no solution
 */
template <typename I>
std::pair<I, I> sqrt_hensel_lift_p2(const I& y, I k) {
    typedef moduloX<I> modx;
    if (y % 2 != 1) return { 0, 0 };
    if (k == 1) return { 1, 1 };
    if (y % 4 != 1) return { 0, 0 };
    if (k == 2) return { 1, 1 };
    if (y % 8 != 1) return { 0, 0 };
    I s[2] = { 1, 3 };
    I w2 = 4;
    for (int i = 4; i <= k; i++) {
        for (int j = 0; j < 2; j++) {
            modx r = modx(s[j], w2 * 4);
            modx v = r * r - modx(y, r.M());
            if (v.v != 0) {
                s[j] += w2;
            }
            else if (s[j] >= w2) {
                s[j] += w2 * 2;
            }
        }
        w2 *= 2;
    }
    return { s[0], s[1] };
}

/**
 * Square root of `y` modulo odd prime power `p^k`
 * All solutions are given by {x, -x}.
 *
 * @param y - integer such that: y = x^2 (mod p^k), (y, p) = 1
 */
template <typename I>
I sqrt_hensel_lift(const I& y, const I& p, I k) {
    typedef moduloX<I> modx;
    // f(r) == r^2 - y; f'(r) == 2r;
    modx r = sqrt_cipolla(modx(y, p));
    for (I i = 1; i < k; i *= 2) {
        I phi = r.M() / p * (p - 1); // euler_phi(r.M)
        modx u = powT(r * 2, phi - 1); // f'(r) ^-1
        r.M() = (i * 2 < k) ? r.M() * r.M() : powT(p, k); // lift modulus
        modx v = r * r - modx(y, r.M()); // f(r)
        r -= v * u;
    }
    return r.v;
}

/**
 * Square roots of `y` modulo `m`
 *
 * @param y - integer such that: y = x^2 (mod m), (y, m) = 1
 * @param vf - factorization of m
 * @return - all square roots
 */
template <typename I, typename P>
std::vector<I> sqrt_mod(I y, const std::vector<std::pair<P, int>>& vf) {
    I m = 1;
    std::unordered_set<I> sr{ 0 };
    auto add = [&](std::unordered_set<I>& sr0, I x, I q) {
        if ((x * x - y) % q != 0) return;
        for (I r0 : sr0) {
            sr.insert(chinese_remainder(r0, m, x, q));
            sr.insert(chinese_remainder(r0, m, modulo_neg(x, q), q));
        }
    };
    for (const auto& f : vf) {
        std::unordered_set<I> sr0;
        sr0.swap(sr);
        I q = powT<I>(f.first, f.second);
        if (f.first == 2) {
            auto x = sqrt_hensel_lift_p2(y, f.second);
            add(sr0, x.first, q);
            add(sr0, x.second, q);
        }
        else {
            auto x = sqrt_hensel_lift(y, f.first, f.second);
            add(sr0, x, q);
        }
        m *= q;
    }
    return std::vector<I>(sr.begin(), sr.end());
}

/**
 * Multiplicative order of `a` in a cyclic group of order `n`
 *
 * @param n
 * @param n_factors - unique prime factors of `n`
 * @return - smallest `k` such that `a^k = a^0`, `k` divides `n`
 */
template<typename G, typename I, typename P>
I multiplicative_order_g(G a, I n, const std::vector<P>& n_factors) {
    G id = identityOf(a);
    I k = n;
    for (const P& p : n_factors) {
        while (k % p == 0 && powT(a, I(k / p)) == id) k /= p;
    }
    return k;
}

/**
 * Multiplicative order of x modulo `m`
 *
 * @param m - modulus
 * @param phi - `euler_phi(m)`; number of coprimes with `m` up to `m`
 * @param phi_factors - unique prime factors of `phi`
 * @return - smallest `k` such that `a^k = 1 (mod m)`
 */
template<typename I, typename P>
I multiplicative_order(I a, I m, I phi, const std::vector<P>& phi_factors) {
    typedef moduloX<I> modx;
    return multiplicative_order_g(modx(a, m), phi, phi_factors);
}

/**
 * Primitive root modulo `m`
 *
 * `m` must be 2, 4, p^k or 2p^k.
 *
 * A root `g` generates all the elements coprime to `m`
 * as `g^k` for `0 <= k < phi`
 *
 * @param m - modulus
 * @param phi - `euler_phi(m)`; number of coprimes with `m` up to `m`
 * @param phi_factors - unique prime factors of `phi`
 * @return - root or 0 if there is none
 */
template<typename I, typename P>
I primitive_root(I m, I phi, const std::vector<P> &phi_factors) {
    typedef moduloX<I> modx;
    for (I g = 1; g < m; g += 1) {
        if (gcd(g, m) > 1) continue;
        bool primitive = true;
        for (const P& p : phi_factors) {
            primitive &= (powT(modx(g, m), I(phi / p)) != 1);
        }
        if (primitive) return g;
    }
    return 0;
}

/**
 * Primitive root of unity modulo `m`
 *
 * `r` with maximal `k` where `r^k=1 (mod m)`
 * 
 * While a primitive root that generates all elements
 * coprime to `m` might not exist, `r` always exists.
 * When primitive root does exist, `r = g` holds.
 *
 * @param m - modulus
 * @param lam - `carmichael_lambda(m)`;
 * @param lam_factors - unique prime factors of `lam`
 * @return - root
 */
template<typename I>
I primitive_root_of_unity(I m, I lam, const std::vector<I> &lam_factors) {
    return primitive_root(m, lam, lam_factors);
}

/**
 * Primitive root modulo `m`
 *
 * `m` must be 2, 4, p^k or 2p^k.
 */
int primitive_root(int m, prime_holder& prim);

/**
 * Primitive root of unity modulo `m`
 */
int primitive_root_of_unity(int m, prime_holder& prim);

/**
 * `k-th` roots of unity modulo `m`
 *
 * `m` must be 2, 4, p^k or 2p^k.
 *
 * @param m - modulus
 * @param k - k-th root
 * @param lam - `carmichael_lambda(m)`
 * @param g - primitive root of unity modulo `m`
 */
template<typename I>
std::set<I> kth_roots_of_unity(I m, I k, I lam, I g) {
    typedef moduloX<I> modx;
    I d = gcd(k, lam);
    modx w = powT(modx(g, m), I(lam / d));
    modx r = identityT<modx>::of(w);
    std::set<I> sr;
    for (I j = 0; j < d; j++) {
        sr.insert(r.v);
        r *= w;
    }
    return sr;
}

/**
 * `k-th` roots of unity modulo `m`
 */
std::set<int> kth_roots_of_unity(int m, int k, prime_holder& prim);

/**
* `k-th` roots of `n` modulo `m`
*
* `m` must be 2, 4, p^k or 2p^k.
*
* @param m - modulus
* @param k - k-th root
* @param phi - `euler_phi(m)`
* @param g - primitive root modulo `m`
* @param l - `l` such that `g ^ l == n  (mod m)`
*/
template<typename I>
std::set<I> kth_roots(I m, I k, I phi, I g, I l) {
    typedef moduloX<I> modx;
    I d = gcd(k, phi);
    if (d == 0 || l % d != 0) return{};
    phi /= d; l /= d; k /= d;
    // g^(l/k) == n^(1/k)
    modx h = modx(l, phi) / k;
    modx r = powT(modx(g, m), h.v);
    modx w = powT(modx(g, m), phi);
    std::set<I> sr;
    for (int i = 0; i < d; i++) {
        sr.insert(r.v);
        r *= w;
    }
    return sr;
}

/**
 * Discrete logarithm in group; brute-force in O(max_iter)
 *
 * `a^x = b`
 *
 * @param max_iter - maximum number of iterations
 * @return smallest `x`
 */
template <typename G, typename I>
I discrete_log_brute_force_g(G a, G b, I max_iter) {
    G a_x = identityOf(a);
    for (I x = 0; x < max_iter; x++) {
        if (a_x == b) return x;
        a_x *= a;
    }
    return -1;
}

/**
 * Discrete logarithm in cyclic group of order `n`; in O(sqrt(n))
 *
 * `a^x = b`, order of `a` is `n` (or its divisor)
 *
 * @return `x`
 */
template <typename G, typename I>
I discrete_log_baby_giant_g(G a, G b, I n) {
    I q = sqrtT(n) + 1;
    std::unordered_map<G, I, hasherT<G>> baby;
    // baby steps: a^j
    G cur = identityOf(a);
    for (I j = 0; j < q; ++j) {
        if (!baby.emplace(cur, j).second) break;
        cur *= a;
    }
    G alpha = powT(a, n - q); // a^-q
    // giant steps: b * alpha^i
    G gamma = b;
    for (I i = 0; i < q; ++i) {
        if (baby.count(gamma)) {
            return i * q + baby[gamma];
        }
        gamma *= alpha;
    }
    return -1;
}

/**
 * Discrete logarithm in cyclic group of order `p^s`, in O(s*sqrt(p))
 *
 * `a^x = b`, order of `a` is `p^s` (or its divisor) with `p` prime
 *
 * @return `x`
 */
template <typename G, typename I>
I discrete_log_order_pp_g(G a, G b, I p, int s) {
    G id = identityOf(a);
    if (s == 0) return 0;
    if (a == id) return 0;
    if (b == id) return 0;
    if (a == b) return 1;
    // Pohlig�Hellman algorithm
    std::vector<I> pp(s + 1);
    pp[0] = 1;
    for (int k = 1; k <= s; k++) {
        // check if the order is a proper divisor of p^s
        if (powT(a, pp[k - 1]) == id) { s = k - 1; break; }
        pp[k] = pp[k - 1] * p;
    }
    G ai = powT(a, pp[s] - 1); // a^-1
    I x = 0;
    G gamma = powT(a, pp[s - 1]);
    for (int k = 0; k < s; k++) {
        G bk = powT(powT(ai, x) * b, pp[s - 1 - k]);
        I d = discrete_log_baby_giant_g(gamma, bk, p);
        x += pp[k] * d;
    }
    return x;
}

/**
 * Discrete logarithm in cyclic group of order `n`, in O(Sum[s_i*sqrt(p_i)])
 * Where:
 *   `n = Product[p_i^s_i]`
 *   `a^x = b`
 *
 * @param n
 * @param n_factors - unique prime factors of `n`
 * @param out_o - optional output `o`, order of `a`, `o` divides `n`
 * @return `x`
 */
template <typename G, typename I, typename P>
I discrete_log_g(G a, G b, I n, const std::vector<P>& n_factors, I* out_o = nullptr) {
    I o = multiplicative_order_g(a, n, n_factors);
    I x = 0, q = 1;
    for (const P& p : n_factors) {
        I d = o, pe = 1; int e = 0;
        while (d % p == 0) d /= p, pe *= p, e++;
        if (e == 0) continue;
        G a_sub = powT(a, d); // generator in subgroup
        G b_sub = powT(b, d); // target in subgroup
        I x_pe = discrete_log_order_pp_g(a_sub, b_sub, p, e);
        chinese_remainder(&x, &q, x_pe, pe);
    }
    if (out_o) *out_o = o; // q == o at this point
    return x;
}

/**
 * Discrete logarithm modulo `m`; brute-force in O(m)
 *
 * `a^x = b (mod m)`
 *
 * @return smallest `x`
 */
template <typename I>
I discrete_log_brute_force(I a, I b, I m) {
    typedef moduloX<I> modx;
    return discrete_log_brute_force_g(modx(a, m), modx(b, m), m);
}

/**
 * Discrete logarithm modulo `m`, in O(sqrt(n))
 *
 * `a^x = b (mod m)`, `a` and `b` comprime to `m`,
 * order of `a` mod `m` is `n` (or its divisor)
 *
 * @return `x`
 */
template <typename I>
I discrete_log_baby_giant(I a, I b, I m, I n) {
    typedef moduloX<I> modx;
    return discrete_log_baby_giant_g(modx(a, m), modx(b, m), n);
}

/**
 * Discrete logarithm modulo odd prime `p`, in O(sqrt(p))
 *
 * `a^x = b (mod p)`, `a` and `b` comprime to `p`
 *
 * @return `x`
 */
template <typename I>
I discrete_log_shanks(I a, I b, I p) {
    return discrete_log_baby_giant(a, b, p, p - 1);
}

/**
 * Discrete logarithm modulo odd prime `p`, in O(sqrt(p))
 *
 * `a^x = b (mod p)`, `a` and `b` comprime to `p`
 *
 * @return `x`
 */
template <typename I>
I discrete_log_oddp(I a, I b, I p) {
    // brute force for small p
    if (p < 1100) return discrete_log_brute_force(a, b, p);
    return discrete_log_shanks(a, b, p);
}

/**
 * Discrete logarithm modulo `2^s`, in O(s)
 *
 * `a^x = b (mod 2^s)`, `a` and `b` odd
 *
 * @return `x`
 */
template <typename I>
I discrete_log_p2(I a, I b, int s) {
    I m = powT<I>(2, s);
    I o = 1;
    I ao = a;
    int t = 0;
    while (ao != 1) {
        ao = modulo_mul(ao, ao, m);
        o *= 2;
        t++;
    }
    I ai = modulo_inv(a, m);
    I x = 0;
    for (int i = 0; i < t; i++) {
        I e = o >> (i + 1);
        I y = modulo_power(modulo_mul(b, modulo_power(ai, x, m), m), e, m);
        if (y != 1) x += I(1) << i;
    }
    return x;
}

/**
 * Discrete logarithm modulo prime power `p^s`, in O(s*sqrt(p))
 *
 * `a^x = b (mod p^s)`, `a` and `b` comprime to `p`
 *
 * Note: this is not necessarily the smallest such `x`!
 * Smallest one is `x % o` where `o = multiplicative_order(a, p^s)`
 *
 * @return `x`
 */
template <typename I>
I discrete_log_pp(I a, I b, I p, int s) {
    if (s == 0) return 0;
    if (a == 1) return 0;
    if (b == 1) return 0;
    if (a == b) return 1;
    if (p == 2) return discrete_log_p2(a, b, s);
    I x1 = discrete_log_oddp(a % p, b % p, p);
    if (s == 1) return x1;
    // Eric Bach - Discrete Logarithms and Factoring
    I p1 = p - 1;
    I phi = p1 * powT(p, s - 1);
    I ps1 = powT(p, s - 1), ps = ps1 * p;
    auto theta = [&](I c) { return (modulo_power(c, phi, ps * ps1) - 1) / ps; };
    I ta = theta(a), tb = theta(b);
    while (ta % p == 0 && tb % p == 0 && ps1 % p == 0) ta /= p, tb /= p, ps1 /= p;
    I xs1 = modulo_mul(tb, modulo_inv(ta, ps1), ps1);
    return chinese_remainder(x1, p1, xs1, ps1);
}

/**
 * Discrete logarithm modulo `m`, in O(Sum[s_i*sqrt(p_i)])
 * Where:
 *   `phi = Product[p_i^s_i]`
 *   `a^x = b (mod m)`, `a` and `b` comprime to `m`
 *
 * @param m - modulus
 * @param phi - `euler_phi(m)`; number of coprimes with `m` up to `m`
 * @param phi_factors - unique prime factors of `phi`
 * @param out_o - optional output `o`, order of `a` mod `m`.
 * @return `x`
 */
template <typename I, typename P>
I discrete_log(I a, I b, I m, I phi, const std::vector<P>& phi_factors, I* out_o = nullptr) {
    typedef moduloX<I> modx;
    return discrete_log_g(modx(a, m), modx(b, m), phi, phi_factors, out_o);
}

} // math
} // altruct
