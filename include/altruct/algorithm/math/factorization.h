#pragma once

#include "base.h"
#include "altruct/structure/math/modulo.h"

#include <vector>
#include <string>

namespace altruct {
namespace math {

/**
 * Miller-Rabin primality test.
 *
 * Probabilistic primality test with accuracy `4^-k`, where `k` is the number of bases tested.
 *
 * Complexity: O(p^(1/2)) <= O(n^(1/4)), where `p` is the smallest prime factor of `n`.
 *
 * @param n - number to test for primality
 * @param bases - a null-terminated array of bases to test against.
 * @return - true means `n` is probably prime, false means `n` is certainly composite.
 */
template<typename T>
bool miller_rabin(const T& n, const T* bases) {
    if (n == 0 || n == 1) return 0;
    if (n == 2 || n == 3) return 1;
    if ((n % 2) == 0) return 0;
    typedef moduloX<T> modx;
    T d = n - 1; int r = 0; // n-1 = 2^r * d
    while (d % 2 == 0) d /= 2, r++;
    for (int i = 0; bases[i] && bases[i] < n; i++) {
        modx x = powT(modx(bases[i], n), d);
        if (x == 1 || x == n - 1) continue;
        for (int i = 1; i < r; i++) {
            x *= x;
            if (x == 1 || x == n - 1) break;
        }
        if (x != n - 1) return false; // composite
    }
    return true; // probably prime
}

/**
 * Miller-Rabin primality test.
 *
 * Selects the appropriate bases based on the input size so that the test is
 * deterministic. See `miller_rabin` above.
 */
template<typename T>
bool miller_rabin(const T& n) {
    // 10^3, 2^10
    static T bases1[] = { 2, 0 };
    if (n < 2047) return miller_rabin(n, bases1);
    // 10^6, 2^23
    static T bases2[] = { 31, 73, 0 };
    if (n < 9080191) return miller_rabin(n, bases2);
    // 10^9, 2^32
    static T bases3[] = { 2, 7, 61, 0 };
    if (n < 4759123141LL) return miller_rabin(n, bases3);
    // 10^12, 2^40
    static T bases4[] = { 2, 13, 23, 1662803, 0 };
    if (n < 1122004669633LL) return miller_rabin(n, bases4);
    // 10^15, 2^48
    static T bases7[] = { 2, 3, 5, 7, 11, 13, 17, 0 };
    if (n < 341550071728321LL) return miller_rabin(n, bases7);
    // 10^18, 2^61
    static T bases9[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 0 };
    if (n < 3825123056546413051LL) return miller_rabin(n, bases9);
    // fallback to bases9 for larger numbers too
    return miller_rabin(n, bases9);
}

/**
 * Pollard's Rho factorization algorithm.
 *
 * Attempts to find a non-trivial, not necessarily prime factor of `n`.
 *
 * Before running this algorithm one should make sure that `n` is not a prime.
 * See the `miller_rabin` primality test above.
 *
 * In case the factorization fails (returned value is same as the input value `n`),
 * one should try with different `k` and `a`. See `pollard_rho_repeated` below.
 *
 * Complexity: O(p^(1/2)) <= O(n^(1/4)), where `p` is the smallest prime factor of `n`.
 *
 * @param n - number to factor
 * @param k - initial value
 * @param a - parameter of the polynomial g(x) = x^2 + a
 * @param max_inner_iter - maximum allowed number of iterations
 * @return d - a nontrivial factor of `n`, or `n` if factorization failed
 */
template<typename I>
I pollard_rho(const I& n, const I& k = 2, const I& a = 1, I max_inner_iter = 1000000) {
    if (n == 0) return 0;
    if (n == 1) return 1;
    if (n % 2 == 0) return 2;
    typedef moduloX<I> modx;
    auto g = [a](const modx& x){ return x*x + a; };
    modx x = modx(k, n), y = modx(k, n); I d = 1;
    while (d == 1 && max_inner_iter-- > 0) {
        x = g(x);
        y = g(g(y));
        d = gcd(absT((x - y).v), n);
    }
    return (d == 1) ? n : d;
}

/**
 * Pollard's Rho factorization algorithm.
 *
 * Pollard's Rho algorithm applied iteratively with `k` and `a` being increased in each iteration.
 * By trying different `k` and `a` parameters, the algorithm significantly reduces a chance of
 * factorization failure.
 */
template<typename I>
I pollard_rho_repeated(const I& n, const I& max_iter = 20, const I& max_inner_iter = 1000000) {
    for (I k = 2; k <= max_iter; k++) {
        I d = pollard_rho(n, k, k, max_inner_iter);
        if (d != n) return d;
    }
    return n;
}

/**
 * Factors integer `n` using a general-purpose factoring algorithm.
 */
template<typename I>
std::vector<std::pair<I, int>> factor_integer(const I& n, int max_iter = 20) {
    std::vector<std::pair<I, int>> vf;
    if (n == 0 || n == 1) return vf;
    std::vector<I> q = { n };
    while (!q.empty()) {
        I a = q.back(); q.pop_back();
        if (a == 1) {
            continue;
        }
        if (miller_rabin<I>(a)) {
            // a prime factor found
            int e = 1;
            for (auto& b : q) {
                while (b % a == 0) b /= a, e++;
            }
            vf.push_back({ a, e });
            continue;
        }
        // `a` is composite
        I d = pollard_rho_repeated<I>(a, max_iter);
        if (d == 1 || d == a) {
            // failed to factor the composite
            vf.push_back({ a, 1 });
            continue;
        }
        // a non-trivial factorization `a = d * e`
        q.push_back(d);
        q.push_back(a / d);
    }
    return vf;
}

/**
 * Factors integer `n` by trial division.
 */
template<typename I>
std::vector<std::pair<I, int>> factor_integer_slow(I n) {
    std::vector<std::pair<I, int>> vf;
    for (I i = 2; i <= n / i; i += 1) {
        if (n % i != 0) continue;
        int e = 0; while (n % i == 0) n /= i, e++;
        vf.push_back({ i, e });
    }
    if (n > 1) vf.push_back({ n, 1 });
    return vf;
}

/**
 * Factors out all factors `p` out of `n`.
 */
template<typename I, typename P, typename E>
I factor_out(I n, P p, E& e) {
    while (n % p == 0) n /= p, e++;
    return n;
}

/**
 * Reconstructs number from its factorization.
 */
template<typename P, typename I = P>
I from_factorization(const std::vector<std::pair<P, int>>& vf) {
    I r = 1;
    for (const auto& f : vf) {
        r *= powT<I>(f.first, f.second);
    }
    return r;
}

/**
 * Jointly reduces the given lists of numerators and denominators.
 *
 *    n     n_0 * ... * n_l1
 *   --- = ------------------
 *    d     d_0 * ... * d_l2
 */
template<typename I, typename P, typename G>
void fraction_reduce(std::vector<I>& numerators, std::vector<P>& denominators, G gcd) {
    for (auto& d : denominators) {
        for (size_t i = 0; d > 1 && i < numerators.size();) {
            auto g = gcd(numerators[i], d);
            if (g > 1) {
                d /= g, numerators[i] /= g;
            } else {
                i++;
            }
        }
    }
}

} // math
} // altruct
