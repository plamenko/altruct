#pragma once

#include <iterator>
#include <type_traits>
#include <vector>

#include "altruct/algorithm/math/base.h"

namespace altruct {
namespace math {

/**
 * Builds a range look-up table up to `n`.
 *
 * `v[i] = i * step`
 */
template<typename It, typename T = typename std::iterator_traits<It>::value_type>
void range(It begin, It end, T step = T(1)) {
    T v = zeroT<T>::of(step);
    for (It it = begin; it != end; ++it) {
        *it = v; v += step;
    }
}
template<typename T>
std::vector<T> range(size_t n, T step = T(1)) {
    std::vector<T> v(n, step);
    range(v.begin(), v.end(), step);
    return v;
}

/**
 * Builds the powers of `base` look-up table up to `n`.
 *
 * `v[i] = base ^ i`
 */
template<typename It, typename T = typename std::iterator_traits<It>::value_type>
void powers(It begin, It end, T base) {
    T v = identityT<T>::of(base);
    for (It it = begin; it != end; ++it) {
        *it = v; v *= base;
    }
}
template<typename T>
std::vector<T> powers(size_t n, T base) {
    std::vector<T> v(n, base);
    powers(v.begin(), v.end(), base);
    return v;
}

/**
 * Builds the factorial look-up table up to `n`.
 *
 * `v[i] = i!`
 */
template<typename It, typename T = typename std::iterator_traits<It>::value_type>
void factorials(It begin, It end, T id = T(1)) {
    T v = id, i = id;
    for (It it = begin; it != end; ++it) {
        *it = v; v *= i; i += id;
    }
}
template<typename T>
std::vector<T> make_factorials(size_t n, T id = T(1)) {
    std::vector<T> v(n, id);
    factorials(v.begin(), v.end(), id);
    return v;
}

/**
 * Builds the inverse factorial look-up table up to `n`.
 *
 * `v[i] = 1 / i!`
 * `k!` and `k` can be provided to avoid computing `n!` from scratch
 */
template<typename It, typename T = typename std::iterator_traits<It>::value_type>
void inv_factorials(It begin, It end, T fact_k = T(1), size_t k = 1) {
    if (k == 0) k = 1;
    T id = identityOf(fact_k), val_k = castOf(fact_k, k);
    for (size_t n = std::distance(begin, end); ++k < n;) {
        val_k += id, fact_k *= val_k;
    }
    T ifact = id / fact_k;
    for (It it = end; it != begin; ) {
        *--it = ifact; ifact *= val_k; val_k -= id;
    }
}
template<typename T>
std::vector<T> make_inv_factorials(size_t n, T fact_k = T(1), size_t k = 1) {
    std::vector<T> v(n, identityOf(fact_k));
    inv_factorials(v.begin(), v.end(), fact_k, k);
    return v;
}

/**
 * Builds the inverses look-up table up to `n`.
 *
 * `v[i] = 1 / i`
 */
template<typename It, typename T = typename std::iterator_traits<It>::value_type>
void inverses_from_ifact(It begin, It end, T id = T(1)) {
    T fact = id, i = id;
    *begin = zeroOf(id);
    for (It it = ++begin; it != end; ++it) {
        *it *= fact; fact *= i; i += id;
    }
}
template<typename It, typename T = typename std::iterator_traits<It>::value_type>
void inverses(It begin, It end, T id = T(1)) {
    inv_factorials(begin, end, id);
    inverses_from_ifact(begin, end, id);
}
template<typename T>
std::vector<T> make_inverses(size_t n, T id = T(1)) {
    std::vector<T> v(n, id);
    inverses(v.begin(), v.end(), id);
    return v;
}

/**
 * Raises the elements to the `n`-th power.
 *
 * `v[i] <- v[i] ^ n`
 */
template<typename It, typename I>
void power(It begin, It end, I n) {
    for (It it = begin; it != end; ++it) {
        *it = powT(*it, n);
    }
}

/**
 * Inverts the table elements with respect to multiplication.
 * Note: Zero values are left zeros.
 *
 * `v[i] <- 1 / v[i]`
 */
template<typename It, typename T = typename std::iterator_traits<It>::value_type>
void invert(It begin, It end, T id = T(1)) {
    T e0 = zeroT<T>::of(id);
    for (It it = begin; it != end; ++it) {
        if (*it != e0) *it = id / *it;
    }
}

/**
 * Negates the table elements.
 *
 * `v[i] <- -v[i]`
 */
template<typename It>
void negate(It begin, It end) {
    for (It it = begin; it != end; ++it) {
        *it = -*it;
    }
}

/**
 * Alternates the sign of the table elements.
 * I.e. negates odd elements.
 *
 * `v[i] <- v[i] * (-1)^i`
 */
template<typename It>
void alternate(It begin, It end) {
    int s = 1;
    for (It it = begin; it != end; ++it, s = -s) {
        if (s < 0) *it = -*it;
    }
}

/**
 * Accumulates the table elements.
 *
 * `v[i] <- Sum[v[j], {j, 0, i}]`
 */
template<typename It, typename T = typename std::iterator_traits<It>::value_type>
void accumulate(It begin, It end) {
    if (begin == end) return;
    It prev = begin++;
    for (It it = begin; it != end; ++it) {
        *it += *prev++;
    }
}

/**
 * Differences between the table elements.
 *
 * `v[i] <- v[i] - v[i - 1]`
 */
template<typename It, typename T = typename std::iterator_traits<It>::value_type>
void differentiate(It begin, It end) {
    if (begin == end) return;
    It prev = --end;
    for (It it = end; it != begin; --it) {
        *it -= *--prev;
    }
}

} // math
} // altruct
