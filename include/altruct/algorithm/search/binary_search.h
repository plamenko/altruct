#pragma once

namespace altruct {
namespace search {


/**
 * Returns the first index for which `predicate` returns `true`.
 *
 * Predicate must be sorted so all `false` come before `true`.
 * `It` must be a random-access iterator, pointer, or integral type.
 *
 * @param lo  - lower bound (inclusive)
 * @param hi  - upper bound (exclusive)
 * @return the smallest index for which `predicate` returns `true`.
 */
template<typename It, typename F>
It binary_search_pred(It lo, It hi, F predicate) {
    while (lo < hi) {
        It mid = lo + (hi - lo) / 2;
        if (predicate(mid)) {
            hi = mid;
        } else {
            lo = ++mid;
        }
    }
    return lo;
}

/**
 * Returns the first value for which `predicate` returns `true`.
 *
 * Predicate must be sorted so all `false` come before `true`.
 * `V` must be a numerical type.
 *
 * @param lo  - lower bound (inclusive)
 * @param hi  - upper bound (exclusive)
 * @return the smallest value for which `predicate` returns `true`.
 */
template<typename X, typename F>
X binary_search_num(X lo, X hi, X eps, F predicate) {
    auto prev_mid = hi;
    while ((hi - lo) >= eps) {
        auto mid = lo + (hi - lo) / 2;
        if (mid == prev_mid) break;
        if (predicate(mid)) {
            hi = mid;
        } else {
            lo = mid + eps;
        }
        prev_mid = mid;
    }
    return lo;
}

/**
 * Finds the smallest argument `x` in the range `[lo, hi)` for which `f(x) >= val` when `f` is not decreasing,
 * or `f(x) <= val` when `f` is decreasing.
 *
 * See `binary_search_num` for more details.
 */
template<typename X, typename Y, typename F>
X lower_bound_num(X lo, X hi, X eps, const Y& val, F f, bool decreasing = false) {
    if (decreasing) {
        return binary_search_num(lo, hi, eps, [&](X mid){ return f(mid) <= val; });
    } else {
        return binary_search_num(lo, hi, eps, [&](X mid){ return f(mid) >= val; });
    }
}

/**
 * Finds the smallest argument `x` in the range `[lo, hi)` for which `f(x) > val` when `f` is not decreasing,
 * or `f(x) < val` when `f` is decreasing.
 *
 * See `binary_search_num` for more details.
 */
template<typename X, typename Y, typename F>
X upper_bound_num(X lo, X hi, X eps, const Y& val, F f, bool decreasing = false) {
    if (decreasing) {
        return binary_search_num(lo, hi, eps, [&](X mid){ return f(mid) < val; });
    } else {
        return binary_search_num(lo, hi, eps, [&](X mid){ return f(mid) > val; });
    }
}

} // search
} // altruct
