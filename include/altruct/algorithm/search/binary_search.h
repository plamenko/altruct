#pragma once

namespace altruct {
namespace search {

/**
 * Binary search
 *
 * @param lo  - lower bound (inclusive)
 * @param hi  - upper bound (exclusive)
 * @param eps - searches until `(hi - lo) < eps`.
 *              Must be 1 for integral, pointer and iterator types.
 * @param f   - `f(mid)` returns `true` iff the result is bigger than `mid`.
 *              The following invariant must hold: `f(mid) = (mid < mid0)`
 *              for some `mid0` to be found.
 * @return    - the smallest `mid0` for which `f(mid0)` is `false`
 */
template<typename T, typename D, typename F>
T binary_search2(T lo, T hi, D eps, F f) {
    T mid = lo, mid0 = hi;
    while (mid != mid0 && (hi - lo) >= eps) {
        mid0 = mid;
        mid = lo + (hi - lo) / 2;
        if (f(mid)) {
            lo = mid + eps;
        } else {
            hi = mid;
        }
    }
    return lo;
}

/**
 * Finds the smallest argument `x` in the range `[lo, hi)` for which `f(x) >= val`.
 *
 * See `binary_search2` for more details.
 */
template<typename T, typename D, typename V, typename F>
T lower_bound(T lo, T hi, D eps, F f, const V& val, bool decreasing = false) {
    if (decreasing) {
        return binary_search2(lo, hi, eps, [&](const T& mid){ return f(mid) > val; });
    } else {
        return binary_search2(lo, hi, eps, [&](const T& mid){ return f(mid) < val; });
    }
}

/**
 * Finds the smallest argument `x` in the range `[lo, hi)` for which `f(x) > val`.
 *
 * See `binary_search2` for more details.
 */
template<typename T, typename D, typename V, typename F>
T upper_bound(T lo, T hi, D eps, F f, const V& val, bool decreasing = false) {
    if (decreasing) {
        return binary_search2(lo, hi, eps, [&](const T& mid){ return f(mid) >= val; });
    } else {
        return binary_search2(lo, hi, eps, [&](const T& mid){ return f(mid) <= val; });
    }
}

} // search
} // altruct
