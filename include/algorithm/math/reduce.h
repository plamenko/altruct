#pragma once

#include "algorithm/math/base.h"

#include <algorithm>
#include <limits>

namespace altruct {
namespace math {

/**
 * Reduces the elements of the container by a given functor.
 *
 * @param c - container to apply the functor over
 * @param f - functor that reduces two elements
 * @param id - identity element
 */
template<typename C, typename F, typename T>
T reduce(const C& c, const F& f, T id) {
	T r = id;
	for (const auto& e : c) {
		r = f(r, e);
	}
	return r;
}

/**
 * Sum of the elements in the container.
 *
 * @param c - container to apply sum over
 * @param id - identity element, default 0
 */
template<typename C, typename T = typename C::value_type>
T reduce_sum(const C& c, T id = T(0)) {
	return reduce(c, [](const T& r, const T& e) { return r + e;  }, id);
}

/**
 * Product of the elements in the container.
 *
 * @param c - container to apply product over
 * @param id - identity element, default 1
 */
template<typename C, typename T = typename C::value_type>
T reduce_product(const C& c, T id = T(1)) {
	return reduce(c, [](const T& r, const T& e) { return r * e;  }, id);
}

/**
 * Minimum of the elements in the container.
 *
 * @param c - container to apply minimum over
 * @param id - identity element, default +infinity
 */
template<typename C, typename T = typename C::value_type>
T reduce_min(const C& c, T id = +std::numeric_limits<T>::max()) {
	return reduce(c, [](const T& r, const T& e) { return e < r ? e : r; }, id);
}

/**
 * Maximum of the elements in the container.
 *
 * @param c - container to apply maximum over
 * @param id - identity element, default -infinity
 */
template<typename C, typename T = typename C::value_type>
T reduce_max(const C& c, T id = -std::numeric_limits<T>::max()) {
	return reduce(c, [](const T& r, const T& e) { return r < e ? e : r; }, id);
}

/**
 * Minimum excludant of the elements in the container.
 *
 * @param c - a sorted container to apply mex over
 * @param id - identity element, default 0
 */
template<typename C, typename T = typename C::value_type>
T reduce_mex(const C& c, T id = 0) {
	return reduce(c, [](const T& r, const T& e) { return r < e ? r : r + 1; }, id);
}

} // math
} // altruct
