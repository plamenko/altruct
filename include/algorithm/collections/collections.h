#pragma once

#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>

namespace altruct {
namespace collections {

/**
 * Returns the vector of the elements in the range [begin, end) that satisfy predicate `p`.
 */
template<
	typename It,
	typename P,
	typename T = typename std::iterator_traits<It>::value_type
>
std::vector<T> filter(It begin, It end, const P& p) {
	std::vector<T> r;
	for (It it = begin; it != end; ++it) {
		if (p(*it)) r.push_back(*it);
	}
	return r;
}

/**
 * Returns the vector of the elements in collection `c` that satisfy predicate `p`.
 */
template<
	typename C,
	typename P,
	typename T = typename C::value_type
>
std::vector<T> filter(const C& c, const P& p) {
	return filter(c.begin(), c.end(), p);
}

/**
 * Returns the vector of the elements in the range [begin, end) transformed by functor `f`.
 */
template<
	typename It,
	typename F,
	typename T = typename std::result_of<F(typename std::iterator_traits<It>::value_type)>::type
>
std::vector<T> transform(It begin, It end, const F& f) {
	std::vector<T> r;
	for (It it = begin; it != end; ++it) {
		r.push_back(f(*it));
	}
	return r;
}

/**
 * Returns the vector of the elements in collection `c` transformed by functor `f`.
 */
template<
	typename C,
	typename F,
	typename T = typename std::result_of<F(typename C::value_type)>::type
>
std::vector<T> transform(const C& c, const F& f) {
	return transform(c.begin(), c.end(), f);
}

/**
 * Returns the Run-Length encoding of the elements in the range [begin, end).
 */
template<
	typename It,
	typename T = typename std::iterator_traits<It>::value_type
>
std::vector<std::pair<T, int>> run_length(It begin, It end) {
	std::vector<std::pair<T, int>> r;
	for (It it = begin; it != end; ++it) {
		if (!r.empty() && r.back().first == *it) {
			r.back().second++;
		} else {
			r.push_back({ *it, 1 });
		}
	}
	return r;
}

/**
 * Returns the Run-Length encoding of the elements in collection `c`.
 */
template<
	typename C,
	typename T = typename C::value_type
>
std::vector<std::pair<T, int>> run_length(const C& c) {
	return run_length(c.begin(), c.end());
}

} // collections
} // altruct
