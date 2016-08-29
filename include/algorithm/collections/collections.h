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

/**
 * Compares two sequences.
 *
 * @param max_len - compares at most `max_len` elements.
 * @return - integer `-1`, `0`, `+1`, based on the comparison result.
 */
template<typename It>
int compare(It b1, It e1, It b2, It e2, size_t max_len = -1) {
	while (b1 != e1 && b2 != e2 && max_len != 0) {
		if (!(*b1 == *b2)) return (*b1 < *b2) ? -1 : +1;
		++b1, ++b2, --max_len;
	}
	if (max_len == 0) return 0;
	if (b2 != e2) return -1;
	if (b1 != e1) return +1;
	return 0;
}

} // collections
} // altruct
