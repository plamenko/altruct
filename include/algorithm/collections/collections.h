#pragma once

#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>

namespace altruct {
namespace collections {

/**
 * Returns a sorted vector of the elements in the range [begin, end).
 */
template<
    typename It,
    typename T = typename std::iterator_traits<It>::value_type
>
std::vector<T> sorted(It begin, It end) {
    std::vector<T> r(begin, end);
    std::sort(r.begin(), r.end());
    return r;
}

/**
 * Returns a sorted vector of the elements in the collection `c`.
 */
template<
    typename C,
    typename T = typename C::value_type
>
std::vector<T> sorted(const C& c) {
    return sorted(c.begin(), c.end());
}

/**
 * Returns a reversed vector of the elements in the range [begin, end).
 */
template<
    typename It,
    typename T = typename std::iterator_traits<It>::value_type
>
std::vector<T> reversed(It begin, It end) {
    std::vector<T> r(begin, end);
    std::reverse(r.begin(), r.end());
    return r;
}

/**
 * Returns a reversed vector of the elements in the collection `c`.
 */
template<
    typename C,
    typename T = typename C::value_type
>
std::vector<T> reversed(const C& c) {
    return reversed(c.begin(), c.end());
}

/**
 * Returns a vector of the the first `n` elements in the range [begin, end).
 */
template<
    typename It,
    typename T = typename std::iterator_traits<It>::value_type
>
std::vector<T> take(It begin, It end, size_t n) {
    std::vector<T> r;
    for (It it = begin; it != end && n-- > 0; ++it) {
        r.push_back(*it);
    }
    return r;
}

/**
 * Returns a vector of the the first `n` elements in the collection `c`.
 */
template<
    typename C,
    typename T = typename C::value_type
>
std::vector<T> take(const C& c, size_t n) {
    return take(c.begin(), c.end(), n);
}

/**
 * Returns a vector of the elements in the range [begin, end) that satisfy predicate `p`.
 */
template<
	typename It,
	typename P,
	typename T = typename std::iterator_traits<It>::value_type
>
std::vector<T> filter(It begin, It end, P p) {
	std::vector<T> r;
	for (It it = begin; it != end; ++it) {
		if (p(*it)) r.push_back(*it);
	}
	return r;
}

/**
 * Returns a vector of the elements in the collection `c` that satisfy predicate `p`.
 */
template<
	typename C,
	typename P,
	typename T = typename C::value_type
>
std::vector<T> filter(const C& c, P p) {
	return filter(c.begin(), c.end(), p);
}

/**
 * Returns a vector of the elements in the range [begin, end) transformed by functor `f`.
 */
template<
	typename It,
	typename F,
	typename T = typename std::result_of<F(typename std::iterator_traits<It>::value_type)>::type
>
std::vector<T> transform(It begin, It end, F f) {
	std::vector<T> r;
	for (It it = begin; it != end; ++it) {
		r.push_back(f(*it));
	}
	return r;
}

/**
 * Returns a vector of the elements in the collection `c` transformed by functor `f`.
 */
template<
	typename C,
	typename F,
	typename T = typename std::result_of<F(typename C::value_type)>::type
>
std::vector<T> transform(const C& c, F f) {
	return transform(c.begin(), c.end(), f);
}

/**
 * Returns Run-Length encoding of the elements in the range [begin, end).
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
 * Returns Run-Length encoding of the elements in collection `c`.
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
template<typename It1, typename It2>
int compare(It1 b1, It1 e1, It2 b2, It2 e2, size_t max_len = -1) {
	while (b1 != e1 && b2 != e2 && max_len != 0) {
		if (!(*b1 == *b2)) return (*b1 < *b2) ? -1 : +1;
		++b1, ++b2, --max_len;
	}
	if (max_len == 0) return 0;
	if (b2 != e2) return -1;
	if (b1 != e1) return +1;
	return 0;
}

/**
 * Reserves space for at least `c.size() + sz` elements,
 * while maintaining exponential growth.
 */
template<typename C>
void reserve_more(C& c, size_t sz) {
	if (c.size() + sz <= c.capacity()) return;
	c.reserve(std::max(c.size() + sz, c.capacity() + c.capacity() / 2));
}

} // collections
} // altruct
