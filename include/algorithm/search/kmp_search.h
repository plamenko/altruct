#pragma once

namespace altruct {
namespace search {

/**
 * Knuth–Morris–Pratt searching algorithm.
 *
 * Searches for and reports all occurences of `pattern` within the `text`.
 *
 * Complexity: `O(n)`
 *
 * @param t - text
 * @param n - length of the text
 * @param p - pattern
 * @param m - length of the pattern
 * @param callback - functor `bool callback(size_t pos);`
 *                   called for each match at position `pos`;
 *                   returns whether the search should continue;
 * @return - the last reported match position, or `n` if no matches found
 */
template<typename RanIt, typename F>
size_t kmp_search(RanIt t, size_t n, RanIt p, size_t m, F callback) {
	if (m == 0) return 0;
	if (m > n) return n;
	// preprocess
	vector<size_t> b(m + 1, -1);
	for (size_t i = 0, j = -1; i < m;) {
		while (j != -1 && p[i] != p[j]) j = b[j];
		i++; j++;
		b[i] = j;
	}
	// search
	size_t r = n;
	for (size_t i = 0, j = 0; i < n;) {
		while (j != -1 && t[i] != p[j]) j = b[j];
		i++; j++;
		if (j == m) {
			r = i - j;
			if (!callback(r)) break;
			j = b[j];
		}
	}
	return r;
}

} // search
} // altruct
