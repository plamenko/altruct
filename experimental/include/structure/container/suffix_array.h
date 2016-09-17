#pragma once

#include "structure/container/range_minimum_query.h"

#include <vector>
#include <iterator>
#include <functional>

namespace altruct {
namespace container {

/**
 * Suffix array structure.
 * Also computes inverse suffix array, lcp and its rmq.
 *
 * Space complexity: `O(n)`.
 * Time complexities:
 *   build: `O(n)`
 *   query: `O(1)`
 */
template<typename ALPHA_T = char, typename INDEX_T = int>
class suffix_array {
public:
	typedef ALPHA_T alpha_t;
	typedef INDEX_T index_t;

	std::vector<alpha_t> _string;     // input string (we only need this for `compare_substrings`)
	std::vector<index_t> _suff_arr;   // position (offset in the string) of the k-th lexicographically smallest suffix
	std::vector<index_t> _suff_ord;   // lexicographical order of the suffix starting at position `i`; inverse of `_suff_arr`
	//std::vector<index_t> _lcp_arr;  // longest-common-prefix for each pair of successive sorted suffixes
	direct_rmq<index_t> _lcp_arr_rmq; // range-minimum-query structure for lcp array

	template<typename It>
	suffix_array(It begin, It end) {
		build_all(begin, end);
	}

	template<typename It>
	void build_all(It begin, It end) {
		build_suffix_array(begin, end);
		build_inverse_suffix_array();
		build_lcp_array();
	}

	template<typename It>
	void build_suffix_array(It begin, It end) {
		alpha_t max_elem = *std::max_element(begin, end);
		build_suffix_array(begin, end, (index_t)max_elem + 1);
	}

	template<typename It>
	void build_suffix_array(It begin, It end, index_t alpha_size) {
		_string.assign(begin, end);
		_suff_arr.assign(_string.size() + 1, 0);
		if (_string.size() > 0) {
			sa_is<alpha_t>(_string.data(), (index_t)_string.size(), alpha_size, _suff_arr.data());
		}
	}

	void build_inverse_suffix_array() {
		index_t n = size();
		_suff_ord.resize(n + 1);
		for (index_t i = 0; i <= n; i++) {
			_suff_ord[_suff_arr[i]] = i;
		}
	}

	// Construction by SA-IS algorithm

	template<typename alphat>
	static void sa_is(const alphat *str, index_t n, index_t alpha_size, index_t *sa) {
		std::vector<index_t> bucket_offsets(std::max(alpha_size, (n + 1) / 2) + 1);
		sa_is<alphat>(str, n, alpha_size, sa, bucket_offsets);
	}

	template<typename alphat>
	static void sa_is(const alphat *str, index_t n, index_t alpha_size, index_t *sa, std::vector<index_t>& bucket_offsets) {
		std::vector<bool> types(n + 1);
		types[n - 1] = 0; types[n] = 1;
		for (index_t i = n - 2; i >= 0; i--) {
			types[i] = (str[i] < str[i + 1]) || (str[i] == str[i + 1]) && types[i + 1];
		}

		count_alphabets(str, n, alpha_size, bucket_offsets);
		get_bucket_offsets(str, n, true, alpha_size, bucket_offsets);
		std::fill(sa, sa + n + 1, -1);
		for (index_t i = 1; i < n; i++) {
			if (types[i] && !types[i - 1]) sa[--bucket_offsets[(index_t)str[i]]] = i;
		}
		sa[0] = n;
		induced_sort(str, n, alpha_size, types, sa, bucket_offsets);

		index_t n1 = 0;
		for (index_t i = 0; i <= n; i++) {
			index_t j = sa[i];
			if (j > 0 && types[j] && !types[j - 1]) sa[n1++] = j;
		}

		index_t *buffer = sa + n1;
		std::fill(buffer, sa + n + 1, -1);
		index_t unique_lms_count = 0, prev_pos = -1;
		buffer[sa[0] / 2] = unique_lms_count++;
		for (index_t i = 1; i < n1; i++) {
			index_t pos = sa[i];
			bool diff = false;
			if (prev_pos == -1) {
				diff = true;
			} else {
				for (index_t j = pos, k = prev_pos;; j++, k++) {
					if (str[j] != str[k] || types[j] != types[k]) {
						diff = true;
						break;
					} else if (j != pos && ((types[j] && !types[j - 1]) || (types[k] && !types[k - 1]))) {
						break;
					}
				}
			}
			if (diff) {
				unique_lms_count++;
				prev_pos = pos;
			}
			buffer[pos / 2] = unique_lms_count - 1;
		}
		for (index_t i = n, j = n; i >= n1; i--) {
			if (sa[i] >= 0) sa[j--] = sa[i];
		}

		index_t *sa1 = sa, *s1 = sa + n + 1 - n1;
		if (unique_lms_count == n1) {
			for (index_t i = 0; i < n1; i++) sa1[s1[i]] = i;
		} else {
			sa_is<index_t>(s1, n1 - 1, unique_lms_count, sa1, bucket_offsets);
		}

		count_alphabets(str, n, alpha_size, bucket_offsets);
		get_bucket_offsets(str, n, true, alpha_size, bucket_offsets);
		for (index_t i = 1, j = 0; i <= n; i++) {
			if (types[i] && !types[i - 1]) s1[j++] = i;
		}
		for (index_t i = 0; i < n1; i++) sa1[i] = s1[sa1[i]];
		std::fill(sa + n1, sa + n + 1, -1);
		for (index_t i = n1 - 1; i >= 1; i--) {
			index_t j = sa[i]; sa[i] = -1;
			sa[--bucket_offsets[(index_t)str[j]]] = j;
		}
		induced_sort(str, n, alpha_size, types, sa, bucket_offsets);
	}

	template<typename alphat>
	static void induced_sort(const alphat *str, index_t n, index_t alpha_size, const std::vector<bool>& types, index_t *sa, std::vector<index_t>& bucket_offsets) {
		get_bucket_offsets(str, n, false, alpha_size, bucket_offsets);
		for (index_t i = 0; i < n; i++) {
			index_t j = sa[i] - 1;
			if (j >= 0 && !types[j]) sa[bucket_offsets[(index_t)str[j]]++] = j;
		}
		get_bucket_offsets(str, n, true, alpha_size, bucket_offsets);
		for (index_t i = n; i >= 1; i--) {
			index_t j = sa[i] - 1;
			if (j >= 0 && types[j]) sa[--bucket_offsets[(index_t)str[j]]] = j;
		}
	}

	template<typename alphat>
	static void get_bucket_offsets(const alphat *str, index_t n, bool dir, index_t alpha_size, std::vector<index_t>& bucket_offsets) {
		typename std::vector<index_t>::iterator alphabet_counts;
		if ((index_t)bucket_offsets.size() / 2 < alpha_size) {
			count_alphabets(str, n, alpha_size, bucket_offsets, true);
			alphabet_counts = bucket_offsets.begin();
		} else {
			alphabet_counts = bucket_offsets.begin() + alpha_size;
		}
		index_t cumsum = 1;
		if (dir) {
			for (int i = 0; i < alpha_size; i++) {
				cumsum += alphabet_counts[i];
				bucket_offsets[i] = cumsum;
			}
		} else {
			for (int i = 0; i < alpha_size; i++) {
				index_t x = alphabet_counts[i];
				bucket_offsets[i] = cumsum;
				cumsum += x;
			}
		}
	}

	template<typename alphat>
	static void count_alphabets(const alphat *str, index_t n, index_t alpha_size, std::vector<index_t>& bucket_offsets, bool b = false) {
		if (b || (index_t)bucket_offsets.size() / 2 >= alpha_size) {
			typename std::vector<index_t>::iterator alphabet_counts = bucket_offsets.begin() + (b ? 0 : alpha_size);
			std::fill(alphabet_counts, alphabet_counts + alpha_size, 0);
			for (index_t i = 0; i < n; i++) {
				alphabet_counts[(index_t)str[i]]++;
			}
		}
	}


	// computes the length of the largest-common-prefix
	// for each pair of successive (sorted) suffixes and
	// preprocesses range-minimum-query for lcp array
	void build_lcp_array() {
		index_t n = size();
		std::vector<index_t> _lcp_arr(n + 2);
		index_t h = 0;
		for (index_t i = 0; i < n; i++) {
			index_t ord = _suff_ord[i];
			index_t j = _suff_arr[ord - 1];
			index_t hmax = std::min(n - j, n - i);
			for (index_t k = 0; h < hmax && _string[i + h] == _string[j + h]; ++h);
			_lcp_arr[ord - 1] = h;
			if (h > 0) --h;
		}
		_lcp_arr_rmq.build(_lcp_arr.begin(), _lcp_arr.end());
	}

	// gets the longest-common-prefix of two suffixes
	// (starting at positions `i` and `j` respectively)
	index_t get_lcp(index_t i, index_t j) const {
		index_t n = size();
		if (i == j) return n - i;
		index_t x = _suff_ord[i], y = _suff_ord[j];
		if (x > y) std::swap(x, y);
		return _lcp_arr_rmq.get_value(x, y);
	}

	// compares two substrings of the input string
	int compare_substrings(index_t begin1, index_t end1, index_t begin2, index_t end2) const {
		if (begin1 == begin2) return _compare(end1, end2);
		index_t len = get_lcp(begin1, begin2);
		if (begin1 + len < end1 && begin2 + len < end2) {
			return _compare(_string[begin1 + len], _string[begin2 + len]);
		} else {
			return _compare(end1 - begin1, end2 - begin2);
		}
	}
	template<typename T>
	static int _compare(T a, T b) { return (a < b) ? -1 : ((b < a) ? +1 : 0); }

	// returns the position (offset in the input string)
	// of the k-th lexicographically smallest suffix
	index_t get_kth_suffix(index_t k) const { return _suff_arr[k]; }

	// returns the size of the input string.
	// note that there is one more suffix, the empty one.
	index_t size() const { return (index_t)_string.size(); }
};

} // container
} // altruct
