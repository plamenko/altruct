#pragma once

#include <vector>

namespace altruct {
namespace container {

/**
 * A map whose keys can only be `k` or `floor(max_key / k)` for `1 <= k <= max_lo_key`.
 *
 * The first `max_lo_key` values are kept in a vector at the index `k`.
 * Values with a key larger than `max_lo_key` are kept in another vector at the index `floor(max_key / k)`.
 *
 * `max_lo_key` must be at least `floor(sqrt(max_key))`.
 *
 * Time complexity for insert/erase/at is O(1).
 * Space complexity is O(max_lo_key).
 */
template<typename I, typename T>
class sqrt_map {
	I max_lo_key, max_key;

	std::vector<int> cnt_lo;
	std::vector<T> tbl_lo;
	std::vector<int> cnt_hi;
	std::vector<T> tbl_hi;


	int get_and_set(int& cnt, int new_cnt) { int prev_cnt = cnt; cnt = new_cnt; return prev_cnt; }

public:
	typedef typename I key_type;
	typedef typename T mapped_type;
	typedef typename std::pair<I, T> value_type;

	sqrt_map(I max_lo_key, I max_key) :
		sqrt_map(std::max(max_lo_key, max_key / max_lo_key)) {
		reset_max(max_key);
	}
	
	sqrt_map(I max_lo_key) :
		max_lo_key(max_lo_key),
		cnt_lo(max_lo_key + 1, 0),
		tbl_lo(max_lo_key + 1) {}

	void reset_max(I max_key = 0) {
		this->max_key = max_key;
		cnt_hi.assign(max_key / (max_lo_key + 1) + 1, 0);
		tbl_hi.resize(max_key / (max_lo_key + 1) + 1);
	}

	int count(const I& k) const {
		return (k <= max_lo_key) ? cnt_lo[k] : cnt_hi[max_key / k];
	}

	const T& at(const I& k) const {
		return (k <= max_lo_key) ? tbl_lo.at(k) : tbl_hi.at(max_key / k);
	}

	T& at(const I& k) {
		return (k <= max_lo_key) ? tbl_lo.at(k) : tbl_hi.at(max_key / k);
	}

	T& operator [] (const I& k) {
		if (k <= max_lo_key) {
			cnt_lo[k] = 1;
			return tbl_lo[k];
		} else {
			cnt_hi[max_key / k] = 1;
			return tbl_hi[max_key / k];
		}
	}

	std::pair<I, bool> insert(const value_type& e) {
		const I& k = e.first; const T& v = e.second;
		if (k <= max_lo_key) {
			if (cnt_lo[k]) return{ k, false };
			tbl_lo[k] = v;
			cnt_lo[k] = 1;
			return{ k, true };
		} else {
			if (cnt_hi[max_key / k]) return{ k, false };
			tbl_hi[max_key / k] = v;
			cnt_hi[max_key / k] = 1;
			return{ k, true };
		}
	}

	int erase(const I& k) {
		return (k <= max_lo_key) ? get_and_set(cnt_lo[k], 0) : get_and_set(cnt_hi[max_key / k], 0);
	}
};

} // container
} // altruct