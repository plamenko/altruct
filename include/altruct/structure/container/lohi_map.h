#pragma once

#include <vector>
#include <unordered_map>

namespace altruct {
namespace container {

/**
 * A map that uses vector for the elements with a small enough key.
 *
 * The first `max_lo_key` values are kept in a vector at the index `k`.
 * Values with a key larger than `max_lo_key` are kept in an unordered_map.
 *
 * Time complexity for insert/erase/at is O(1) amortized.
 * Space complexity is O(max_lo_key + tbl_hi_size).
 */
template<typename I, typename T>
class lohi_map {
    I max_lo_key;

    std::vector<int> cnt_lo;
    std::vector<T> tbl_lo;
    std::unordered_map<I, T> tbl_hi;

    int get_and_set(int& cnt, int new_cnt) { int prev_cnt = cnt; cnt = new_cnt; return prev_cnt; }

public:
    typedef I key_type;
    typedef T mapped_type;
    typedef std::pair<I, T> value_type;

    lohi_map(I max_lo_key) :
        max_lo_key(max_lo_key),
        cnt_lo(max_lo_key + 1, 0),
        tbl_lo(max_lo_key + 1) {}

    int count(const I& k) const {
        return (k <= max_lo_key) ? cnt_lo[k] : (int)tbl_hi.count(k);
    }

    const T& at(const I& k) const {
        return (k <= max_lo_key) ? tbl_lo.at(k) : tbl_hi.at(k);
    }

    T& at(const I& k) {
        return (k <= max_lo_key) ? tbl_lo.at(k) : tbl_hi.at(k);
    }

    T& operator [] (const I& k) {
        return (k <= max_lo_key) ? (cnt_lo[k] = 1, tbl_lo[k]) : (tbl_hi[k]);
    }

    std::pair<I, bool> insert(const value_type& e) {
        const I& k = e.first; const T& v = e.second;
        if (k <= max_lo_key) {
            if (cnt_lo[k]) return{ k, false };
            tbl_lo[k] = v;
            cnt_lo[k] = 1;
            return{ k, true };
        } else {
            bool b = tbl_hi.insert(e).second;
            return{ k, b };
        }
    }

    int erase(const I& k) {
        return (k <= max_lo_key) ? get_and_set(cnt_lo[k], 0) : (int)tbl_hi.erase(k);
    }
};

} // container
} // altruct
