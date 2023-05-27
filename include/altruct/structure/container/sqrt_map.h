#pragma once

#include <vector>
#include <stdexcept>

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

    std::vector<char> cnt_lo;
    std::vector<T> tbl_lo;
    std::vector<char> cnt_hi;
    std::vector<T> tbl_hi;


    int get_and_set(char& cnt, char new_cnt) { int prev_cnt = cnt; cnt = new_cnt; return prev_cnt; }

    void throw_oor() const {
        throw std::out_of_range("invalid sqrt_map<I, T> key");
    }

public:
    typedef I key_type;
    typedef T mapped_type;
    typedef std::pair<I, T> value_type;

    sqrt_map(I max_lo_key, I max_key) :
        sqrt_map(max_lo_key) {
        reset_max(max_key);
    }

    sqrt_map(I max_lo_key) :
        max_lo_key(max_lo_key),
        cnt_lo(max_lo_key + 1, 0),
        tbl_lo(max_lo_key + 1) {}

    void swap(sqrt_map& rhs) {
        std::swap(max_lo_key, rhs.max_lo_key);
        std::swap(max_key, rhs.max_key);
        cnt_lo.swap(rhs.cnt_lo);
        tbl_lo.swap(rhs.tbl_lo);
        cnt_hi.swap(rhs.cnt_hi);
        tbl_hi.swap(rhs.tbl_hi);
    }

    void reset_max(I max_key = 0) {
        this->max_key = max_key;
        cnt_hi.assign(max_key / max_lo_key + 1, 0);
        tbl_hi.resize(max_key / max_lo_key + 1);
    }

    int count(const I& k) const {
        return (k <= max_lo_key) ? cnt_lo[k] : cnt_hi[max_key / k];
    }

    // unchecked low element access
    const T& lo(const I& k) const {
        return tbl_lo[k];
    }

    // unchecked low element access
    T& lo(const I& k) {
        return tbl_lo[k];
    }

    // unchecked high element access
    const T& hi(const I& k) const {
        return tbl_hi[k];
    }

    // unchecked high element access
    T& hi(const I& k) {
        return tbl_hi[k];
    }

    // unchecked element access
    const T& el(const I& k) const {
        return (k <= max_lo_key) ? tbl_lo[k] : tbl_hi[max_key / k];
    }

    // unchecked element access
    T& el(const I& k) {
        return (k <= max_lo_key) ? tbl_lo[k] : tbl_hi[max_key / k];
    }

    // unchecked element access as function
    const T& operator () (const I& k) const {
        return (k <= max_lo_key) ? tbl_lo[k] : tbl_hi[max_key / k];
    }

    // unchecked element access
    const T& operator [] (const I& k) const {
        return (k <= max_lo_key) ? tbl_lo[k] : tbl_hi[max_key / k];
    }

    // unchecked element access, creates default if it doesn't exist
    T& operator [] (const I& k) {
        if (k <= max_lo_key) {
            cnt_lo[k] = 1;
            return tbl_lo[k];
        } else {
            I i = max_key / k;
            cnt_hi[i] = 1;
            return tbl_hi[i];
        }
    }

    // checked element access, throws if it doesn't exist
    const T& at(const I& k) const {
        if (k <= max_lo_key) {
            if (!cnt_lo.at(k)) throw_oor();
            return tbl_lo[k];
        } else {
            I i = max_key / k;
            if (!cnt_hi.at(i)) throw_oor();
            return tbl_hi[i];
        }
    }

    // checked element access, throws if it doesn't exist
    T& at(const I& k) {
        if (k <= max_lo_key) {
            if (!cnt_lo.at(k)) throw_oor();
            return tbl_lo[k];
        } else {
            I i = max_key / k;
            if (!cnt_hi.at(i)) throw_oor();
            return tbl_hi[i];
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
            I i = max_key / k;
            if (cnt_hi[i]) return{ k, false };
            tbl_hi[i] = v;
            cnt_hi[i] = 1;
            return{ k, true };
        }
    }

    int erase(const I& k) {
        return (k <= max_lo_key) ? get_and_set(cnt_lo[k], 0) : get_and_set(cnt_hi[max_key / k], 0);
    }
};

} // container
} // altruct
