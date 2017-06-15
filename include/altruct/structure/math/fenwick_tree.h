#pragma once

#include <vector>
#include <functional>

namespace altruct {
namespace math {

/**
 * Fenwick tree.
 */
template<typename T, typename F = std::function<T(T, T)>>
struct fenwick_tree {
    std::vector<T> v;
    F f;

    fenwick_tree(size_t sz, const F& f, T id = 0) : v(sz + 1, id), f(f) {
    }

    void reset(T id = 0) {
        v.assign(v.size(), id);
    }

    void add(size_t index, const T& val) {
        for (index++; index < v.size(); index += lo_bit(index)) v[index] = f(v[index], val);
    }

    T get_sum(size_t index, T id = 0) {
        T r = id;
        for (index++; index > 0; index -= lo_bit(index)) r = f(r, v[index]);
        return r;
    }

    size_t lo_bit(size_t index) {
        return index & -index;
    }
};

} // math
} // altruct
