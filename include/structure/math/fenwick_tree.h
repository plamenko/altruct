#pragma once

#include <vector>
#include <functional>

namespace altruct {
namespace math {

/**
 * Fenwick tree.
 */
template<typename T>
struct fenwick_tree {
	std::vector<T> v;
	std::function<T(T,T)> f;

    fenwick_tree(size_t sz, std::function<T(T, T)> f, T id = 0) : v(sz + 1, id), f(f) {
    }

    void reset(T id = 0) {
        v.assign(v.size(), id);
    }

	void add(size_t index, const T& val) {
		for (index++; index < v.size(); index += index & -index) v[index] = f(v[index], val);
	}

	T get_sum(size_t index, T id = 0) {
		T r = id;
		for (index++; index > 0; index -= index & -index) r = f(r, v[index]);
		return r;
	}
};

} // math
} // altruct
