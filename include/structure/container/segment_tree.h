#pragma once

#include <vector>
#include <iterator>
#include <functional>

namespace altruct {
namespace container {

/**
 * Segment tree that supports range operations.
 *
 * Space complexity: `O(n)`.
 * Time complexities:
 *   build: `O(n)`
 *   set: `O(log n)`
 *   get: `O(log n)`, or more precisely `O(log dist)`
 *
 * param T  - element type
 * param f  - associative functor; i.e. `f(f(a, b), c) = f(a, f(b, c))`
 *            commutativity is not required.
 * param id - neutral element with respect to `f`; i.e. `f(e, id) = f(id, e) = e`.
 *            e.g. `0` for addition, `1` for multiplication, `+inf` for minimum, etc.
 */
template<typename T>
class segment_tree {
public:
	std::vector<T> v;
	std::function<T(T, T)> f;

	segment_tree() {}
	
	segment_tree(size_t sz, std::function<T(T, T)> f, T id = T()) : f(f) {
		v.resize(make_pow2(sz) * 2, id);
		rebuild();
	}

	template<typename It>
	segment_tree(It begin, It end, std::function<T(T, T)> f, T id = T()) : f(f) {
		auto sz = std::distance(begin, end);
		v.resize(make_pow2(sz) * 2, id);
		std::copy(begin, end, v.begin() + size());
		rebuild();
	}

	void set(size_t index, const T& t) {
		index += size();
		v[index] = t;
		while ((index /= 2) > 0) {
			update(index);
		}
	}

	// If the returned element is being modified,
	// the index won't be updated automatically.
	// Use `set` instead to update immediately,
	// or `rebuild` after all modifications.
	T& operator[] (size_t index) {
		index += size();
		return v[index];
	}

	const T& operator[] (size_t index) const {
		index += size();
		return v[index];
	}

	T get(size_t index) const {
		index += size();
		return v[index];
	}
	
	T get(size_t begin, size_t end) const {
		T tl = v[0], tr = v[0]; // id
		size_t i = size();
		while (begin < end) {
			if (begin & 1) tl = f(tl, v[i + begin++]);
			if (end & 1) tr = f(v[i + --end], tr);
			begin /= 2, end /= 2, i /= 2;
		}
		return f(tl, tr);
	}

	size_t size() const {
		return v.size() / 2;
	}

	void rebuild() {
		for (size_t i = size() - 1; i > 0; i--) {
			update(i);
		}
	}

private:
	void update(size_t i) {
		v[i] = f(v[2 * i + 0], v[2 * i + 1]);
	}

	static size_t make_pow2(size_t sz) {
		size_t w = 1; while (w < sz) w *= 2;
		return w;
	}
};

} // container
} // altruct
