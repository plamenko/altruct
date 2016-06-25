#pragma once

#include <vector>
#include <iterator>
#include <functional>

namespace altruct {
namespace container {

/**
 * Interval tree that supports range operations.
 *
 * Space complexity: `O(n)`.
 * Time complexities:
 *   build:  `O(n)`
 *   update: `O(log n)`
 *   get:    `O(log n)`
 *
 * param T      - node type
 * param f_up   - associative functor for upward propagation; commutativity is not required.
 * param f_down - associative functor for lazy downward propagation; commutativity is not required.
 * param id     - neutral element with respect to `f`; i.e. `f(e, id) = f(id, e) = e`.
 *                e.g. `0` for addition, `1` for multiplication, `+inf` for minimum, etc.
 */
template<typename T>
class interval_tree {
public:
	typedef std::function<void(T& node)> update_functor;
	typedef std::function<void(T& parent, const T& left, const T& right)> up_functor;
	typedef std::function<void(T& parent, T& left, T& right)> down_functor;

	std::vector<T> v;
	up_functor f_up;     // for up-propagation on update
	down_functor f_down; // for down-propagation for lazy updating

	interval_tree(size_t sz, up_functor f_up, down_functor f_down, T id = T()) : f_up(f_up), f_down(f_down) {
		v.resize(calc_pow2(sz) * 2, id);
		//rebuild(); // must be called manually
	}

	template<typename It>
	interval_tree(It begin, It end, up_functor f_up, down_functor f_down, T id = T()) : f_up(f_up), f_down(f_down) {
		auto sz = std::distance(begin, end);
		v.resize(calc_pow2(sz) * 2, id);
		std::copy(begin, end, v.begin() + size());
		rebuild();
	}

	T get(size_t begin, size_t end) {
		propagate_down(begin, end);
		T tl = v[0], tr = v[0]; // id
		size_t b = begin, e = end, i = size();
		while (b < e) {
			if (b & 1) f_up(tl, tl, v[i + b++]);
			if (e & 1) f_up(tr, v[i + --e], tr);
			b /= 2, e /= 2, i /= 2;
		}
		f_up(tl, tl, tr);
		return tl;
	}

	void update(size_t begin, size_t end, update_functor f) {
		propagate_down(begin, end);
		size_t b = begin, e = end, i = size();
		while (b < e) {
			if (b & 1) f(v[i + b++]);
			if (e & 1) f(v[i + --e]);
			b /= 2, e /= 2, i /= 2;
		}
		propagate_up(begin, end);
	}

	// `rebuild` must be called after all modifications are made.
	T& operator[] (size_t index) {
		index += size();
		return v[index];
	}

	void rebuild() {
		for (size_t i = size() - 1; i > 0; i--) {
			update_up(i);
		}
	}

	size_t size() {
		return v.size() / 2;
	}

private:
	void propagate_down(size_t begin, size_t end) {
		update_from_root(top(begin));
		update_from_root(top(end) - 1);
	}

	void update_from_root(size_t i) {
		for (int h = calc_height(i); h >= 1; h--) {
			update_down(i >> h);
		}
	}

	void update_down(size_t i) {
		f_down(v[i], v[2 * i + 0], v[2 * i + 1]);
	}

	void propagate_up(size_t begin, size_t end) {
		update_to_root(top(begin));
		update_to_root(top(end) - 1);
	}

	void update_to_root(size_t i) {
		while ((i /= 2) > 0) {
			update_up(i);
		}
	}

	void update_up(size_t i) {
		f_up(v[i], v[2 * i + 0], v[2 * i + 1]);
	}

	size_t top(size_t begin) {
		size_t i = size() + begin;
		while (i % 2 == 0) i /= 2;
		return i;
	}

	static int calc_height(size_t sz) {
		int h = 0; while (sz >= (size_t)1 << h) h++;
		return h - 1;
	}

	static size_t calc_pow2(size_t sz) {
		size_t w = 1; while (w < sz) w *= 2;
		return w;
	}
};

} // container
} // altruct
