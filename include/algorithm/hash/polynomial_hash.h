#pragma once

#include <stdint.h>
#include <type_traits>
#include <array>
#include <vector>

namespace altruct {
namespace hash {

/**
 * Polynomial hash with `K` bases.
 *
 * `M`, `B`, `BI` must be defined by the client.
 */
template<size_t K, typename I = int64_t>
class polynomial_hash {
public:
	static I M[K];  // modulos
	static I B[K];  // bases
	static I BI[K]; // base inverses; B * BI == 1 (mod M)

	static std::vector<I> W[K];
	static std::vector<I> WI[K];

	static void _ensure(std::vector<I>& w, size_t sz, I b, I m) {
		size_t i0 = w.size();
		w.resize(sz);
		if (i0 == 0) {
			w[i0++] = 1;
		}
		for (size_t i = i0; i < sz; i++) {
			w[i] = (w[i - 1] * b) % m;
		}
	}

	static void _ensure(size_t size) {
		size_t curr_size = W[0].size();
		if (curr_size >= size) return;
		for (int k = 0; k < K; k++) {
			_ensure(W[k], size, B[k], M[k]);
			_ensure(WI[k], size, BI[k], M[k]);
		}
	}


	std::array<I, K> h;

	polynomial_hash() {
		for (int k = 0; k < K; k++) {
			h[k] = 0;
		}
	}

	polynomial_hash(std::initializer_list<I> list) {
		size_t sz = std::min(h.size(), list.size());
		std::copy(list.begin(), list.begin() + sz, h.begin());
		std::fill_n(h.begin() + sz, h.size() - sz, 0);
	}

	bool operator < (const polynomial_hash& rhs) const {
		for (int k = 0; k < K; k++) {
			if (h[k] != rhs.h[k]) return h[k] < rhs.h[k];
		}
		return false;
	}

	bool operator == (const polynomial_hash& rhs) const {
		for (int k = 0; k < K; k++) {
			if (h[k] != rhs.h[k]) return false;
		}
		return true;
	}

	polynomial_hash& add(I rhs, size_t pos) {
		_ensure(pos + 1);
		for (int k = 0; k < K; k++) {
			h[k] = (h[k] + rhs * W[k][pos]) % M[k];
		}
		return *this;
	}

	// H = H + (RHS << pos)
	polynomial_hash& add(const polynomial_hash& rhs, size_t pos) {
		_ensure(pos + 1);
		for (int k = 0; k < K; k++) {
			h[k] = (h[k] + rhs.h[k] * W[k][pos]) % M[k];
		}
		return *this;
	}

	// H = (H - RHS) >> pos
	polynomial_hash& subtract(const polynomial_hash& rhs, size_t pos) {
		_ensure(pos + 1);
		for (int k = 0; k < K; k++) {
			h[k] = (h[k] + M[k] - rhs.h[k]) * WI[k][pos] % M[k];
		}
		return *this;
	}

	size_t hash() const {
		return (size_t)h[0];
	}
};

template<size_t K, typename I>
std::vector<I> polynomial_hash<K, I>::W[K];
template<size_t K, typename I>
std::vector<I> polynomial_hash<K, I>::WI[K];

/**
 * Cumulative hashes of a sequence (e.g. of a string).
 *
 * @param HASH - the underlying hash type that supports `add` and `subtract`
 *               (e.g. `polynomial_hash`).
 * 
 * Space complexity: `O(n)`.
 * Time complexities:
 *   build: `O(n)`
 *   get: `O(1)`
 */
template<typename HASH>
class cumulative_hash {
	std::vector<HASH> vh;

public:
	cumulative_hash() {
	}

	template<typename It>
	cumulative_hash(It begin, It end) {
		HASH h;
		size_t pos = 0;
		for (It it = begin; it != end; ++it) {
			h.add(*it, pos++);
			vh.push_back(h);
		}
	}

	size_t size() {
		return vh.size();
	}

	template<typename I>
	void push_back(I rhs) {
		HASH h = vh.empty() ? HASH() : vh.back();
		h.add(rhs, vh.size());
		vh.push_back(h);
	}

	void pop_back() {
		vh.pop_back();
	}

	HASH get(size_t begin, size_t end) {
		HASH r;
		if (end > 0) r.add(vh.at(end - 1), 0);
		if (begin > 0) r.subtract(vh.at(begin - 1), begin);
		return r;
	}
};

}
}
