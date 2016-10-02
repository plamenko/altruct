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
template<size_t K, typename I = int32_t, typename IT = int64_t>
class polynomial_hash {
public:
	static I M[K];  // modulos
	static I B[K];  // bases
	static I BI[K]; // base inverses; B * BI == 1 (mod M)

	static std::vector<I> W[K];
	static std::vector<I> WI[K];

	static I mul(I a, I b, I m) { return (I)((a * (IT)b) % m); }
	static I mul(I s, I a, I b, I m) { return (I)((s + a * (IT)b) % m); }

	static void _ensure(std::vector<I>& w, size_t sz, I b, I m) {
		size_t i0 = w.size();
		w.resize(sz);
		if (i0 == 0) {
			w[i0++] = 1;
		}
		for (size_t i = i0; i < sz; i++) {
			w[i] = mul(w[i - 1], b, m);
		}
	}

	static void ensure(size_t size) {
		size_t curr_size = W[0].size();
		if (curr_size >= size) return;
		size = std::max(size, curr_size + curr_size / 2);
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

	polynomial_hash clone() const {
		return polynomial_hash(*this);
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

	polynomial_hash operator * (I rhs) const { return clone() *= rhs; }
	polynomial_hash& operator *= (I rhs) {
		for (int k = 0; k < K; k++) {
			h[k] = mul(h[k], rhs, M[k]);
		}
		return *this;
	}

	polynomial_hash operator * (const polynomial_hash& rhs) const { return clone() *= rhs; }
	polynomial_hash& operator *= (const polynomial_hash& rhs) {
		for (int k = 0; k < K; k++) {
			h[k] = mul(h[k], rhs.h[k], M[k]);
		}
		return *this;
	}

	polynomial_hash operator + (I rhs) const { return clone() += rhs; }
	polynomial_hash& operator += (I rhs) {
		for (int k = 0; k < K; k++) {
			h[k] = (h[k] + rhs) % M[k];
		}
		return *this;
	}

	polynomial_hash operator + (const polynomial_hash& rhs) const { return clone() += rhs; }
	polynomial_hash& operator += (const polynomial_hash& rhs) {
		for (int k = 0; k < K; k++) {
			h[k] = (h[k] + rhs.h[k]) % M[k];
		}
		return *this;
	}

	polynomial_hash operator - (I rhs) const { return clone() -= rhs; }
	polynomial_hash& operator -= (I rhs) {
		for (int k = 0; k < K; k++) {
			h[k] = (h[k] + M[k] - rhs) % M[k];
		}
		return *this;
	}

	polynomial_hash operator - (const polynomial_hash& rhs) const { return clone() -= rhs; }
	polynomial_hash& operator -= (const polynomial_hash& rhs) {
		for (int k = 0; k < K; k++) {
			h[k] = (h[k] + M[k] - rhs.h[k]) % M[k];
		}
		return *this;
	}

	// divide by B^cnt
	polynomial_hash operator >> (int cnt) const { return clone() >>= cnt; }
	polynomial_hash& operator >>= (int cnt) {
		ensure(cnt + 1);
		for (int k = 0; k < K; k++) {
			h[k] = mul(h[k], WI[k][cnt], M[k]);
		}
		return *this;
	}

	// multiply by B^cnt
	polynomial_hash operator << (int cnt) const { return clone() <<= cnt; }
	polynomial_hash& operator <<= (int cnt) {
		ensure(cnt + 1);
		for (int k = 0; k < K; k++) {
			h[k] = mul(h[k], W[k][cnt], M[k]);
		}
		return *this;
	}


	// H = H + (RHS << pos)
	polynomial_hash& add(I rhs, size_t pos) {
		ensure(pos + 1);
		for (int k = 0; k < K; k++) {
			h[k] = mul(h[k], rhs, W[k][pos], M[k]);
		}
		return *this;
	}

	// H = H + (RHS << pos)
	polynomial_hash& add(const polynomial_hash& rhs, size_t pos) {
		ensure(pos + 1);
		for (int k = 0; k < K; k++) {
			h[k] = mul(h[k], rhs.h[k], W[k][pos], M[k]);
		}
		return *this;
	}

	// H = (H - RHS) >> pos
	polynomial_hash& sub_shr(I rhs, size_t pos) {
		ensure(pos + 1);
		for (int k = 0; k < K; k++) {
			h[k] = mul(h[k] + M[k] - rhs, WI[k][pos], M[k]);
		}
		return *this;
	}

	// H = (H - RHS) >> pos
	polynomial_hash& sub_shr(const polynomial_hash& rhs, size_t pos) {
		ensure(pos + 1);
		for (int k = 0; k < K; k++) {
			h[k] = mul(h[k] + M[k] - rhs.h[k], WI[k][pos], M[k]);
		}
		return *this;
	}

	size_t hash() const {
		return (size_t)h[0];
	}
};

template<size_t K, typename I, typename IT>
std::vector<I> polynomial_hash<K, I, IT>::W[K];
template<size_t K, typename I, typename IT>
std::vector<I> polynomial_hash<K, I, IT>::WI[K];

/**
 * Cumulative hashes of a sequence (e.g. of a string).
 *
 * @param HASH - the underlying hash type that supports `add` and `sub_shr`
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
		assign(begin, end);
	}
	
	template<typename It>
	void assign(It begin, It end) {
		HASH h;
		size_t pos = 0;
		h.ensure(std::distance(begin, end));
		for (It it = begin; it != end; ++it) {
			h.add(*it, pos++);
			vh.push_back(h);
		}
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

	HASH get(size_t begin, size_t end) const {
		HASH r;
		if (end > 0) r.add(vh.at(end - 1), 0);
		if (begin > 0) r.sub_shr(vh.at(begin - 1), begin);
		return r;
	}

	size_t size() const {
		return vh.size();
	}
};

}
}
