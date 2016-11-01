#pragma once

#include "algorithm/math/bits.h"

#include <string>
#include <vector>
#include <iterator>
#include <functional>

namespace altruct {
namespace container {

/**
 * An efficient bit-vector implementation.
 *
 * Note: vector<bool> provides something similar,
 * but lacks many important operations.
 */
template<typename W = uint64_t>
class bit_vector {
public:
	typedef W word_t;
	static const W W1 = 1;
	static const W W0 = 0;
	static const W WW = ~W0;
	static const int L = math::bit_size<W>::value;

	// a proxy used as a bit getter/setter
	class bit_proxy {
	public:
		size_t pos;
		bit_vector* owner;
		bit_proxy(bit_vector* owner, size_t pos) : pos(pos), owner(owner) {}
		operator int() { return owner->bit_at(pos); }
		bit_proxy& operator = (int val) { owner->set(pos, val & 1); return *this; }
		bit_proxy& operator ^= (int val) { if (val == 1) owner->set(pos, -1); return *this; }
		bit_proxy& operator |= (int val) { if (val == 1) owner->set(pos, 1); return *this; }
		bit_proxy& operator &= (int val) { if (val == 0) owner->set(pos, 0); return *this; }
	};

	size_t sz;
	std::vector<W> words;

	// constructs a new bit-vector initialized to all zeros
	bit_vector(size_t sz = 0) : sz(sz), words(num_words(sz) + 1) {}

	// constructs a new bit-vector from the elements in the given iterator range
	template<typename It>
	bit_vector(It begin, It end) {
		size_t pos = 0;
		resize(std::distance(begin, end));
		for (auto it = begin; it != end; ++it) {
			set(pos++, (*it) & 1);
		}
	}

	// constructs a new bit-vector from an initializer list
	template<typename I>
	bit_vector(std::initializer_list<I> list) : bit_vector(list.begin(), list.end()) {}

	// returns a clone of this vector
	bit_vector clone() const {
		return bit_vector(*this);
	}

	// comparison operators
	bool operator == (const bit_vector& that) const { return compare(*this, 0, sz, that, 0, that.sz) == 0; }
	bool operator != (const bit_vector& that) const { return compare(*this, 0, sz, that, 0, that.sz) != 0; }
	bool operator <  (const bit_vector& that) const { return compare(*this, 0, sz, that, 0, that.sz) <  0; }
	bool operator >  (const bit_vector& that) const { return compare(*this, 0, sz, that, 0, that.sz) >  0; }
	bool operator <= (const bit_vector& that) const { return compare(*this, 0, sz, that, 0, that.sz) <= 0; }
	bool operator >= (const bit_vector& that) const { return compare(*this, 0, sz, that, 0, that.sz) >= 0; }

	// logical operators
	bit_vector operator & (const bit_vector& that) const { return clone() &= that; }
	bit_vector operator | (const bit_vector& that) const { return clone() |= that; }
	bit_vector operator ^ (const bit_vector& that) const { return clone() ^= that; }

	bit_vector operator ~ () const { return clone().apply(0, sz, op_flip); }

	bit_vector& operator &= (const bit_vector& that) { reserve(that.sz); return apply(*this, 0, that, 0, that.sz, op_and), apply(that.sz, sz, op_set0); }
	bit_vector& operator |= (const bit_vector& that) { reserve(that.sz); return apply(*this, 0, that, 0, that.sz, op_or); }
	bit_vector& operator ^= (const bit_vector& that) { reserve(that.sz); return apply(*this, 0, that, 0, that.sz, op_xor); }

	// scans the range from left to right
	// i.e.: `visitor(this[begin:end])`
	// `bool visitor(W w, size_t pos, int l)` operates on words;
	template<typename F>
	bool scan(size_t begin, size_t end, F visitor) const {
		size_t ib = begin / L, ie = end / L;
		int lb = begin % L, le = end % L;
		bool r = visitor(words[ib] >> lb, 0, ((ib != ie) ? L : le) - lb);
		for (size_t i = ib + 1; r && i < ie; i++) r &= visitor(words[i], i * L - begin, L);
		if (r && ib != ie) r &= visitor(words[ie], end - begin - le, le);
		return r;
	}

	// scans the two ranges from left to right
	// i.e.: `visitor(this[begin1:end1], that[begin2:end2])`
	// `bool visitor(W w1, W w2, int l)` operates on words;
	// it is allowed for `bv` to be `this`;
	template<typename F>
	static bool scan(const bit_vector& v1, size_t begin1, const bit_vector& v2, size_t begin2, size_t len, F visitor) {
		auto visitor2 = [&](W w, size_t pos, int l){
			return visitor(w, v2.word_at(begin2 + pos), l);
		};
		return v1.scan(begin1, begin1 + len, visitor2);
	}

	// applies the operation on the range
	// i.e.: `this[begin:end] = op(this[begin:end])`
	// `op(W w, size_t pos, int len)` operates on words;
	template<typename F>
	bit_vector& apply(size_t begin, size_t end, F op, bool backwards = false) {
		size_t ib = begin / L, ie = end / L;
		int lb = begin % L, le = end % L;
		// calculate the boundary words
		W mb = WW << lb, me = first_bits(le);
		if (ib == ie) mb = me = mb & me;
		W wb = (words[ib] & ~mb) | ((op(words[ib] >> lb, 0, L - lb) << lb) & mb);
		if (ib == ie) { words[ib] = wb; return *this; }
		W we = (words[ie] & ~me) | (op(words[ie], end - begin - le, le) & me);
		// perform the operation on the whole words
		if (backwards) {
			for (size_t i = ie - 1; i > ib; i--) words[i] = op(words[i], i * L - begin, L);
		} else {
			for (size_t i = ib + 1; i < ie; i++) words[i] = op(words[i], i * L - begin, L);
		}
		// set the boundary words
		words[ib] = wb;
		words[ie] = we;
		return *this;
	}

	// applies the operation on the two ranges and stores the result to the first range
	// i.e.: `this[begin1:end1] = op(this[begin1:end1], that[begin2:end2])`
	// `op(W w1, W w2, int len)` operates on words;
	// it is allowed for `bv` to be `this`;
	template<typename F>
	static bit_vector& apply(bit_vector& v1, size_t begin1, const bit_vector& v2, size_t begin2, size_t len, F op) {
		auto op2 = [&](W w, size_t pos, int len){
			return op(w, v2.word_at(begin2 + pos), len);
		};
		// `bv` can be `this`, so loop forward / backwards accordingly
		return v1.apply(begin1, begin1 + len, op2, begin2 < begin1);
	}

	// reverses the given range (position-wise)
	bit_vector& reverse(size_t begin, size_t end) {
		if (begin >= end) return *this;
		size_t ib = begin / L, ie = end / L;
		int lb = begin % L, le = end % L, l = lb + le;
		// calculate the boundary words
		W mb = WW << lb, me = first_bits(le);
		if (ib == ie) mb = me = mb & me;
		W we = words[ie] & ~me;
		if (le != 0) we |= (math::bit_reverse(word_at(begin)) >> (L - le)) & me;
		if (ib == ie) { words[ie] = we; return *this; }
		W wb = words[ib] & ~mb;
		wb |= math::bit_reverse(word_at(end - (L - lb))) & mb;
		// reverse bits within the words
		if (l < L) {
			for (size_t i = ie - 1; i > ib; i--) words[i] = math::bit_reverse(word_at(i - 1, l));
		} else {
			for (size_t i = ib + 1; i < ie; i++) words[i] = math::bit_reverse(word_at(i, l - L));
		}
		// reverse whole words
		for (size_t i = ib + 1, j = ie - 1; i < j; i++, j--) std::swap(words[i], words[j]);
		// set the boundary words
		words[ib] = wb;
		words[ie] = we;
		return *this;
	}

	// swaps the two ranges
	// result is undefined if the two ranges overlap
	bit_vector& swap(size_t begin1, size_t end1, size_t begin2, size_t end2) {
		if (begin1 > begin2) return swap(begin2, end2, begin1, end1);
		reverse(begin1, end1); reverse(end1, begin2); reverse(begin2, end2); reverse(begin1, end2);
		return *this;
	}

	// rotates the given range left
	bit_vector& rotate_left(size_t begin, size_t end, size_t cnt) {
		if (begin >= end) return *this;
		size_t mid = begin + cnt % (end - begin);
		reverse(begin, mid); reverse(mid, end); reverse(begin, end);
		return *this;
	}

	// rotates the given range right
	bit_vector& rotate_right(size_t begin, size_t end, size_t cnt) {
		if (begin >= end) return *this;
		size_t mid = end - cnt % (end - begin);
		reverse(begin, mid); reverse(mid, end); reverse(begin, end);
		return *this;
	}

	// returns the number of bits that differ in the two ranges
	static size_t hamming_distance(const bit_vector& v1, size_t begin1, const bit_vector& v2, size_t begin2, size_t len) {
		size_t dist = 0;
		scan(v1, begin1, v2, begin2, len, [&](W w1, W w2, int l) {
			W w = w1 ^ w2;
			if (l < L) w &= first_bits(l);
			dist += math::bit_cnt1(w);
			return true;
		});
		return dist;
	}

	// copares the two ranges
	static int compare(const bit_vector& v1, size_t begin1, size_t end1, const bit_vector& v2, size_t begin2, size_t end2) {
		size_t len1 = end1 - begin1, len2 = end2 - begin2, len = std::min(len1, len2);
		int ret = 0;
		scan(v1, begin1, v2, begin2, len, [&](W w1, W w2, int l) {
			if (l < L) w1 &= first_bits(l);
			if (l < L) w2 &= first_bits(l);
			if (w1 != w2) ret = (math::bit_reverse(w1) < math::bit_reverse(w2)) ? -1 : +1;
			return w1 == w2;
		});
		if (ret != 0) return ret;
		if (len1 != len2) return (len1 < len2) ? -1 : +1;
		return 0;
	}

	// converts the given range to a string
	std::string to_string(size_t begin, size_t end) const {
		std::string s; s.reserve(end - begin + L);
		for (size_t pos = begin; pos < end; pos += L) {
			W w = word_at(pos);
			for (int k = 0; k < L; k++) {
				s.push_back('0' + (w & 1)), w >>= 1;
			}
		}
		s.resize(end - begin);
		return s;
	}

	// sets the bit at the given bit-position
	void set(size_t pos, int val) {
		size_t i = pos / L;
		int l = pos % L;
		W m = W1 << l;
		if (val == 0) {
			words[i] &= ~m;
		} else if (val == 1) {
			words[i] |= m;
		} else if (val == -1) {
			words[i] ^= m;
		}
	}

	// returns the bit at the given bit-position
	int bit_at(size_t pos) const {
		size_t i = pos / L; int l = pos % L;
		return int((words[i] >> l) & 1);
	}

	// returns the word at the given bit-position
	W word_at(size_t pos) const {
		return word_at(pos / L, pos % L);
	}

	// returns the word at the given word-index and offset within it
	W word_at(size_t i, int l) const {
		return (l == 0) ? words[i] : ((words[i] >> l) | (words[i + 1] << (L - l)));
	}

	// returns the getter/setter bit proxy
	bit_proxy operator[] (size_t pos) {
		return bit_proxy(this, pos);
	}

	// pushes the new element to the end
	void push_back(int val) {
		resize(sz + 1);
		set(sz - 1, val);
	}

	// pops the last element from the end
	void pop_back() {
		resize(sz - 1);
	}

	// reserves space for at least new_size bits
	void reserve(size_t new_size) {
		if (sz < new_size) resize(new_size);
	}

	// resizes the container to the new size
	void resize(size_t new_size) {
		sz = new_size;
		words.resize(num_words(sz) + 1);
		words[sz / L] &= first_bits(sz % L);
		words.back() = W0;
	}

	// returns the number of bits in this bit-array
	size_t size() const {
		return sz;
	}

	// returns the number of words required to store `sz` bits
	static size_t num_words(size_t sz) {
		return (sz + L - 1) / L;
	}

	// returns the word with first l` bits set; `0 <= l < L`
	static W first_bits(int l) {
		return (W1 << l) - 1;
	}

	// operations
	static W op_set0(W w, size_t pos, int len){ return W0; }
	static W op_set1(W w, size_t pos, int len){ return WW; }
	static W op_nop(W w, size_t pos, int len){ return w; }
	static W op_flip(W w, size_t pos, int len){ return ~w; }
	static W op_set(W w1, W w2, int len){ return w2; }
	static W op_setn(W w1, W w2, int len){ return ~w2; }
	static W op_and(W w1, W w2, int len){ return w1 & w2; }
	static W op_or(W w1, W w2, int len){ return w1 | w2; }
	static W op_xor(W w1, W w2, int len){ return w1 ^ w2; }
	static W op_nand(W w1, W w2, int len){ return ~(w1 & w2); }
	static W op_nor(W w1, W w2, int len){ return ~(w1 | w2); }
	static W op_eq(W w1, W w2, int len){ return ~(w1 ^ w2); }
};

} // container
} // altruct
