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

	// TODO: add shifts
	// TODO: add operators
	//	bit_vector ^ bit_vector
	//	bit_vector | bit_vector
	//	bit_vector & bit_vector
	//  comparison
	// TODO: unit tests

	// a proxy used as a bit getter/setter
	class bit_proxy {
	public:
		size_t pos;
		bit_vector* owner;
		bit_proxy(bit_vector* owner, size_t pos) : owner(owner), pos(pos) {}
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
		for (auto it = begin; it != end; i++) {
			set(pos++, *it);
		}
	}

	// sets the bits in the given range; 0=clear, 1=set, -1=flip
	void set(size_t begin, size_t end, int val) {
		size_t ib = begin / L, ie = end / L;
		W mb = WW << (begin % L), me = (W1 << (end % L)) - 1;
		if (ib == ie) mb = me = mb & me;
		if (val == 0) {
			if (ib != ie) words[ib] &= ~mb;
			for (size_t i = ib + 1; i < ie; i++) words[i] = W0;
			words[ie] &= ~me;
		} else if (val == 1) {
			if (ib != ie) words[ib] |= mb;
			for (size_t i = ib + 1; i < ie; i++) words[i] = WW;
			words[ie] |= me;
		} else if (val == -1) {
			if (ib != ie) words[ib] ^= mb;
			for (size_t i = ib + 1; i < ie; i++) words[i] ^= WW;
			words[ie] ^= me;
		}
	}

	// flips the bits in the given range
	void flip(size_t begin, size_t end) {
		set(begin, end, -1);
	}

	// reverses the given range (position-wise)
	void reverse(size_t begin, size_t end) {
		size_t ib = begin / L, ie = end / L;
		size_t lb = begin % L, le = end % L;
		// calculate the boundary words
		W mb = WW << (begin % L), me = (W1 << (end % L)) - 1;
		if (ib == ie) mb = me = mb & me;
		W we = words[ie] & ~me;
		if (le != 0) we |= math::bit_reverse(word_at(begin) << (L - le)) & me;
		if (ib == ie) { words[ie] = we; return; }
		W wb = words[ib] & ~mb;
		wb |= math::bit_reverse(word_at(end - (L - lb))) & mb;
		// reverse bits within the words
		int rot_cnt = int(le - (L - lb));
		if (rot_cnt < 0) {
			for (size_t i = ie - 1; i > ib; i--) words[i] = math::bit_reverse(word_at(i * L + rot_cnt));
		} else {
			for (size_t i = ib + 1; i < ie; i++) words[i] = math::bit_reverse(word_at(i * L + rot_cnt));
		}
		// reverse whole words
		for (size_t i = ib + 1, j = ie - 1; i < j; i++, j--) std::swap(words[i], words[j]);
		// set the boundary words
		words[ib] = wb;
		words[ie] = we;
	}
	
	// result is undefined if the two ranges overlap
	void swap(size_t begin1, size_t end1, size_t begin2, size_t end2) {
		if (begin1 > begin2) return swap(begin2, end2, begin1, end1);
		reverse(begin1, end1); reverse(end1, begin2); reverse(begin2, end2); reverse(begin1, end2);
	}

	// rotates the given range left
	void rotate_left(size_t begin, size_t end, size_t cnt) {
		size_t mid = begin + cnt % (end - begin);
		reverse(begin, mid); reverse(mid, end); reverse(begin, end);
	}

	// rotates the given range right
	void rotate_right(size_t begin, size_t end, size_t cnt) {
		size_t mid = end - cnt % (end - begin);
		reverse(begin, mid); reverse(mid, end); reverse(begin, end);
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

	// returns the number of bits that differ in the two ranges
	size_t hamming_distance(size_t begin1, size_t begin2, size_t len) const {
		size_t dist = 0;
		for (size_t pos = 0; pos + L <= len; pos += L) {
			W w1 = word_at(begin1 + pos);
			W w2 = word_at(begin2 + pos);
			dist += math::bit_cnt1(w1 ^ w2);
		}
		size_t l = len % L;
		if (l > 0) {
			W m = WW >> (L - l);
			W w1 = word_at(begin1 + len - l) & m;
			W w2 = word_at(begin2 + len - l) & m;
			dist += math::bit_cnt1(w1 ^ w2);
		}
		return dist;
	}

	// sets the bit at the given bit-position
	void set(size_t pos, int val) {
		size_t i = pos / L, l = pos % L;
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
		size_t i = pos / L, l = pos % L;
		return int((words[i] >> l) & 1);
	}

	// returns the word at the given bit-position
	W word_at(size_t pos) const {
		size_t i = pos / L, l = pos % L;
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

	// resizes the container to the new size
	void resize(size_t new_size) {
		sz = new_size;
		words.resize(num_words(sz) + 1);
	}

	// returns the number of bits in this bit-array
	size_t size() const {
		return sz;
	}

	// returns the number of words required to store `sz` bits
	static size_t num_words(size_t sz) {
		return (sz + L - 1) / L;
	}
};

} // container
} // altruct
