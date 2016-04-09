#include "algorithm/random/random.h"
#include "algorithm/random/xorshift.h"

namespace altruct {
namespace random {

xorshift_1024star::xorshift_1024star() { seed(uint64_t(0)); }

xorshift_1024star::xorshift_1024star(uint64_t state) { seed(state); }

xorshift_1024star::xorshift_1024star(uint64_t* state) { seed(state); }

void xorshift_1024star::seed(uint64_t state) {
	xorshift_64star xs64(state);
	for (int i = 0; i < 16; i++) {
		s_[i] = xs64.next();
	}
	p_ = 0;
}

void xorshift_1024star::seed(uint64_t* state) {
	for (int i = 0; i < 16; i++) {
		s_[i] = state[i];
	}
	p_ = 0;
}

uint64_t xorshift_1024star::next() {
	// The constants used in this implementation are as suggested by the author.
	uint64_t s0 = s_[p_];
	p_ = (p_ + 1) & 15;
	uint64_t s1 = s_[p_];
	s1 ^= s1 << 31; // a
	s1 ^= s1 >> 11; // b
	s0 ^= s0 >> 30; // c
	s_[p_] = s0 ^ s1;
	return s_[p_] * 1181783497276652981ULL;
}

uint64_t xorshift_1024star::next(uint64_t min, uint64_t max) {
	return integer_to_range<uint64_t>(next(), min, max);
}

uint64_t xorshift_1024star::next_uniform(uint64_t min, uint64_t max) {
	return uniform_next<uint64_t>([&](){ return next(); }, min, max);
}

double xorshift_1024star::next_0_1() {
	return integer_to_double_0_1<uint64_t>(next());
}

xorshift_64star::xorshift_64star() { seed(uint64_t(0)); }

xorshift_64star::xorshift_64star(uint64_t state) { seed(state); }

void xorshift_64star::seed(uint64_t state) { x_ = state; }

uint64_t xorshift_64star::next() {
	// The constants used in this implementation are as suggested by the author.
	x_ ^= x_ >> 12; // a
	x_ ^= x_ << 25; // b
	x_ ^= x_ >> 27; // c
	return x_ * 2685821657736338717ULL;
}

} // altruct
} // random
