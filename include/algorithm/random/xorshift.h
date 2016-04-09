#pragma once

#include <cstdint>

namespace altruct {
namespace random {

/**
 * This is an adaptation of the public domain xorshift1024star generator
 * from http://xorshift.di.unimi.it/xorshift1024star.c
 */
class xorshift_1024star {
public:
	/**
	 * Constructs a new instance of this class.
	 * The state must be initialized by calling {@code seed} before consuming
	 * the values produced by {@code next}.
	 */
	xorshift_1024star();

	/**
	 * Constructs a new instance of this class seeded with the provided state.
	 * The values produced by {@code next} are ready to be consumed.
	 */
	explicit xorshift_1024star(uint64_t state);

	/**
	 * Constructs a new instance of this class seeded with the provided state.
	 * The provided state should consist of <b>exactly</b> 16 uint64_t values.
	 * The values produced by {@code next} are ready to be consumed.
	 */
	explicit xorshift_1024star(uint64_t* state);

	virtual ~xorshift_1024star() {}

	/**
	 * Seeds this instance with the provided 64bit state.
	 * The provided 64bit state will be used as a seed for xorshift_64star
	 * generator which will then be used to generate the 16 uint64_t values
	 * for the new state of this generator.
	 */
	virtual void seed(uint64_t state64);

	/**
	 * Seeds this instance with the provided state.
	 * The provided state should consist of <b>exactly</b> 16 uint64_t values.
	 */
	virtual void seed(uint64_t* state);

	/**
	 * Gets the next random number.
	 */
	virtual uint64_t next();

	/**
	 * Gets the next random number in [min, max] range, both inclusive.
	 */
	virtual uint64_t next(uint64_t min, uint64_t max);

	/*
	 * Gets the next random number in range [min, max] inclusive, with stronger
	 * uniformity guarantees at expense of decreased performance.
	 * In most cases {@code next(min, max)} will suffice.
	 */
	virtual uint64_t next_uniform(uint64_t min, uint64_t max);

	/**
	 * Gets the next random number as a double in [0, 1] range, both inclusive.
	 */
	virtual double next_0_1();

private:
	/**
	 * The state must be seeded so that it is not everywhere zero. If you have
	 * a 64-bit seed, we suggest to seed a xorshift64* generator and use its
	 * output to fill the state s.
	 */
	uint64_t s_[16];
	int p_;
};

/**
 * A 64 bit variant of the xorshift* pseudo-random number generator family.
 */
class xorshift_64star {
public:
	/**
	 * Constructs a new instance of this class.
	 * The state must be initialized by calling {@code seed} before consuming
	 * the values produced by {@code next}.
	 */
	xorshift_64star();

	/**
	 * Constructs a new instance of this class seeded with the provided state.
	 * The values produced by {@code next} are ready to be consumed.
	 */
	explicit xorshift_64star(uint64_t state);

	virtual ~xorshift_64star() {}

	/**
	 * Seeds this instance with the provided 64bit state.
	 */
	virtual void seed(uint64_t state);

	/**
	 * Gets the next random number.
	 */
	virtual uint64_t next();

private:
	/* The state. Must be seeded with a nonzero value. */
	uint64_t x_;
};

} // random
} // altruct
