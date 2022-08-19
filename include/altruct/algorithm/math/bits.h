#pragma once

#include <cmath>
#include <stdint.h>
#include <type_traits>

namespace altruct {
namespace math {

// TODO: Use compiler built-ins when available:
//   GCC: __builtin_clz, __builtin_ctz, __builtin_popcount
//   MSVS: _BitScanForward, _BitScanReverse

/**
 * Size in bits of the given type.
 */
template<typename T>
struct bit_size {
    static const int value = sizeof(T) * 8;
};
template<typename T>
const int bit_size<T>::value;


/** Makes the bit at the given position. */
template<typename T> T make_bit(int pos) { return T(1) << pos; }
/** Makes `cnt` ones in a row. */
template<typename T> T make_ones(int cnt) { return make_bit<T>(cnt) - 1; }
/** Gets the bit at the given position. */
template<typename T> T get_bit(T val, int pos) { return (val >> pos) & 1; }
/** Sets the bit at the given position. */
template<typename T> T set_bit(T val, int pos) { return val | (T(1) << pos); }
/** Flips the bit at the given position. */
template<typename T> T flip_bit(T val, int pos) { return val ^ (T(1) << pos); }
/** Clearsthe bit at the given position. */
template<typename T> T clear_bit(T val, int pos) { return val & ~(T(1) << pos); }
/** Erases the bit at the given position. The higher bits get shifted right. */
template<typename T> T erase_bit(T val, int pos) {
    T lo_mask = (T(1) << pos) - 1;
    return ((val >> 1) & ~lo_mask) | (val & lo_mask);
}
/** If the corresponding mask bit is 0, bit from val0 is used, otherwise bit from val1 is used. */
template<typename T> T mix_bits(T val0, T val1, T mask) { return (val0 & ~mask) | (val1 & mask); }


/**
 * Base-2 Logarithm.
 * Note: `ilog2(0) = 0` for simpler implementation.
 */
int ilog2(uint8_t x);
int ilog2(uint16_t x);
int ilog2(uint32_t x);
int ilog2(uint64_t x);
int ilog2(int8_t x);
int ilog2(int16_t x);
int ilog2(int32_t x);
int ilog2(int64_t x);

/**
 * Number of bits set to 1.
 */
int bit_cnt1(uint8_t x);
int bit_cnt1(uint16_t x);
int bit_cnt1(uint32_t x);
int bit_cnt1(uint64_t x);

/**
 * Reverse bits (position-wise).
 */
uint8_t bit_reverse(uint8_t x);
uint16_t bit_reverse(uint16_t x);
uint32_t bit_reverse(uint32_t x);
uint64_t bit_reverse(uint64_t x);

/**
 * or_down recursive template.
 *
 *  x |= (x >> 32);
 *  x |= (x >> 16);
 *  x |= (x >> 8);
 *  x |= (x >> 4);
 *  x |= (x >> 2);
 *  x |= (x >> 1);
 */
template<typename I, int SHIFT>
struct _or_down {
    static I of(I x) {
        return _or_down<I, SHIFT / 2>::of(x | (x >> SHIFT));
    }
};
template<typename I>
struct _or_down<I, 0> {
    static I of(I x) {
        return x;
    }
};
template<typename I>
I or_down(I x) {
    return _or_down<I, bit_size<I>::value / 2>::of(x);
}

/**
 * xor_down recursive template.
 *
 *  x ^= (x >> 32);
 *  x ^= (x >> 16);
 *  x ^= (x >> 8);
 *  x ^= (x >> 4);
 *  x ^= (x >> 2);
 *  x ^= (x >> 1);
 */
template<typename I, int SHIFT>
struct _xor_down {
    static I of(I x) {
        return _xor_down<I, SHIFT / 2>::of(x ^ (x >> SHIFT));
    }
};
template<typename I>
struct _xor_down<I, 0> {
    static I of(I x) {
        return x;
    }
};
template<typename I>
I xor_down(I x) {
    return _xor_down<I, bit_size<I>::value / 2>::of(x);
}

/**
 * Performs two's complement negation.
 * Without a compiler warning for unsigned types.
 */
template<typename I>
I neg(I x) {
    return I(~x) + I(1);
}

/**
 * Gray-code to binary number conversion.
 */
template<typename I>
I gray_to_bin(I x) {
    return xor_down(x);
}

/**
 * Binary number to Gray-code conversion.
 */
template<typename I>
I bin_to_gray(I x) {
    return (x ^ (x >> 1));
}

/**
 * Leaves only the highest bit set.
 */
template<typename I>
I hi_bit(I x) {
    x = or_down(x);
    return (x ^ (x >> 1));
}

/**
 * Leaves only the lowest bit set.
 */
template<typename I>
I lo_bit(I x) {
    return (x & neg(x));
}

/**
 * Whether the number is not power of two.
 * Note: `0` is considered to be a power of two.
 */
template<typename I>
bool is_not_pow2(I x) {
    return (x & (x - 1)) != 0;
}

/**
 * Whether the number is power of two.
 * Note: `0` is considered to be a power of two.
 */
template<typename I>
bool is_pow2(I x) {
    return (x & (x - 1)) == 0;
}

/**
 * The smallest power of two bigger than `x`.
 * Note: `0` is considered to be a power of two.
 */
template<typename I>
I next_pow2(I x) {
    return I(or_down(x) + 1);
}

/**
 * Leading zeros count. (Zeros from MSB)
 */
template<typename I>
int lzc(I x) {
    return bit_size<I>::value - bit_cnt1(or_down(x));
}

/**
 * Trailing zeros count. (Zeros from LSB)
 */
template<typename I>
int tzc(I x) {
    return bit_cnt1(I(lo_bit(x) - 1));
}

/**
 * Two's complement <==> Sign & Magnitude.
 * The conversion procedure is same both ways.
 */
template<typename I>
I sign_mag(I x) {
    const I HI_BIT = I(1) << (bit_size<I>::value - 1);
    return (x & HI_BIT) ? neg(x) ^ HI_BIT : x;
}

/**
 * Next lexicographic combination with the same number of ones.
 * Modifies `x` and returns false if last combination was passed.
 */
template<typename I>
bool next_combination(I& x, int size) {
    I v = (x + lo_bit(x)) & make_ones<I>(size);
    I w = (v ^ x) >> tzc(x);
    I y = (v | (w >> 2));
    x = (v != 0) ? y : w;
    return (v != 0);
}

} // math
} // altruct
