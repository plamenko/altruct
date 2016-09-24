#include "algorithm/math/bits.h"

#include <cmath>
#include <stdint.h>

namespace altruct {
namespace math {

/**
 * Base-2 Logarithm.
 */
namespace {
int8_t _tbl_ilog2_16_full[1 << 16];
int8_t _tbl_ilog2_16[16] = {
	0, 7, 1, 13, 8, 10, 2, 14, 6, 12, 9, 5, 11, 4, 3, 15 };
int8_t _tbl_ilog2_32[32] = {
	0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
	8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31 };
int8_t _tbl_ilog2_64[64] = {
	0, 58, 1, 59, 47, 53, 2, 60, 39, 48, 27, 54, 33, 42, 3, 61,
	51, 37, 40, 49, 18, 28, 20, 55, 30, 34, 11, 43, 14, 22, 4, 62,
	57, 46, 52, 38, 26, 32, 41, 50, 36, 17, 19, 29, 10, 13, 21, 56,
	45, 25, 31, 35, 16, 9, 12, 44, 24, 15, 8, 23, 7, 6, 5, 63 };
}
int ilog2(uint8_t x) { return _tbl_ilog2_16_full[x]; }
int ilog2(uint16_t x) { return _tbl_ilog2_16_full[x]; }
//int ilog2(uint16_t x) { return _tbl_ilog2_16[uint16_t(or_down(x) * (uint16_t)0x0F2D) >> 12]; }
int ilog2(uint32_t x) { return _tbl_ilog2_32[uint32_t(or_down(x) * (uint32_t)0x07C4ACDD) >> 27]; }
int ilog2(uint64_t x) { return _tbl_ilog2_64[uint64_t(or_down(x) * (uint64_t)0x03f6eaf2cd271461) >> 58]; }

/**
 * Number of bits set to 1.
 */
namespace {
uint8_t _tbl_cnt1_16[1 << 16];
}
int bit_cnt1(uint8_t x) { return _tbl_cnt1_16[x]; }
int bit_cnt1(uint16_t x) { return _tbl_cnt1_16[x]; }
int bit_cnt1(uint32_t x) { return bit_cnt1(uint16_t(x)) + bit_cnt1(uint16_t(x >> 16)); }
int bit_cnt1(uint64_t x) { return bit_cnt1(uint32_t(x)) + bit_cnt1(uint32_t(x >> 32)); }

/**
 * Reverse bits (position-wise).
 */
namespace {
uint8_t _tbl_rev_4[1 << 4] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };
uint8_t _tbl_rev_8[1 << 8];
uint16_t _tbl_rev_16[1 << 16];
}
uint8_t bit_reverse(uint8_t x) { return _tbl_rev_8[x]; }
uint16_t bit_reverse(uint16_t x) { return _tbl_rev_16[x]; }
uint32_t bit_reverse(uint32_t x) { return (uint32_t(bit_reverse(uint16_t(x))) << 16) | bit_reverse(uint16_t(x >> 16)); }
uint64_t bit_reverse(uint64_t x) { return (uint64_t(bit_reverse(uint32_t(x))) << 32) | bit_reverse(uint32_t(x >> 32)); }
uint8_t _bit_reverse_8(uint8_t x) { return uint8_t((_tbl_rev_4[x & 0xF] << 4) | _tbl_rev_4[x >> 4]); }
uint16_t _bit_reverse_16(uint16_t x) { return (uint16_t(bit_reverse(uint8_t(x))) << 8) | bit_reverse(uint8_t(x >> 8)); }

/**
 * Global initialization.
 * We could do static local initialization in methods,
 * but that would impose an additional overhead for
 * each method call. Given that the methods here are
 * branch-free and perform a constant number of ALU
 * operations, this overhead would be unacceptable.
 */
bool _altruct_math_bits_init() {
	for (int x = 2; x < (1 << 16); x++) _tbl_ilog2_16_full[x] = _tbl_ilog2_16_full[x / 2] + 1;
	for (int x = 1; x < (1 << 16); x++) _tbl_cnt1_16[x] = _tbl_cnt1_16[x / 2] + (x & 1);
	for (int x = 1; x < (1 << 8); x++) _tbl_rev_8[x] = _bit_reverse_8(uint8_t(x));
	for (int x = 1; x < (1 << 16); x++) _tbl_rev_16[x] = _bit_reverse_16(uint16_t(x));
	return true;
}
bool _altruct_math_bits_initialized = _altruct_math_bits_init();

} // math
} // altruct
