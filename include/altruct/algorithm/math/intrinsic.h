#pragma once

#include "altruct/algorithm/math/base.h"

#if defined(__clang__)
// builtin, no header needed
#elif defined(__GNUC__)
// builtin, no header needed
#elif defined(_MSC_VER)
#include <intrin.h>
#endif

namespace altruct {
namespace math {

#if defined(__clang__) || defined(__GNUC__)
template<typename I>
bool add_overflow(I x, I y, I* r) {
    return __builtin_add_overflow(x, y, r);
}
#elif defined(_MSC_VER)
inline bool add_overflow(uint64_t x, uint64_t y, uint64_t* r) { return _addcarry_u64(0, x, y, r); }
inline bool add_overflow(uint32_t x, uint32_t y, uint32_t* r) { return _addcarry_u32(0, x, y, r); }
inline bool add_overflow(uint16_t x, uint16_t y, uint16_t* r) { return _addcarry_u16(0, x, y, r); }
inline bool add_overflow(uint8_t x, uint8_t y, uint8_t* r) { return _addcarry_u8(0, x, y, r); }
//inline bool add_overflow(int64_t x, int64_t y, int64_t* r) { auto c = _addcarry_u64(0, x, y, (uint64_t*)r); auto ci = (x ^ y ^ *r) >> 63; return c ^ ci; }
//inline bool add_overflow(int32_t x, int32_t y, int32_t* r) { auto c = _addcarry_u32(0, x, y, (uint32_t*)r); auto ci = (x ^ y ^ *r) >> 31; return c ^ ci; }
//inline bool add_overflow(int16_t x, int16_t y, int16_t* r) { auto c = _addcarry_u16(0, x, y, (uint16_t*)r); auto ci = (x ^ y ^ *r) >> 15; return c ^ ci; }
//inline bool add_overflow(int8_t x, int8_t y, int8_t* r) { auto c = _addcarry_u8(0, x, y, (uint8_t*)r); auto ci = (x ^ y ^ *r) >> 7; return c ^ ci; }
#endif

} // math
} // altruct
