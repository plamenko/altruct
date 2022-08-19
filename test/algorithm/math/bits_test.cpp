#include "altruct/algorithm/math/bits.h"

#include "gtest/gtest.h"
#include <ctime>

using namespace std;
using namespace altruct::math;

namespace {
uint64_t or_down_64(uint64_t x) {
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);
    x |= (x >> 32);
    return x;
}

uint64_t bit_rev_64(uint64_t x) {
    x = ((x >> 1) & 0x5555555555555555ULL) | ((x & 0x5555555555555555ULL) << 1);
    x = ((x >> 2) & 0x3333333333333333ULL) | ((x & 0x3333333333333333ULL) << 2);
    x = ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((x & 0x0F0F0F0F0F0F0F0FULL) << 4);
    x = ((x >> 8) & 0x00FF00FF00FF00FFULL) | ((x & 0x00FF00FF00FF00FFULL) << 8);
    x = ((x >> 16) & 0x0000FFFF0000FFFFULL) | ((x & 0x0000FFFF0000FFFFULL) << 16);
    x = ((x >> 32)) | ((x) << 32);
    return x;
}

int bit_cnt1_64(uint64_t x) {
    x = (x & 0x5555555555555555ULL) + ((x >> 1) & 0x5555555555555555ULL);
    x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
    x = (x & 0x0F0F0F0F0F0F0F0FULL) + ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL);
    x += x >> 8;
    x += x >> 16;
    x += x >> 32;
    return (int)(x & 0x7F);
}

// ilog2_64(0) = -1
int ilog2_64(uint64_t x) {
    return (bit_cnt1_64(or_down_64(x)) - 1);
}
}

TEST(bits_test, bit_size) {
    ASSERT_EQ(8, bit_size<uint8_t>::value);
    ASSERT_EQ(16, bit_size<uint16_t>::value);
    ASSERT_EQ(32, bit_size<uint32_t>::value);
    ASSERT_EQ(64, bit_size<uint64_t>::value);
    ASSERT_EQ(8, bit_size<int8_t>::value);
    ASSERT_EQ(16, bit_size<int16_t>::value);
    ASSERT_EQ(32, bit_size<int32_t>::value);
    ASSERT_EQ(64, bit_size<int64_t>::value);
}

TEST(bits_test, make_bit) {
    for (int i = 0; i < 64; i++) ASSERT_EQ((uint64_t(1) << i), make_bit<uint64_t>(i));
    for (int i = 0; i < 32; i++) ASSERT_EQ((uint32_t(1) << i), make_bit<uint32_t>(i));
    for (int i = 0; i < 16; i++) ASSERT_EQ((uint16_t(1) << i), make_bit<uint16_t>(i));
    for (int i = 0; i < 8; i++) ASSERT_EQ((uint8_t(1) << i), make_bit<uint8_t>(i));
}

TEST(bits_test, make_ones) {
    for (int i = 0; i < 64; i++) ASSERT_EQ((uint64_t(1) << i) - 1, make_ones<uint64_t>(i));
    for (int i = 0; i < 32; i++) ASSERT_EQ((uint32_t(1) << i) - 1, make_ones<uint32_t>(i));
    for (int i = 0; i < 16; i++) ASSERT_EQ((uint16_t(1) << i) - 1, make_ones<uint16_t>(i));
    for (int i = 0; i < 8; i++) ASSERT_EQ((uint8_t(1) << i) - 1, make_ones<uint8_t>(i));
}

TEST(bits_test, get_bit) {
    uint64_t x = 0x7BD152B330F0A777; uint64_t y = ~x;
    for (int i = 0; i < 64; i++) ASSERT_EQ((x >> i) & 1, get_bit(uint64_t(x), i)) << i;
    for (int i = 0; i < 64; i++) ASSERT_EQ((y >> i) & 1, get_bit(uint64_t(y), i)) << i;
    for (int i = 0; i < 32; i++) ASSERT_EQ((x >> i) & 1, get_bit(uint32_t(x), i)) << i;
    for (int i = 0; i < 32; i++) ASSERT_EQ((y >> i) & 1, get_bit(uint32_t(y), i)) << i;
    for (int i = 0; i < 16; i++) ASSERT_EQ((x >> i) & 1, get_bit(uint16_t(x), i)) << i;
    for (int i = 0; i < 16; i++) ASSERT_EQ((y >> i) & 1, get_bit(uint16_t(y), i)) << i;
    for (int i = 0; i < 8; i++) ASSERT_EQ((x >> i) & 1, get_bit(uint8_t(x), i)) << i;
    for (int i = 0; i < 8; i++) ASSERT_EQ((y >> i) & 1, get_bit(uint8_t(y), i)) << i;
}

TEST(bits_test, set_bit) {
    uint64_t x = 0x7BD152B330F0A777; uint64_t y = ~x;
    for (int i = 0; i < 64; i++) ASSERT_EQ(uint64_t(x) | (uint64_t(1) << i), set_bit(uint64_t(x), i)) << i;
    for (int i = 0; i < 64; i++) ASSERT_EQ(uint64_t(y) | (uint64_t(1) << i), set_bit(uint64_t(y), i)) << i;
    for (int i = 0; i < 32; i++) ASSERT_EQ(uint32_t(x) | (uint32_t(1) << i), set_bit(uint32_t(x), i)) << i;
    for (int i = 0; i < 32; i++) ASSERT_EQ(uint32_t(y) | (uint32_t(1) << i), set_bit(uint32_t(y), i)) << i;
    for (int i = 0; i < 16; i++) ASSERT_EQ(uint16_t(x) | (uint16_t(1) << i), set_bit(uint16_t(x), i)) << i;
    for (int i = 0; i < 16; i++) ASSERT_EQ(uint16_t(y) | (uint16_t(1) << i), set_bit(uint16_t(y), i)) << i;
    for (int i = 0; i < 8; i++) ASSERT_EQ(uint8_t(x) | (uint8_t(1) << i), set_bit(uint8_t(x), i)) << i;
    for (int i = 0; i < 8; i++) ASSERT_EQ(uint8_t(y) | (uint8_t(1) << i), set_bit(uint8_t(y), i)) << i;
}

TEST(bits_test, flip_bit) {
    uint64_t x = 0x7BD152B330F0A777; uint64_t y = ~x;
    for (int i = 0; i < 64; i++) ASSERT_EQ(uint64_t(x) ^ (uint64_t(1) << i), flip_bit(uint64_t(x), i)) << i;
    for (int i = 0; i < 64; i++) ASSERT_EQ(uint64_t(y) ^ (uint64_t(1) << i), flip_bit(uint64_t(y), i)) << i;
    for (int i = 0; i < 32; i++) ASSERT_EQ(uint32_t(x) ^ (uint32_t(1) << i), flip_bit(uint32_t(x), i)) << i;
    for (int i = 0; i < 32; i++) ASSERT_EQ(uint32_t(y) ^ (uint32_t(1) << i), flip_bit(uint32_t(y), i)) << i;
    for (int i = 0; i < 16; i++) ASSERT_EQ(uint16_t(x) ^ (uint16_t(1) << i), flip_bit(uint16_t(x), i)) << i;
    for (int i = 0; i < 16; i++) ASSERT_EQ(uint16_t(y) ^ (uint16_t(1) << i), flip_bit(uint16_t(y), i)) << i;
    for (int i = 0; i < 8; i++) ASSERT_EQ(uint8_t(x) ^ (uint8_t(1) << i), flip_bit(uint8_t(x), i)) << i;
    for (int i = 0; i < 8; i++) ASSERT_EQ(uint8_t(y) ^ (uint8_t(1) << i), flip_bit(uint8_t(y), i)) << i;
}

TEST(bits_test, clear_bit) {
    uint64_t x = 0x7BD152B330F0A777; uint64_t y = ~x;
    for (int i = 0; i < 64; i++) ASSERT_EQ(uint64_t(x) & ~(uint64_t(1) << i), clear_bit(uint64_t(x), i)) << i;
    for (int i = 0; i < 64; i++) ASSERT_EQ(uint64_t(y) & ~(uint64_t(1) << i), clear_bit(uint64_t(y), i)) << i;
    for (int i = 0; i < 32; i++) ASSERT_EQ(uint32_t(x) & ~(uint32_t(1) << i), clear_bit(uint32_t(x), i)) << i;
    for (int i = 0; i < 32; i++) ASSERT_EQ(uint32_t(y) & ~(uint32_t(1) << i), clear_bit(uint32_t(y), i)) << i;
    for (int i = 0; i < 16; i++) ASSERT_EQ(uint16_t(x) & ~(uint16_t(1) << i), clear_bit(uint16_t(x), i)) << i;
    for (int i = 0; i < 16; i++) ASSERT_EQ(uint16_t(y) & ~(uint16_t(1) << i), clear_bit(uint16_t(y), i)) << i;
    for (int i = 0; i < 8; i++) ASSERT_EQ(uint8_t(x) & ~(uint8_t(1) << i), clear_bit(uint8_t(x), i)) << i;
    for (int i = 0; i < 8; i++) ASSERT_EQ(uint8_t(y) & ~(uint8_t(1) << i), clear_bit(uint8_t(y), i)) << i;
}

namespace {
template<typename T>
T erase_slow(T val, int pos) {
    int sz = bit_size<T>::value;
    vector<int> v(sz);
    for (int i = 0; i < sz; i++) {
        v[i] = (val >> i) & 1;
    }
    v.erase(v.begin() + pos);
    val = 0;
    for (int i = 0; i < sz - 1; i++) {
        val |= T(v[i]) << i;
    }
    return val;
}
}

TEST(bits_test, erase_bit) {
    uint64_t x = 0x7BD152B330F0A777; uint64_t y = ~x;
    for (int i = 0; i < 64; i++) ASSERT_EQ(erase_slow(uint64_t(x), i), erase_bit(uint64_t(x), i)) << i;
    for (int i = 0; i < 64; i++) ASSERT_EQ(erase_slow(uint64_t(y), i), erase_bit(uint64_t(y), i)) << i;
    for (int i = 0; i < 32; i++) ASSERT_EQ(erase_slow(uint32_t(x), i), erase_bit(uint32_t(x), i)) << i;
    for (int i = 0; i < 32; i++) ASSERT_EQ(erase_slow(uint32_t(y), i), erase_bit(uint32_t(y), i)) << i;
    for (int i = 0; i < 16; i++) ASSERT_EQ(erase_slow(uint16_t(x), i), erase_bit(uint16_t(x), i)) << i;
    for (int i = 0; i < 16; i++) ASSERT_EQ(erase_slow(uint16_t(y), i), erase_bit(uint16_t(y), i)) << i;
    for (int i = 0; i < 8; i++) ASSERT_EQ(erase_slow(uint8_t(x), i), erase_bit(uint8_t(x), i)) << i;
    for (int i = 0; i < 8; i++) ASSERT_EQ(erase_slow(uint8_t(y), i), erase_bit(uint8_t(y), i)) << i;
}

TEST(bits_test, mix_bits) {
    uint64_t x = 0x7BD152B330F0A777;
    uint64_t y = 0x1234567890ABCDEF;
    ASSERT_EQ(0x7BD152B330F0A777, mix_bits<uint64_t>(x, y, 0x0000000000000000));
    ASSERT_EQ(0x1234567890ABCDEF, mix_bits<uint64_t>(x, y, 0xFFFFFFFFFFFFFFFF));
    ASSERT_EQ(0x1B31527390A0C7E7, mix_bits<uint64_t>(x, y, 0xF0F0F0F0F0F0F0F0));
    for (int i = 0; i < 64; i++) {
        ASSERT_EQ((clear_bit(x, i) | (get_bit(y, i) << i)), mix_bits(x, y, make_bit<uint64_t>(i)));
    }
}

TEST(bits_test, log2) {
    // for x = 0 the two implementations differ
    vector<int> v; for (int x = 0; x <= 32; x++) v.push_back(ilog2_64(x));
    ASSERT_EQ((vector<int>{-1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5}), v);

    ASSERT_EQ(0, ilog2(uint8_t(0)));
    ASSERT_EQ(0, ilog2(uint16_t(0)));
    ASSERT_EQ(0, ilog2(uint32_t(0)));
    ASSERT_EQ(0, ilog2(uint64_t(0)));

    ASSERT_EQ(6, ilog2(uint8_t(0x7F)));
    ASSERT_EQ(7, ilog2(uint8_t(0x80)));
    ASSERT_EQ(7, ilog2(uint8_t(0xFF)));

    ASSERT_EQ(14, ilog2(uint16_t(0x7FFF)));
    ASSERT_EQ(15, ilog2(uint16_t(0x8000)));
    ASSERT_EQ(15, ilog2(uint16_t(0xFFFF)));

    ASSERT_EQ(30, ilog2(uint32_t(0x7FFFFFFF)));
    ASSERT_EQ(31, ilog2(uint32_t(0x80000000)));
    ASSERT_EQ(31, ilog2(uint32_t(0xFFFFFFFF)));

    ASSERT_EQ(62, ilog2(uint64_t(0x7FFFFFFFFFFFFFFF)));
    ASSERT_EQ(63, ilog2(uint64_t(0x8000000000000000)));
    ASSERT_EQ(63, ilog2(uint64_t(0xFFFFFFFFFFFFFFFF)));

    for (int x = 1; x < (1 << 8); x++) ASSERT_EQ(ilog2_64(uint8_t(x)), ilog2(uint8_t(x))) << "8bit x: " << x;
    for (int x = 1; x < (1 << 16); x++) ASSERT_EQ(ilog2_64(uint16_t(x)), ilog2(uint16_t(x))) << "16bit x: " << x;
    for (int x = 1; x < (1 << 20); x++) ASSERT_EQ(ilog2_64(uint32_t(x)), ilog2(uint32_t(x))) << "32bit x: " << x;
    for (int x = 1; x < (1 << 20); x++) ASSERT_EQ(ilog2_64(uint64_t(x) << 20), ilog2(uint64_t(x) << 20)) << "64bit x: " << x;
}

TEST(bits_test, bit_cnt1) {
    vector<int> v; for (int x = 0; x <= 32; x++) v.push_back(bit_cnt1_64(x));
    ASSERT_EQ((vector<int>{0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 1}), v);

    ASSERT_EQ(7, bit_cnt1(uint8_t(0x7F)));
    ASSERT_EQ(1, bit_cnt1(uint8_t(0x80)));
    ASSERT_EQ(8, bit_cnt1(uint8_t(0xFF)));

    ASSERT_EQ(15, bit_cnt1(uint16_t(0x7FFF)));
    ASSERT_EQ( 1, bit_cnt1(uint16_t(0x8000)));
    ASSERT_EQ(16, bit_cnt1(uint16_t(0xFFFF)));

    ASSERT_EQ(31, bit_cnt1(uint32_t(0x7FFFFFFF)));
    ASSERT_EQ( 1, bit_cnt1(uint32_t(0x80000000)));
    ASSERT_EQ(32, bit_cnt1(uint32_t(0xFFFFFFFF)));

    ASSERT_EQ(63, bit_cnt1(uint64_t(0x7FFFFFFFFFFFFFFF)));
    ASSERT_EQ( 1, bit_cnt1(uint64_t(0x8000000000000000)));
    ASSERT_EQ(64, bit_cnt1(uint64_t(0xFFFFFFFFFFFFFFFF)));

    for (int x = 0; x < (1 << 8); x++) ASSERT_EQ(bit_cnt1_64(uint8_t(x)), bit_cnt1(uint8_t(x))) << "8bit x: " << x;
    for (int x = 0; x < (1 << 16); x++) ASSERT_EQ(bit_cnt1_64(uint16_t(x)), bit_cnt1(uint16_t(x))) << "16bit x: " << x;
    for (int x = 0; x < (1 << 20); x++) ASSERT_EQ(bit_cnt1_64(uint32_t(x)), bit_cnt1(uint32_t(x))) << "32bit x: " << x;
    for (int x = 0; x < (1 << 20); x++) ASSERT_EQ(bit_cnt1_64(uint64_t(x) << 20), bit_cnt1(uint64_t(x) << 20)) << "64bit x: " << x;
}

uint64_t shl(uint64_t x, int shift) { return x << shift; }

TEST(bits_test, bit_reverse) {
    vector<uint64_t> v; for (int x = 0; x <= 7; x++) v.push_back(bit_rev_64(x));
    ASSERT_EQ((vector<uint64_t>{0, shl(1, 63), shl(1, 62), shl(3, 62), shl(1, 61), shl(5, 61), shl(3, 61), shl(7, 61) }), v);

    ASSERT_EQ(uint8_t(0x95), bit_reverse(uint8_t(0xA9)));
    ASSERT_EQ(uint16_t(0x2395), bit_reverse(uint16_t(0xA9C4)));
    ASSERT_EQ(uint32_t(0x16402395), bit_reverse(uint32_t(0xA9C40268)));
    ASSERT_EQ(uint64_t(0xF7BDEAC816402395), bit_reverse(uint64_t(0xA9C402681357BDEF)));

    for (int x = 0; x < (1 << 8); x++) ASSERT_EQ(bit_rev_64(uint8_t(x)) >> 56, bit_reverse(uint8_t(x))) << "8bit x: " << x;
    for (int x = 0; x < (1 << 16); x++) ASSERT_EQ(bit_rev_64(uint16_t(x)) >> 48, bit_reverse(uint16_t(x))) << "16bit x: " << x;
    for (int x = 0; x < (1 << 20); x++) ASSERT_EQ(bit_rev_64(uint32_t(x)) >> 32, bit_reverse(uint32_t(x))) << "32bit x: " << x;
    for (int x = 0; x < (1 << 20); x++) ASSERT_EQ(bit_rev_64(uint64_t(x) << 20) >> 0, bit_reverse(uint64_t(x) << 20)) << "64bit x: " << x;
}

TEST(bits_test, or_down) {
    vector<uint64_t> v; for (int x = 0; x <= 32; x++) v.push_back(or_down_64(x));
    ASSERT_EQ((vector<uint64_t>{0, 1, 3, 3, 7, 7, 7, 7, 15, 15, 15, 15, 15, 15, 15, 15, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 63}), v);

    ASSERT_EQ(uint8_t(0x3F), or_down(uint8_t(0x39)));
    ASSERT_EQ(uint16_t(0x3FFF), or_down(uint16_t(0x29C4)));
    ASSERT_EQ(uint32_t(0xFFFFFFFF), or_down(uint32_t(0x89C40268)));
    ASSERT_EQ(uint64_t(0x7FFFFFFFFFFFFFFF), or_down(uint64_t(0x59C402681357BDEF)));

    for (int x = 0; x < (1 << 8); x++) ASSERT_EQ(or_down_64(uint8_t(x)), or_down(uint8_t(x))) << "8bit x: " << x;
    for (int x = 0; x < (1 << 16); x++) ASSERT_EQ(or_down_64(uint16_t(x)), or_down(uint16_t(x))) << "16bit x: " << x;
    for (int x = 0; x < (1 << 20); x++) ASSERT_EQ(or_down_64(uint32_t(x)), or_down(uint32_t(x))) << "32bit x: " << x;
    for (int x = 0; x < (1 << 20); x++) ASSERT_EQ(or_down_64(uint64_t(x) << 20), or_down(uint64_t(x) << 20)) << "64bit x: " << x;
}

TEST(bits_test, xor_down) {
    vector<int> v; for (int x = 0; x <= 32; x++) v.push_back(xor_down(uint32_t(x)));
    ASSERT_EQ((vector<int>{0, 1, 3, 2, 7, 6, 4, 5, 15, 14, 12, 13, 8, 9, 11, 10, 31, 30, 28, 29, 24, 25, 27, 26, 16, 17, 19, 18, 23, 22, 20, 21, 63}), v);

    ASSERT_EQ(uint8_t(0x2E), xor_down(uint8_t(0x39)));
    ASSERT_EQ(uint16_t(0x3178), xor_down(uint16_t(0x29C4)));
    ASSERT_EQ(uint32_t(0xF17803B0), xor_down(uint32_t(0x89C40268)));
    ASSERT_EQ(uint64_t(0x6E87FC4FE265294A), xor_down(uint64_t(0x59C402681357BDEF)));
}

TEST(bits_test, neg) {
    for (int x = 0; x < (1 << 8); x++) ASSERT_EQ(uint8_t(-x), neg(uint8_t(x))) << "8bit x: " << x;
    for (int x = 0; x < (1 << 16); x++) ASSERT_EQ(uint16_t(-x), neg(uint16_t(x))) << "16bit x: " << x;
    for (int x = 0; x < (1 << 20); x++) ASSERT_EQ(uint32_t(-x), neg(uint32_t(x))) << "32bit x: " << x;
    for (int x = 0; x < (1 << 20); x++) ASSERT_EQ(uint64_t(-(int64_t(x) << 20)), neg(uint64_t(x) << 20)) << "64bit x: " << x;
    for (int x = 0; x < (1 << 8); x++) ASSERT_EQ(uint8_t(x), neg(uint8_t(-x))) << "8bit x: " << x;
    for (int x = 0; x < (1 << 16); x++) ASSERT_EQ(uint16_t(x), neg(uint16_t(-x))) << "16bit x: " << x;
    for (int x = 0; x < (1 << 20); x++) ASSERT_EQ(uint32_t(x), neg(uint32_t(-x))) << "32bit x: " << x;
    for (int x = 0; x < (1 << 20); x++) ASSERT_EQ(uint64_t(int64_t(x) << 20), neg(uint64_t(-x) << 20)) << "64bit x: " << x;
}

TEST(bits_test, bin_to_gray_to_bin) {
    vector<int> vg; for (int x = 0; x <= 32; x++) vg.push_back(bin_to_gray(uint32_t(x)));
    ASSERT_EQ((vector<int>{0, 1, 3, 2, 6, 7, 5, 4, 12, 13, 15, 14, 10, 11, 9, 8, 24, 25, 27, 26, 30, 31, 29, 28, 20, 21, 23, 22, 18, 19, 17, 16, 48}), vg);

    ASSERT_EQ(uint8_t(0x39), bin_to_gray(uint8_t(0x2E)));
    ASSERT_EQ(uint16_t(0x29C4), bin_to_gray(uint16_t(0x3178)));
    ASSERT_EQ(uint32_t(0x89C40268), bin_to_gray(uint32_t(0xF17803B0)));
    ASSERT_EQ(uint64_t(0x59C402681357BDEF), bin_to_gray(uint64_t(0x6E87FC4FE265294A)));

    vector<int> vb; for (int x = 0; x <= 32; x++) vb.push_back(gray_to_bin(uint32_t(x)));
    ASSERT_EQ((vector<int>{0, 1, 3, 2, 7, 6, 4, 5, 15, 14, 12, 13, 8, 9, 11, 10, 31, 30, 28, 29, 24, 25, 27, 26, 16, 17, 19, 18, 23, 22, 20, 21, 63}), vb);

    ASSERT_EQ(uint8_t(0x2E), gray_to_bin(uint8_t(0x39)));
    ASSERT_EQ(uint16_t(0x3178), gray_to_bin(uint16_t(0x29C4)));
    ASSERT_EQ(uint32_t(0xF17803B0), gray_to_bin(uint32_t(0x89C40268)));
    ASSERT_EQ(uint64_t(0x6E87FC4FE265294A), gray_to_bin(uint64_t(0x59C402681357BDEF)));

    for (int x = 0; x < (1 << 8); x++) ASSERT_EQ(uint8_t(x), gray_to_bin(bin_to_gray(uint8_t(x)))) << "8bit x: " << x;
    for (int x = 0; x < (1 << 16); x++) ASSERT_EQ(uint16_t(x), gray_to_bin(bin_to_gray(uint16_t(x)))) << "16bit x: " << x;
    for (int x = 0; x < (1 << 20); x++) ASSERT_EQ(uint32_t(x), gray_to_bin(bin_to_gray(uint32_t(x)))) << "32bit x: " << x;
    for (int x = 0; x < (1 << 20); x++) ASSERT_EQ(uint64_t(x) << 20, gray_to_bin(bin_to_gray(uint64_t(x) << 20))) << "64bit x: " << x;

    for (int x = 0; x < (1 << 8); x++) ASSERT_EQ(uint8_t(x), bin_to_gray(gray_to_bin(uint8_t(x)))) << "8bit x: " << x;
    for (int x = 0; x < (1 << 16); x++) ASSERT_EQ(uint16_t(x), bin_to_gray(gray_to_bin(uint16_t(x)))) << "16bit x: " << x;
    for (int x = 0; x < (1 << 20); x++) ASSERT_EQ(uint32_t(x), bin_to_gray(gray_to_bin(uint32_t(x)))) << "32bit x: " << x;
    for (int x = 0; x < (1 << 20); x++) ASSERT_EQ(uint64_t(x) << 20, bin_to_gray(gray_to_bin(uint64_t(x) << 20))) << "64bit x: " << x;
}

TEST(bits_test, hi_bit) {
    vector<uint32_t> v; for (int x = 0; x <= 32; x++) v.push_back(hi_bit(uint32_t(x)));
    ASSERT_EQ((vector<uint32_t>{0, 1, 2, 2, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 8, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 32 }), v);

    ASSERT_EQ(uint8_t(0x20), hi_bit(uint8_t(0x39)));
    ASSERT_EQ(uint16_t(0x2000), hi_bit(uint16_t(0x29C4)));
    ASSERT_EQ(uint32_t(0x80000000), hi_bit(uint32_t(0x89C40268)));
    ASSERT_EQ(uint64_t(0x4000000000000000), hi_bit(uint64_t(0x59C402681357BDEF)));
}

TEST(bits_test, lo_bit) {
    vector<uint32_t> v; for (int x = 0; x <= 32; x++) v.push_back(lo_bit(uint32_t(x)));
    ASSERT_EQ((vector<uint32_t>{0, 1, 2, 1, 4, 1, 2, 1, 8, 1, 2, 1, 4, 1, 2, 1, 16, 1, 2, 1, 4, 1, 2, 1, 8, 1, 2, 1, 4, 1, 2, 1, 32 }), v);

    ASSERT_EQ(uint8_t(0x1), lo_bit(uint8_t(0x39)));
    ASSERT_EQ(uint16_t(0x4), lo_bit(uint16_t(0x29C4)));
    ASSERT_EQ(uint32_t(0x8), lo_bit(uint32_t(0x89C40268)));
    ASSERT_EQ(uint64_t(0x1), lo_bit(uint64_t(0x59C402681357BDEF)));
}

TEST(bits_test, is_pow_2) {
    vector<int> e((1 << 20) + 1); e[0] = 1; for (int k = 0; k <= 20; k++) e[1 << k] = 1;
    for (int x = 0; x < (1 << 8); x++) ASSERT_EQ(e[x], (int)is_pow2(uint8_t(x))) << "8bit x: " << x;
    for (int x = 0; x < (1 << 16); x++) ASSERT_EQ(e[x], (int)is_pow2(uint16_t(x))) << "16bit x: " << x;
    for (int x = 0; x < (1 << 20); x++) ASSERT_EQ(e[x], (int)is_pow2(uint32_t(x))) << "32bit x: " << x;
    for (int x = 0; x < (1 << 20); x++) ASSERT_EQ(e[x], (int)is_pow2(uint64_t(x))) << "64bit x: " << x;

    for (int x = 0; x < (1 << 8); x++) ASSERT_EQ(1 - e[x], (int)is_not_pow2(uint8_t(x))) << "8bit x: " << x;
    for (int x = 0; x < (1 << 16); x++) ASSERT_EQ(1 - e[x], (int)is_not_pow2(uint16_t(x))) << "16bit x: " << x;
    for (int x = 0; x < (1 << 20); x++) ASSERT_EQ(1 - e[x], (int)is_not_pow2(uint32_t(x))) << "32bit x: " << x;
    for (int x = 0; x < (1 << 20); x++) ASSERT_EQ(1 - e[x], (int)is_not_pow2(uint64_t(x))) << "64bit x: " << x;

    for (int k = 0; k < 8; k++) ASSERT_EQ(true, is_pow2(uint8_t(1) << k)) << "8bit k: " << k;
    for (int k = 0; k < 16; k++) ASSERT_EQ(true, is_pow2(uint16_t(1) << k)) << "16bit k: " << k;
    for (int k = 0; k < 32; k++) ASSERT_EQ(true, is_pow2(uint32_t(1) << k)) << "32bit k: " << k;
    for (int k = 0; k < 64; k++) ASSERT_EQ(true, is_pow2(uint64_t(1) << k)) << "64bit k: " << k;
}

TEST(bits_test, next_pow2) {
    vector<uint32_t> v; for (int x = 0; x <= 32; x++) v.push_back(next_pow2(uint32_t(x)));
    ASSERT_EQ((vector<uint32_t>{1, 2, 4, 4, 8, 8, 8, 8, 16, 16, 16, 16, 16, 16, 16, 16, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 64 }), v);

    ASSERT_EQ(uint8_t(0x40), next_pow2(uint8_t(0x39)));
    ASSERT_EQ(uint16_t(0x4000), next_pow2(uint16_t(0x29C4)));
    ASSERT_EQ(uint32_t(0x20000000), next_pow2(uint32_t(0x19C40268)));
    ASSERT_EQ(uint64_t(0x8000000000000000), next_pow2(uint64_t(0x59C402681357BDEF)));
}


TEST(bits_test, lzc) {
    vector<int> v; for (int x = 0; x <= 32; x++) v.push_back(lzc(uint32_t(x)));
    ASSERT_EQ((vector<int>{32, 31, 30, 30, 29, 29, 29, 29, 28, 28, 28, 28, 28, 28, 28, 28, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 26 }), v);

    ASSERT_EQ(2, lzc(uint8_t(0x30)));
    ASSERT_EQ(2, lzc(uint16_t(0x2000)));
    ASSERT_EQ(0, lzc(uint32_t(0x80000000)));
    ASSERT_EQ(1, lzc(uint64_t(0x4000000000000000)));
}

TEST(bits_test, tzc) {
    vector<uint32_t> v; for (int x = 0; x <= 32; x++) v.push_back(tzc(uint32_t(x)));
    ASSERT_EQ((vector<uint32_t>{32, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5 }), v);

    ASSERT_EQ(4, tzc(uint8_t(0x30)));
    ASSERT_EQ(13, tzc(uint16_t(0x2000)));
    ASSERT_EQ(31, tzc(uint32_t(0x80000000)));
    ASSERT_EQ(62, tzc(uint64_t(0x4000000000000000)));
}

TEST(bits_test, sign_mag) {
    // positive
    for (int x = 0; x < (1 << 7); x++) ASSERT_EQ(uint8_t(x), sign_mag(uint8_t(x))) << "8bit x: " << x;
    for (int x = 0; x < (1 << 15); x++) ASSERT_EQ(uint16_t(x), sign_mag(uint16_t(x))) << "16bit x: " << x;
    for (int x = 0; x < (1 << 20); x++) ASSERT_EQ(uint32_t(x), sign_mag(uint32_t(x))) << "32bit x: " << x;
    for (int x = 0; x < (1 << 20); x++) ASSERT_EQ(uint64_t(x), sign_mag(uint64_t(x))) << "64bit x: " << x;
    // negative: Sign-Magnitude to Two's complement
    for (int x = 0; x < (1 << 7); x++) ASSERT_EQ(uint8_t(-x), sign_mag(uint8_t(x | uint8_t(0x80)))) << "8bit x: " << x;
    for (int x = 0; x < (1 << 15); x++) ASSERT_EQ(uint16_t(-x), sign_mag(uint16_t(x | uint16_t(0x8000)))) << "16bit x: " << x;
    for (int x = 0; x < (1 << 20); x++) ASSERT_EQ(uint32_t(-x), sign_mag(uint32_t(x | uint32_t(0x80000000)))) << "32bit x: " << x;
    for (int x = 0; x < (1 << 20); x++) ASSERT_EQ(uint64_t(-x), sign_mag(uint64_t(x | uint64_t(0x8000000000000000)))) << "64bit x: " << x;
    // negative: Two's complement to Sign-Magnitude
    for (int x = 1; x < (1 << 7); x++) ASSERT_EQ(uint8_t(x | uint8_t(0x80)), sign_mag(uint8_t(-x))) << "8bit x: " << x;
    for (int x = 1; x < (1 << 15); x++) ASSERT_EQ(uint16_t(x | uint16_t(0x8000)), sign_mag(uint16_t(-x))) << "16bit x: " << x;
    for (int x = 1; x < (1 << 20); x++) ASSERT_EQ(uint32_t(x | uint32_t(0x80000000)), sign_mag(uint32_t(-x))) << "32bit x: " << x;
    for (int x = 1; x < (1 << 20); x++) ASSERT_EQ(uint64_t(x | uint64_t(0x8000000000000000)), sign_mag(uint64_t(-x))) << "64bit x: " << x;
}

TEST(bits_test, next_combination) {
    for (int n = 0; n <= 10; n++) {
        vector<vector<uint32_t>> v(n + 1);
        for (uint32_t w = 0; w < make_bit<uint32_t>(n); w++) {
            v[bit_cnt1(w)].push_back(w);
        }
        for (int k = 0; k <= n; k++) {
            for (int i = 0; i < v[k].size() - 1; i++) {
                uint32_t x = v[k][i];
                ASSERT_EQ(true, next_combination(x, n));
                ASSERT_EQ(v[k][i + 1], x);
            }
            uint32_t x = v[k].back();
            ASSERT_EQ(false, next_combination(x, n));
            ASSERT_EQ(v[k].front(), x);
        }
    }
}

TEST(bits_test, perf) {
    return;
    int r = 0;
    auto T0 = clock();
    for (uint64_t x = 1; x < (uint64_t(1) << 33); x++) r += ilog2(x); // 33759 ms
    double dT0 = double(clock() - T0) / CLOCKS_PER_SEC;
    printf("%0.2lf s   %p\n", dT0, &r);

    auto T1 = clock();
    for (uint64_t x = 1; x < (uint64_t(1) << 33); x++) r += ilog2_64(x); // 57907 ms
    double dT1 = double(clock() - T1) / CLOCKS_PER_SEC;
    printf("%0.2lf s   %p\n", dT1, &r);
}
