#include "altruct/structure/math/modulo.h"
#include "structure_test_util.h"

#include "gtest/gtest.h"

#include <functional>
#include <utility>

using namespace std;
using namespace altruct::math;
using namespace altruct::test_util;

typedef moduloX<int32_t> modx;
typedef std::pair<int32_t, int32_t> pii;

namespace {
pii to_pair(const modx& m) { return { m.v, m.M() }; }
} // namespace

TEST(modulox_int32_test, constructor) {
    const int32_t M = INT32_C(2147483629);
    // default
    const modx m1;
    EXPECT_EQ(INT32_C(0), m1.v);
    EXPECT_EQ(INT32_C(1), m1.M());
    // value only
    const modx m2(INT32_C(10));
    EXPECT_EQ(INT32_C(10), m2.v);
    EXPECT_EQ(INT32_C(1), m2.M());
    // value + modulus
    const modx m3(INT32_C(13), M);
    EXPECT_EQ(INT32_C(13), m3.v);
    EXPECT_EQ(M, m3.M());

    // from different integral type: uint32_t
    const modx mu32_0(UINT32_C(0), M);
    EXPECT_EQ(INT32_C(0), mu32_0.v);
    EXPECT_EQ(M, mu32_0.M());
    const modx mu32_1(UINT32_C(10), M);
    EXPECT_EQ(INT32_C(10), mu32_1.v);
    EXPECT_EQ(M, mu32_1.M());
    const modx mu32_2(UINT32_C(2147483628), M); // -1
    EXPECT_EQ(INT32_C(2147483628), mu32_2.v);
    EXPECT_EQ(M, mu32_2.M());
    const modx mu32_3(UINT32_C(2147483630), M); // +1
    EXPECT_EQ(INT32_C(1), mu32_3.v);
    EXPECT_EQ(M, mu32_3.M());

    // from same integral type: int32_t
    const modx mi32_0(INT32_C(0), M);
    EXPECT_EQ(INT32_C(0), mi32_0.v);
    EXPECT_EQ(M, mi32_0.M());
    const modx mi32_1(INT32_C(20), M);
    EXPECT_EQ(INT32_C(20), mi32_1.v);
    EXPECT_EQ(M, mi32_1.M());
    const modx mi32_2(INT32_C(-2), M);
    EXPECT_EQ(INT32_C(2147483627), mi32_2.v);
    EXPECT_EQ(M, mi32_2.M());
    const modx mi32_3(INT32_C(-102), M);
    EXPECT_EQ(INT32_C(2147483527), mi32_3.v);
    EXPECT_EQ(M, mi32_3.M());

    // from different integral type: uint64_t
    const modx mu64_0(UINT64_C(0), M);
    EXPECT_EQ(INT32_C(0), mu64_0.v);
    EXPECT_EQ(M, mu64_0.M());
    const modx mu64_1(UINT64_C(40), M);
    EXPECT_EQ(INT32_C(40), mu64_1.v);
    EXPECT_EQ(M, mu64_1.M());
    const modx mu64_2(UINT64_C(4294967254), M); // -4
    EXPECT_EQ(INT32_C(2147483625), mu64_2.v);
    EXPECT_EQ(M, mu64_2.M());
    const modx mu64_3(UINT64_C(4294967154), M); // -104
    EXPECT_EQ(INT32_C(2147483525), mu64_3.v);
    EXPECT_EQ(M, mu64_3.M());
    const modx mu64_4(UINT64_C(4294967262), M); // 4
    EXPECT_EQ(INT32_C(4), mu64_4.v);
    EXPECT_EQ(M, mu64_4.M());
    const modx mu64_5(UINT64_C(1000000000000), M);
    EXPECT_EQ(INT32_C(1420112515), mu64_5.v);
    EXPECT_EQ(M, mu64_5.M());

    // from different integral type: int64_t
    const modx mi64_0(INT64_C(0), M);
    EXPECT_EQ(INT32_C(0), mi64_0.v);
    EXPECT_EQ(M, mi64_0.M());
    const modx mi64_1(INT64_C(50), M);
    EXPECT_EQ(INT32_C(50), mi64_1.v);
    EXPECT_EQ(M, mi64_1.M());
    const modx mi64_2(INT64_C(-5), M);
    EXPECT_EQ(INT32_C(2147483624), mi64_2.v);
    EXPECT_EQ(M, mi64_2.M());
    const modx mi64_3(INT64_C(-105), M);
    EXPECT_EQ(INT32_C(2147483524), mi64_3.v);
    EXPECT_EQ(M, mi64_3.M());
    const modx mi64_4(INT64_C(4294967296), M);
    EXPECT_EQ(INT32_C(38), mi64_4.v);
    EXPECT_EQ(M, mi64_4.M());
    const modx mi64_5(INT64_C(1000000000000), M);
    EXPECT_EQ(INT32_C(1420112515), mi64_5.v);
    EXPECT_EQ(M, mi64_5.M());
    const modx mi64_6(INT64_C(-1000000000000), M);
    EXPECT_EQ(INT32_C(727371114), mi64_6.v);
    EXPECT_EQ(M, mi64_6.M());

    // from different integral type: uint64_t
    // value + modulus
    const modx mu64_7(UINT64_C(1000000000000), M);
    EXPECT_EQ(INT32_C(1420112515), mu64_7.v);
    EXPECT_EQ(M, mu64_7.M());

    // copy constructor
    const modx mi32_c(mi32_1);
    EXPECT_EQ(INT32_C(20), mi32_c.v);
    EXPECT_EQ(M, mi32_c.M());
    // move constructor
    const modx mi32_m(std::move(mi32_2));
    EXPECT_EQ(INT32_C(2147483627), mi32_m.v);
    EXPECT_EQ(M, mi32_m.M());
    // assignment
    modx mi32_a; mi32_a = mi32_1;
    EXPECT_EQ(INT32_C(20), mi32_a.v);
    EXPECT_EQ(M, mi32_a.M());
    // move assignment
    mi32_a = std::move(mi32_3);
    EXPECT_EQ(INT32_C(2147483527), mi32_a.v);
    EXPECT_EQ(M, mi32_a.M());
}

TEST(modulox_int32_test, operators_comparison) {
    const int32_t M = INT32_C(2147483629);
    const modx m1(10, M);
    const modx m2(20, M);
    ASSERT_COMPARISON_OPERATORS(0, m1, m1);
    ASSERT_COMPARISON_OPERATORS(0, m2, m2);
    ASSERT_COMPARISON_OPERATORS(-1, m1, m2);
    ASSERT_COMPARISON_OPERATORS(+1, m2, m1);
}

TEST(modulox_int32_test, operators_arithmetic) {
    const int32_t M = INT32_C(2147483629);
    const modx m1(-7, M);
    const modx m2(9, M);
    const modx m3(-21, M);
    EXPECT_EQ(pii(M - 7, M), to_pair(m1));
    EXPECT_EQ(pii(9, M), to_pair(m2));
    EXPECT_EQ(pii(M - 21, M), to_pair(m3));
    EXPECT_EQ(pii(2, M), to_pair(m1 + m2));
    EXPECT_EQ(pii(M - 16, M), to_pair(m1 - m2));
    EXPECT_EQ(pii(7, M), to_pair(-m1));
    EXPECT_EQ(pii(M - 63, M), to_pair(m1 * m2));
    EXPECT_EQ(pii(INT32_C(1670265044), M), to_pair(m1 / m2));
    EXPECT_EQ(pii(3, M), to_pair(m1 % m2));
    EXPECT_EQ(pii(2, M), to_pair(m2 + m1));
    EXPECT_EQ(pii(16, M), to_pair(m2 - m1));
    EXPECT_EQ(pii(M - 9, M), to_pair(-m2));
    EXPECT_EQ(pii(M - 63, M), to_pair(m2 * m1));
    EXPECT_EQ(pii(INT32_C(1227133501), M), to_pair(m2 / m1));
    EXPECT_EQ(pii(9, M), to_pair(m2 % m1));
    EXPECT_EQ(pii(3, M), to_pair(m3 / m1));
    EXPECT_EQ(pii(INT32_C(1431655753), M), to_pair(m1 / m3));
}

TEST(modulox_int32_test, operators_inplace) {
    const int32_t M = INT32_C(2147483629);
    const modx m1(-7, M);
    const modx m2(9, M);
    const modx m3(-21, M);
    modx mr;
    mr = m1; mr += m2;
    EXPECT_EQ(pii(2, M), to_pair(mr));
    mr = m1; mr -= m2;
    EXPECT_EQ(pii(M-16, M), to_pair(mr));
    mr = m1; mr *= m2;
    EXPECT_EQ(pii(M-63, M), to_pair(mr));
    mr = m1; mr /= m2;
    EXPECT_EQ(pii(INT32_C(1670265044), M), to_pair(mr));
    mr = m1; mr %= m2;
    EXPECT_EQ(pii(3, M), to_pair(mr));
    mr = m2; mr += m1;
    EXPECT_EQ(pii(2, M), to_pair(mr));
    mr = m2; mr -= m1;
    EXPECT_EQ(pii(16, M), to_pair(mr));
    mr = m2; mr *= m1;
    EXPECT_EQ(pii(M-63, M), to_pair(mr));
    mr = m2; mr /= m1;
    EXPECT_EQ(pii(INT32_C(1227133501), M), to_pair(mr));
    mr = m2; mr %= m1;
    EXPECT_EQ(pii(9, M), to_pair(mr));
    mr = m3; mr /= m1;
    EXPECT_EQ(pii(3, M), to_pair(m3 / m1));
    mr = m1; mr /= m3;
    EXPECT_EQ(pii(INT32_C(1431655753), M), to_pair(m1 / m3));
}

TEST(modulox_int32_test, operators_inplace_self) {
    const int32_t M = INT32_C(2147483629);
    const modx m1(-7, M);
    modx mr;
    mr = m1; mr += mr;
    EXPECT_EQ(pii(M-14,M), to_pair(mr));
    mr = m1; mr -= mr;
    EXPECT_EQ(pii(0, M), to_pair(mr));
    mr = m1; mr *= mr;
    EXPECT_EQ(pii(49, M), to_pair(mr));
    mr = m1; mr /= mr;
    EXPECT_EQ(pii(1, M), to_pair(mr));
    mr = m1; mr %= mr;
    EXPECT_EQ(pii(0, M), to_pair(mr));
}

TEST(modulox_int32_test, casts) {
    const int32_t M = INT32_C(2147483629);
    const modx m1(-7, M);
    const modx e0 = zeroOf(m1);
    const modx e1 = identityOf(m1);
    EXPECT_EQ(INT32_C(0), e0.v);
    EXPECT_EQ(M, e0.M());
    EXPECT_EQ(INT32_C(1), e1.v);
    EXPECT_EQ(M, e1.M());
    const modx m3 = castOf(m1, INT64_C(1000000000000));
    EXPECT_EQ(INT32_C(1420112515), m3.v);
    EXPECT_EQ(M, m3.M());
    const modx m5 = castOf<modx>(5);
    EXPECT_EQ(INT32_C(5), m5.v);
    EXPECT_EQ(INT32_C(1), m5.M());
    const modx m6 = castOf(m1, m3);
    EXPECT_EQ(INT32_C(1420112515), m6.v);
    EXPECT_EQ(M, m6.M());
    const modx m7 = castOf<modx>(m3);
    EXPECT_EQ(INT32_C(1420112515), m7.v);
    EXPECT_EQ(M, m7.M());
    EXPECT_EQ(INT32_C(4), modT(INT32_C(2147483633), M));
    const modx m8 = powT(m1, 100);
    EXPECT_EQ(INT32_C(681305249), m8.v);
    EXPECT_EQ(M, m8.M());
}

TEST(modulox_int32_test, division) {
    // 18 directly divisible by 6
    EXPECT_EQ(modx(3, 1000), modx(18, 1000) / modx(6, 1000));
    EXPECT_EQ(modx(18, 1000), modx(3, 1000) * modx(6, 1000));

    // 7 is invertible modulo 1000
    EXPECT_EQ(modx(430, 1000), modx(10, 1000) / modx(7, 1000));
    EXPECT_EQ(modx(10, 1000), modx(430, 1000) * modx(7, 1000));

    // 48 is not invertible modulo 1000,
    // but after dividing all three (56, 48 and 1000)
    // by their GCD 8,  48/8=6 is invertible modulo 1000/8=125
    EXPECT_EQ(modx(147, 1000), modx(56, 1000) / modx(48, 1000));
    EXPECT_EQ(modx(56, 1000), modx(147, 1000) * modx(48, 1000));

    // 48 is not invertible modulo 1000,
    // and even after dividing all three (28, 48 and 1000)
    // by their GCD 4,  48/4=12 is still not invertible modulo 1000/4=250
    EXPECT_EQ(modx(0, 1000), modx(28, 1000) / modx(48, 1000));

    EXPECT_EQ(modx(53, 100), modx(17, 100).inv());
}
