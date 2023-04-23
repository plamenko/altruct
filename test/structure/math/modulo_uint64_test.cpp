#include "altruct/structure/math/modulo.h"
#include "structure_test_util.h"

#include "gtest/gtest.h"

#include <functional>

using namespace std;
using namespace altruct::math;
using namespace altruct::test_util;

// largest prime that fits uint64_t: 18446744073709551557 = 2^64 - 59
using mod = modulo<uint64_t, UINT64_C(18446744073709551557), modulo_storage::CONSTANT>;

TEST(modulo_uint64_test, standalone_functions_1000000000000000003) {
    const uint64_t O = UINT64_C(0);
    const uint64_t M = UINT64_C(1000000000000000003);
    EXPECT_EQ(O + 0u, modulo_normalize(INT64_C(-1000000000000000003), M));
    EXPECT_EQ(O + 0u, modulo_normalize(O + 0u, M));
    EXPECT_EQ(O + 12u, modulo_normalize(M + 12, M));
    EXPECT_EQ(M - 1u, modulo_add(M - 3u, O + 2u, M));
    EXPECT_EQ(M - 5u, modulo_add(M - 3u, M - 2u, M));
    EXPECT_EQ(O + 9u, modulo_add(O + 9u, M, M));
    EXPECT_EQ(M - 5u, modulo_sub(M - 3u, O + 2, M));
    EXPECT_EQ(M - 1u, modulo_sub(M - 3u, M - 2u, M));
    EXPECT_EQ(O + 14u, modulo_sub(O + 14u, M, M));
    EXPECT_EQ(O + 0u, modulo_neg(O + 0u, M));
    EXPECT_EQ(M - 2u, modulo_neg(O + 2u, M));
    EXPECT_EQ(O + 3u, modulo_neg(M - 3u, M));
    EXPECT_EQ(O + 15u, modulo_mul(O + 3u, O + 5u, M));
    EXPECT_EQ(M - 6u, modulo_mul(O + 3u, M - 2u, M));
    EXPECT_EQ(M - 6u, modulo_mul(M - 3u, O + 2u, M));
    EXPECT_EQ(O + 18u, modulo_mul(M - 3u, M - 6u, M));
    EXPECT_EQ(O + 1u, modulo_inv(O + 1u, M));
    EXPECT_EQ(M - 1u, modulo_inv(M - 1u, M));
    EXPECT_EQ(UINT64_C(500000000000000002), modulo_inv(O + 2u, M));
    EXPECT_EQ(O + 2u, modulo_inv(UINT64_C(500000000000000002), M));
    EXPECT_EQ(UINT64_C(666666666666666669), modulo_inv(O + 3u, M));
    EXPECT_EQ(O + 3u, modulo_inv(UINT64_C(666666666666666669), M));
    EXPECT_EQ(O + 0u, modulo_div(O + 0u, O + 7u, M));
    EXPECT_EQ(O + 7u, modulo_div(O + 7u, O + 1u, M));
    EXPECT_EQ(UINT64_C(142857142857142858), modulo_div(O + 3u, O + 7u, M));
    EXPECT_EQ(O + 7u, modulo_div(O + 3u, UINT64_C(142857142857142858), M));
}

TEST(modulo_uint64_test, constructor) {
    const uint64_t M = UINT64_C(18446744073709551557);
    // default
    const mod m1;
    EXPECT_EQ(UINT64_C(0), m1.v);
    EXPECT_EQ(M, m1.M());
    // value only
    const mod m2(UINT64_C(10));
    EXPECT_EQ(UINT64_C(10), m2.v);
    EXPECT_EQ(M, m2.M());
    // value + modulus, modulus ignored
    const mod m3(UINT64_C(13), UINT64_C(12345));
    EXPECT_EQ(UINT64_C(13), m3.v);
    EXPECT_EQ(M, m3.M());

    // from different integral type: uint32_t
    const mod mu32_0(UINT32_C(0));
    EXPECT_EQ(UINT64_C(0), mu32_0.v);
    EXPECT_EQ(M, mu32_0.M());
    const mod mu32_1(UINT32_C(10));
    EXPECT_EQ(UINT64_C(10), mu32_1.v);
    EXPECT_EQ(M, mu32_1.M());
    const mod mu32_2(UINT32_C(4294967290));
    EXPECT_EQ(UINT64_C(4294967290), mu32_2.v);
    EXPECT_EQ(M, mu32_2.M());
    const mod mu32_3(UINT32_C(4294967292));
    EXPECT_EQ(UINT64_C(4294967292), mu32_3.v);
    EXPECT_EQ(M, mu32_3.M());

    // from different integral type: int32_t
    const mod mi32_0(INT32_C(0));
    EXPECT_EQ(UINT64_C(0), mi32_0.v);
    EXPECT_EQ(M, mi32_0.M());
    const mod mi32_1(INT32_C(20));
    EXPECT_EQ(UINT64_C(20), mi32_1.v);
    EXPECT_EQ(M, mi32_1.M());
    const mod mi32_2(INT32_C(-2));
    EXPECT_EQ(UINT64_C(18446744073709551555), mi32_2.v);
    EXPECT_EQ(M, mi32_2.M());
    const mod mi32_3(INT32_C(-102));
    EXPECT_EQ(UINT64_C(18446744073709551455), mi32_3.v);
    EXPECT_EQ(M, mi32_3.M());

    // from same integral type: uint64_t
    const mod mu64_0(UINT64_C(0));
    EXPECT_EQ(UINT64_C(0), mu64_0.v);
    EXPECT_EQ(M, mu64_0.M());
    const mod mu64_1(UINT64_C(40));
    EXPECT_EQ(UINT64_C(40), mu64_1.v);
    EXPECT_EQ(M, mu64_1.M());
    const mod mu64_2(UINT64_C(18446744073709551553)); // -4
    EXPECT_EQ(UINT64_C(18446744073709551553), mu64_2.v);
    EXPECT_EQ(M, mu64_2.M());
    const mod mu64_3(UINT64_C(18446744073709551453)); // -104
    EXPECT_EQ(UINT64_C(18446744073709551453), mu64_3.v);
    EXPECT_EQ(M, mu64_3.M());
    const mod mu64_4(UINT64_C(18446744073709551561)); // 4
    EXPECT_EQ(UINT64_C(4), mu64_4.v);
    EXPECT_EQ(M, mu64_4.M());
    const mod mu64_5(UINT64_C(1000000000000));
    EXPECT_EQ(UINT64_C(1000000000000), mu64_5.v);
    EXPECT_EQ(M, mu64_5.M());

    // from different integral type: int64_t
    const mod mi64_0(INT64_C(0));
    EXPECT_EQ(UINT64_C(0), mi64_0.v);
    EXPECT_EQ(M, mi64_0.M());
    const mod mi64_1(INT64_C(50));
    EXPECT_EQ(UINT64_C(50), mi64_1.v);
    EXPECT_EQ(M, mi64_1.M());
    const mod mi64_2(INT64_C(-5));
    EXPECT_EQ(UINT64_C(18446744073709551552), mi64_2.v);
    EXPECT_EQ(M, mi64_2.M());
    const mod mi64_3(INT64_C(-105));
    EXPECT_EQ(UINT64_C(18446744073709551452), mi64_3.v);
    EXPECT_EQ(M, mi64_3.M());
    const mod mi64_4(INT64_C(4294967296));
    EXPECT_EQ(UINT64_C(4294967296), mi64_4.v);
    EXPECT_EQ(M, mi64_4.M());
    const mod mi64_5(INT64_C(1000000000000));
    EXPECT_EQ(UINT64_C(1000000000000), mi64_5.v);
    EXPECT_EQ(M, mi64_5.M());
    const mod mi64_6(INT64_C(-1000000000000));
    EXPECT_EQ(UINT64_C(18446743073709551557), mi64_6.v);
    EXPECT_EQ(M, mi64_6.M());

    // from different integral type: int64_t
    // value + modulus, modulus ignored
    const mod mi64_7(INT64_C(-1000000000000), UINT64_C(12345));
    EXPECT_EQ(UINT64_C(18446743073709551557), mi64_6.v);
    EXPECT_EQ(M, mi64_6.M());

    // copy constructor
    const mod mu64_c(mu64_1);
    EXPECT_EQ(UINT64_C(40), mu64_c.v);
    EXPECT_EQ(M, mu64_c.M());
    // move constructor
    const mod mu64_m(std::move(mu64_2));
    EXPECT_EQ(UINT64_C(18446744073709551553), mu64_m.v);
    EXPECT_EQ(M, mu64_m.M());
    // assignment
    mod mu64_a; mu64_a = mu64_1;
    EXPECT_EQ(UINT64_C(40), mu64_a.v);
    EXPECT_EQ(M, mu64_a.M());
    // move assignment
    mu64_a = std::move(mu64_3);
    EXPECT_EQ(UINT64_C(18446744073709551453), mu64_a.v);
    EXPECT_EQ(M, mu64_a.M());
}

TEST(modulo_uint64_test, operators_comparison) {
    const mod m1 = 10;
    const mod m2 = 20;
    ASSERT_COMPARISON_OPERATORS(0, m1, m1);
    ASSERT_COMPARISON_OPERATORS(0, m2, m2);
    ASSERT_COMPARISON_OPERATORS(-1, m1, m2);
    ASSERT_COMPARISON_OPERATORS(+1, m2, m1);
}

TEST(modulo_uint64_test, operators_arithmetic) {
    const uint64_t M = UINT64_C(18446744073709551557);
    const mod m1 = -7;
    const mod m2 = 9;
    const mod m3 = -21;
    EXPECT_EQ(mod(-7), m1);
    EXPECT_EQ(mod(9), m2);
    EXPECT_EQ(mod(-21), m3);
    EXPECT_EQ(mod(2), m1 + m2);
    EXPECT_EQ(mod(-16), m1 - m2);
    EXPECT_EQ(mod(7), -m1);
    EXPECT_EQ(mod(-63), m1 * m2);
    EXPECT_EQ(mod(UINT64_C(16397105843297379161)), m1 / m2);
    EXPECT_EQ(mod(4), m1 % m2);
    EXPECT_EQ(mod(2), m2 + m1);
    EXPECT_EQ(mod(16), m2 - m1);
    EXPECT_EQ(mod(-9), -m2);
    EXPECT_EQ(mod(-63), m2 * m1);
    EXPECT_EQ(mod(UINT64_C(13176245766935393968)), m2 / m1);
    EXPECT_EQ(mod(9), m2 % m1);
    EXPECT_EQ(mod(3), m3 / m1);
    EXPECT_EQ(mod(UINT64_C(6148914691236517186)), m1 / m3);
}

TEST(modulo_uint64_test, operators_inplace) {
    const uint64_t M = UINT64_C(18446744073709551557);
    const mod m1 = -7;
    const mod m2 = 9;
    const mod m3 = -21;
    mod mr;
    mr = m1; mr += m2;
    EXPECT_EQ(mod(2), mr);
    mr = m1; mr -= m2;
    EXPECT_EQ(mod(-16), mr);
    mr = m1; mr *= m2;
    EXPECT_EQ(mod(-63), mr);
    mr = m1; mr /= m2;
    EXPECT_EQ(mod(UINT64_C(16397105843297379161)), mr);
    mr = m1; mr %= m2;
    EXPECT_EQ(mod(4), mr);
    mr = m2; mr += m1;
    EXPECT_EQ(mod(2), mr);
    mr = m2; mr -= m1;
    EXPECT_EQ(mod(16), mr);
    mr = m2; mr *= m1;
    EXPECT_EQ(mod(-63), mr);
    mr = m2; mr /= m1;
    EXPECT_EQ(mod(UINT64_C(13176245766935393968)), mr);
    mr = m2; mr %= m1;
    EXPECT_EQ(mod(9), mr);
    mr = m3; mr /= m1;
    EXPECT_EQ(mod(3), m3 / m1);
    mr = m1; mr /= m3;
    EXPECT_EQ(mod(UINT64_C(6148914691236517186)), m1 / m3);
}

TEST(modulo_uint64_test, operators_inplace_self) {
    const uint64_t M = UINT64_C(18446744073709551557);
    const mod m1 = -7;
    mod mr;
    mr = m1; mr += mr;
    EXPECT_EQ(mod(-14), mr);
    mr = m1; mr -= mr;
    EXPECT_EQ(mod(0), mr);
    mr = m1; mr *= mr;
    EXPECT_EQ(mod(49), mr);
    mr = m1; mr /= mr;
    EXPECT_EQ(mod(1), mr);
    mr = m1; mr %= mr;
    EXPECT_EQ(mod(0), mr);
}

TEST(modulo_uint64_test, casts) {
    const uint64_t M = UINT64_C(18446744073709551557);
    const mod m1 = -7;
    const mod e0 = zeroOf(m1);
    const mod e1 = identityOf(m1);
    EXPECT_EQ(UINT64_C(0), e0.v);
    EXPECT_EQ(M, e0.M());
    EXPECT_EQ(UINT64_C(1), e1.v);
    EXPECT_EQ(M, e1.M());
    const mod m3 = castOf<mod>(INT64_C(-1000000000000));
    EXPECT_EQ(UINT64_C(18446743073709551557), m3.v);
    EXPECT_EQ(M, m3.M());
    const mod m5 = castOf(m1, -5);
    EXPECT_EQ(UINT64_C(18446744073709551552), m5.v);
    EXPECT_EQ(M, m5.M());
    const mod m6 = castOf(m1, m5);
    EXPECT_EQ(UINT64_C(18446744073709551552), m6.v);
    EXPECT_EQ(M, m6.M());
    const mod m7 = castOf<mod>(m5);
    EXPECT_EQ(UINT64_C(18446744073709551552), m7.v);
    EXPECT_EQ(M, m7.M());
    EXPECT_EQ(UINT64_C(4), modT(UINT64_C(18446744073709551561), M));
    const mod m8 = powT(m1, 100);
    EXPECT_EQ(UINT64_C(6708427641812857077), m8.v);
    EXPECT_EQ(M, m8.M());
}
