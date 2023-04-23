#include "altruct/structure/math/modulo.h"
#include "structure_test_util.h"

#include "gtest/gtest.h"

#include <functional>

using namespace std;
using namespace altruct::math;
using namespace altruct::test_util;

// largest prime that fits int64_t: 9223372036854775783 = 2^63 - 25
using mod = modulo<int64_t, INT64_C(9223372036854775783), modulo_storage::CONSTANT>;

TEST(modulo_int64_test, standalone_functions_1000000000000000003) {
    const int64_t O = INT64_C(0);
    const int64_t M = INT64_C(1000000000000000003);
    EXPECT_EQ(O + 0, modulo_normalize(INT64_C(-1000000000000000003), M));
    EXPECT_EQ(O + 0, modulo_normalize(O + 0, M));
    EXPECT_EQ(O + 12, modulo_normalize(M + 12, M));
    EXPECT_EQ(M - 1, modulo_add(M - 3, O + 2, M));
    EXPECT_EQ(M - 5, modulo_add(M - 3, M - 2, M));
    EXPECT_EQ(O + 9, modulo_add(O + 9, M, M));
    EXPECT_EQ(M - 5, modulo_sub(M - 3, O + 2, M));
    EXPECT_EQ(M - 1, modulo_sub(M - 3, M - 2, M));
    EXPECT_EQ(O + 14, modulo_sub(O + 14, M, M));
    EXPECT_EQ(O + 0, modulo_neg(O + 0, M));
    EXPECT_EQ(M - 2, modulo_neg(O + 2, M));
    EXPECT_EQ(O + 3, modulo_neg(M - 3, M));
    EXPECT_EQ(O + 15, modulo_mul(O + 3, O + 5, M));
    EXPECT_EQ(M - 6, modulo_mul(O + 3, M - 2, M));
    EXPECT_EQ(M - 6, modulo_mul(M - 3, O + 2, M));
    EXPECT_EQ(O + 18, modulo_mul(M - 3, M - 6, M));
    EXPECT_EQ(O + 1, modulo_inv(O + 1, M));
    EXPECT_EQ(M - 1, modulo_inv(M - 1, M));
    EXPECT_EQ(INT64_C(500000000000000002), modulo_inv(O + 2, M));
    EXPECT_EQ(O + 2, modulo_inv(INT64_C(500000000000000002), M));
    EXPECT_EQ(INT64_C(666666666666666669), modulo_inv(O + 3, M));
    EXPECT_EQ(O + 3, modulo_inv(INT64_C(666666666666666669), M));
    EXPECT_EQ(O + 0, modulo_div(O + 0, O + 7, M));
    EXPECT_EQ(O + 7, modulo_div(O + 7, O + 1, M));
    EXPECT_EQ(INT64_C(142857142857142858), modulo_div(O + 3, O + 7, M));
    EXPECT_EQ(O + 7, modulo_div(O + 3, INT64_C(142857142857142858), M));
}

TEST(modulo_int64_test, constructor) {
    const int64_t M = INT64_C(9223372036854775783);
    // default
    const mod m1;
    EXPECT_EQ(INT64_C(0), m1.v);
    EXPECT_EQ(M, m1.M());
    // value only
    const mod m2(INT64_C(10));
    EXPECT_EQ(INT64_C(10), m2.v);
    EXPECT_EQ(M, m2.M());
    // value + modulus, modulus ignored
    const mod m3(INT64_C(13), INT64_C(12345));
    EXPECT_EQ(INT64_C(13), m3.v);
    EXPECT_EQ(M, m3.M());

    // from different integral type: uint32_t
    const mod mu32_0(UINT32_C(0));
    EXPECT_EQ(INT64_C(0), mu32_0.v);
    EXPECT_EQ(M, mu32_0.M());
    const mod mu32_1(UINT32_C(10));
    EXPECT_EQ(INT64_C(10), mu32_1.v);
    EXPECT_EQ(M, mu32_1.M());
    const mod mu32_2(UINT32_C(4294967290));
    EXPECT_EQ(INT64_C(4294967290), mu32_2.v);
    EXPECT_EQ(M, mu32_2.M());
    const mod mu32_3(UINT32_C(4294967292));
    EXPECT_EQ(INT64_C(4294967292), mu32_3.v);
    EXPECT_EQ(M, mu32_3.M());

    // from same integral type: int32_t
    const mod mi32_0(INT32_C(0));
    EXPECT_EQ(INT64_C(0), mi32_0.v);
    EXPECT_EQ(M, mi32_0.M());
    const mod mi32_1(INT32_C(20));
    EXPECT_EQ(INT64_C(20), mi32_1.v);
    EXPECT_EQ(M, mi32_1.M());
    const mod mi32_2(INT32_C(-2));
    EXPECT_EQ(INT64_C(9223372036854775781), mi32_2.v);
    EXPECT_EQ(M, mi32_2.M());
    const mod mi32_3(INT32_C(-102));
    EXPECT_EQ(INT64_C(9223372036854775681), mi32_3.v);
    EXPECT_EQ(M, mi32_3.M());

    // from different integral type: uint64_t
    const mod mu64_0(UINT64_C(0));
    EXPECT_EQ(INT64_C(0), mu64_0.v);
    EXPECT_EQ(M, mu64_0.M());
    const mod mu64_1(UINT64_C(40));
    EXPECT_EQ(INT64_C(40), mu64_1.v);
    EXPECT_EQ(M, mu64_1.M());
    const mod mu64_2(UINT64_C(9223372036854775779)); // -4
    EXPECT_EQ(INT64_C(9223372036854775779), mu64_2.v);
    EXPECT_EQ(M, mu64_2.M());
    const mod mu64_3(UINT64_C(9223372036854775679)); // -104
    EXPECT_EQ(INT64_C(9223372036854775679), mu64_3.v);
    EXPECT_EQ(M, mu64_3.M());
    const mod mu64_4(UINT64_C(9223372036854775787)); // 4
    EXPECT_EQ(INT64_C(4), mu64_4.v);
    EXPECT_EQ(M, mu64_4.M());
    const mod mu64_5(UINT64_C(1000000000000));
    EXPECT_EQ(INT64_C(1000000000000), mu64_5.v);
    EXPECT_EQ(M, mu64_5.M());

    // from same integral type: int64_t
    const mod mi64_0(INT64_C(0));
    EXPECT_EQ(INT64_C(0), mi64_0.v);
    EXPECT_EQ(M, mi64_0.M());
    const mod mi64_1(INT64_C(50));
    EXPECT_EQ(INT64_C(50), mi64_1.v);
    EXPECT_EQ(M, mi64_1.M());
    const mod mi64_2(INT64_C(-5));
    EXPECT_EQ(INT64_C(9223372036854775778), mi64_2.v);
    EXPECT_EQ(M, mi64_2.M());
    const mod mi64_3(INT64_C(-105));
    EXPECT_EQ(INT64_C(9223372036854775678), mi64_3.v);
    EXPECT_EQ(M, mi64_3.M());
    const mod mi64_4(INT64_C(4294967296));
    EXPECT_EQ(INT64_C(4294967296), mi64_4.v);
    EXPECT_EQ(M, mi64_4.M());
    const mod mi64_5(INT64_C(1000000000000));
    EXPECT_EQ(INT64_C(1000000000000), mi64_5.v);
    EXPECT_EQ(M, mi64_5.M());
    const mod mi64_6(INT64_C(-1000000000000));
    EXPECT_EQ(INT64_C(9223371036854775783), mi64_6.v);
    EXPECT_EQ(M, mi64_6.M());

    // from different integral type: uint64_t
    // value + modulus, modulus ignored
    const mod mi64_7(UINT64_C(10000000000000000000), INT64_C(12345));
    EXPECT_EQ(INT64_C(776627963145224217), mi64_7.v);
    EXPECT_EQ(M, mi64_7.M());

    // copy constructor
    const mod mi64_c(mi64_1);
    EXPECT_EQ(INT64_C(50), mi64_c.v);
    EXPECT_EQ(M, mi64_c.M());
    // move constructor
    const mod mi64_m(std::move(mi64_2));
    EXPECT_EQ(INT64_C(9223372036854775778), mi64_m.v);
    EXPECT_EQ(M, mi64_m.M());
    // assignment
    mod mi64_a; mi64_a = mi64_1;
    EXPECT_EQ(INT64_C(50), mi64_a.v);
    EXPECT_EQ(M, mi64_a.M());
    // move assignment
    mi64_a = std::move(mi64_3);
    EXPECT_EQ(INT64_C(9223372036854775678), mi64_a.v);
    EXPECT_EQ(M, mi64_a.M());
}

TEST(modulo_int64_test, operators_comparison) {
    const mod m1 = 10;
    const mod m2 = 20;
    ASSERT_COMPARISON_OPERATORS(0, m1, m1);
    ASSERT_COMPARISON_OPERATORS(0, m2, m2);
    ASSERT_COMPARISON_OPERATORS(-1, m1, m2);
    ASSERT_COMPARISON_OPERATORS(+1, m2, m1);
}

TEST(modulo_int64_test, operators_arithmetic) {
    const int64_t M = INT64_C(9223372036854775783);
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
    EXPECT_EQ(mod(INT64_C(7173733806442603386)), m1 / m2);
    EXPECT_EQ(mod(3), m1 % m2);
    EXPECT_EQ(mod(2), m2 + m1);
    EXPECT_EQ(mod(16), m2 - m1);
    EXPECT_EQ(mod(-9), -m2);
    EXPECT_EQ(mod(-63), m2 * m1);
    EXPECT_EQ(mod(INT64_C(5270498306774157589)), m2 / m1);
    EXPECT_EQ(mod(9), m2 % m1);
    EXPECT_EQ(mod(3), m3 / m1);
    EXPECT_EQ(mod(INT64_C(6148914691236517189)), m1 / m3);
}

TEST(modulo_int64_test, operators_inplace) {
    const int64_t M = INT64_C(9223372036854775783);
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
    EXPECT_EQ(mod(INT64_C(7173733806442603386)), mr);
    mr = m1; mr %= m2;
    EXPECT_EQ(mod(3), mr);
    mr = m2; mr += m1;
    EXPECT_EQ(mod(2), mr);
    mr = m2; mr -= m1;
    EXPECT_EQ(mod(16), mr);
    mr = m2; mr *= m1;
    EXPECT_EQ(mod(-63), mr);
    mr = m2; mr /= m1;
    EXPECT_EQ(mod(INT64_C(5270498306774157589)), mr);
    mr = m2; mr %= m1;
    EXPECT_EQ(mod(9), mr);
    mr = m3; mr /= m1;
    EXPECT_EQ(mod(3), m3 / m1);
    mr = m1; mr /= m3;
    EXPECT_EQ(mod(INT64_C(6148914691236517189)), m1 / m3);
}

TEST(modulo_int64_test, operators_inplace_self) {
    const int64_t M = INT64_C(9223372036854775783);
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

TEST(modulo_int64_test, casts) {
    const int64_t M = INT64_C(9223372036854775783);
    const mod m1 = -7;
    const mod e0 = zeroOf(m1);
    const mod e1 = identityOf(m1);
    EXPECT_EQ(INT64_C(0), e0.v);
    EXPECT_EQ(M, e0.M());
    EXPECT_EQ(INT64_C(1), e1.v);
    EXPECT_EQ(M, e1.M());
    const mod m3 = castOf<mod>(UINT64_C(10000000000000000000));
    EXPECT_EQ(INT64_C(776627963145224217), m3.v);
    EXPECT_EQ(M, m3.M());
    const mod m5 = castOf(m1, -5);
    EXPECT_EQ(INT64_C(9223372036854775778), m5.v);
    EXPECT_EQ(M, m5.M());
    const mod m6 = castOf(m1, m5);
    EXPECT_EQ(INT64_C(9223372036854775778), m6.v);
    EXPECT_EQ(M, m6.M());
    const mod m7 = castOf<mod>(m5);
    EXPECT_EQ(INT64_C(9223372036854775778), m7.v);
    EXPECT_EQ(M, m7.M());
    EXPECT_EQ(INT64_C(4), modT(INT64_C(9223372036854775787), M));
    const mod m8 = powT(m1, 100);
    EXPECT_EQ(INT64_C(9175964761415298625), m8.v);
    EXPECT_EQ(M, m8.M());
}
