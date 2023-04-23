#include "altruct/structure/math/modulo.h"
#include "structure_test_util.h"

#include "gtest/gtest.h"

#include <functional>

using namespace std;
using namespace altruct::math;
using namespace altruct::test_util;

// largest prime that fits uint32_t: 4294967291 = 2^32 - 5
typedef modulo<uint32_t, 4294967291, modulo_storage::CONSTANT> mod;

TEST(modulo_uint32_test, standalone_functions_1000000007) {
    const uint32_t O = UINT32_C(0);
    const uint32_t M = UINT32_C(1000000007);
    EXPECT_EQ(O + 0u, modulo_normalize(INT32_C(-2000000014), M));
    EXPECT_EQ(O + 0u, modulo_normalize(O + 0u, M));
    EXPECT_EQ(O + 12u, modulo_normalize(UINT32_C(4000000040), M));
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
    EXPECT_EQ(UINT32_C(500000004), modulo_inv(O + 2u, M));
    EXPECT_EQ(O + 2u, modulo_inv(UINT32_C(500000004), M));
    EXPECT_EQ(UINT32_C(333333336), modulo_inv(O + 3u, M));
    EXPECT_EQ(O + 3u, modulo_inv(UINT32_C(333333336), M));
    EXPECT_EQ(O + 0u, modulo_div(O + 0u, O + 7u, M));
    EXPECT_EQ(O + 7u, modulo_div(O + 7u, O + 1u, M));
    EXPECT_EQ(UINT32_C(428571432), modulo_div(O + 3u, O + 7u, M));
    EXPECT_EQ(O + 7u, modulo_div(O + 3u, UINT32_C(428571432), M));
}

TEST(modulo_uint32_test, standalone_functions_4000000007) {
    EXPECT_EQ(UINT32_C(0), modulo_normalize(UINT32_C(0), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(3000000030), modulo_normalize(UINT32_C(3000000030), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(4000000006), modulo_add(UINT32_C(4000000004), UINT32_C(2), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(4000000002), modulo_add(UINT32_C(4000000004), UINT32_C(4000000005), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(13), modulo_add(UINT32_C(13), UINT32_C(4000000007), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(4000000002), modulo_sub(UINT32_C(4000000004), UINT32_C(2), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(4000000006), modulo_sub(UINT32_C(4000000004), UINT32_C(4000000005), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(14), modulo_sub(UINT32_C(14), UINT32_C(4000000007), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(0), modulo_neg(UINT32_C(0), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(4000000005), modulo_neg(UINT32_C(2), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(3), modulo_neg(UINT32_C(4000000004), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(15), modulo_mul(UINT32_C(3), UINT32_C(5), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(4000000001), modulo_mul(UINT32_C(3), UINT32_C(4000000005), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(4000000001), modulo_mul(UINT32_C(4000000004), UINT32_C(2), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(18), modulo_mul(UINT32_C(4000000004), UINT32_C(4000000001), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(1), modulo_inv(UINT32_C(1), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(4000000006), modulo_inv(UINT32_C(4000000006), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(2000000004), modulo_inv(UINT32_C(2), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(2), modulo_inv(UINT32_C(2000000004), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(1333333336), modulo_inv(UINT32_C(3), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(3), modulo_inv(UINT32_C(1333333336), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(0), modulo_div(UINT32_C(0), UINT32_C(7), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(7), modulo_div(UINT32_C(7), UINT32_C(1), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(3428571435), modulo_div(UINT32_C(3), UINT32_C(7), UINT32_C(4000000007)));
    EXPECT_EQ(UINT32_C(7), modulo_div(UINT32_C(3), UINT32_C(3428571435), UINT32_C(4000000007)));
}

TEST(modulo_uint32_test, modulo_gcd_ex) {
    uint32_t ni1, ni2;
    modulo_gcd_ex(UINT32_C(2971215073), UINT32_C(4294930221), ni1, ni2);
    EXPECT_EQ(UINT32_C(367514362), ni1);
    EXPECT_EQ(UINT32_C(2716970148), ni2);
    modulo_gcd_ex(UINT32_C(4294930221), UINT32_C(2971215073), ni1, ni2);
    EXPECT_EQ(UINT32_C(2716970148), ni1);
    EXPECT_EQ(UINT32_C(367514362), ni2);
}

TEST(modulo_uint32_test, constructor) {
    const uint32_t M = UINT32_C(4294967291);
    // default
    const mod m1;
    EXPECT_EQ(UINT32_C(0), m1.v);
    EXPECT_EQ(M, m1.M());
    // value only
    const mod m2(UINT32_C(10));
    EXPECT_EQ(UINT32_C(10), m2.v);
    EXPECT_EQ(M, m2.M());
    // value + modulus, modulus ignored
    const mod m3(UINT32_C(13), UINT32_C(12345));
    EXPECT_EQ(UINT32_C(13), m3.v);
    EXPECT_EQ(M, m3.M());

    // from same integral type: uint32_t
    const mod mu32_0(UINT32_C(0));
    EXPECT_EQ(UINT32_C(0), mu32_0.v);
    EXPECT_EQ(M, mu32_0.M());
    const mod mu32_1(UINT32_C(10));
    EXPECT_EQ(UINT32_C(10), mu32_1.v);
    EXPECT_EQ(M, mu32_1.M());
    const mod mu32_2(UINT32_C(4294967290)); // -1
    EXPECT_EQ(UINT32_C(4294967290), mu32_2.v);
    EXPECT_EQ(M, mu32_2.M());
    const mod mu32_3(UINT32_C(4294967292)); // +1
    EXPECT_EQ(UINT32_C(1), mu32_3.v);
    EXPECT_EQ(M, mu32_3.M());

    // from different integral type: int32_t
    const mod mi32_0(INT32_C(0));
    EXPECT_EQ(UINT32_C(0), mi32_0.v);
    EXPECT_EQ(M, mi32_0.M());
    const mod mi32_1(INT32_C(20));
    EXPECT_EQ(UINT32_C(20), mi32_1.v);
    EXPECT_EQ(M, mi32_1.M());
    const mod mi32_2(INT32_C(-2));
    EXPECT_EQ(UINT32_C(4294967289), mi32_2.v);
    EXPECT_EQ(M, mi32_2.M());
    const mod mi32_3(INT32_C(-102));
    EXPECT_EQ(UINT32_C(4294967189), mi32_3.v);
    EXPECT_EQ(M, mi32_3.M());

    // from different integral type: uint64_t
    const mod mu64_0(UINT64_C(0));
    EXPECT_EQ(UINT32_C(0), mu64_0.v);
    EXPECT_EQ(M, mu64_0.M());
    const mod mu64_1(UINT64_C(40));
    EXPECT_EQ(UINT32_C(40), mu64_1.v);
    EXPECT_EQ(M, mu64_1.M());
    const mod mu64_2(UINT64_C(4294967287)); // -4
    EXPECT_EQ(UINT32_C(4294967287), mu64_2.v);
    EXPECT_EQ(M, mu64_2.M());
    const mod mu64_3(UINT64_C(4294967187)); // -104
    EXPECT_EQ(UINT32_C(4294967187), mu64_3.v);
    EXPECT_EQ(M, mu64_3.M());
    const mod mu64_4(UINT64_C(4294967295)); // 4
    EXPECT_EQ(UINT32_C(4), mu64_4.v);
    EXPECT_EQ(M, mu64_4.M());
    const mod mu64_5(UINT64_C(1000000000000));
    EXPECT_EQ(UINT32_C(3567588488), mu64_5.v);
    EXPECT_EQ(M, mu64_5.M());

    // from different integral type: int64_t
    const mod mi64_0(INT64_C(0));
    EXPECT_EQ(UINT32_C(0), mi64_0.v);
    EXPECT_EQ(M, mi64_0.M());
    const mod mi64_1(INT64_C(50));
    EXPECT_EQ(UINT32_C(50), mi64_1.v);
    EXPECT_EQ(M, mi64_1.M());
    const mod mi64_2(INT64_C(-5));
    EXPECT_EQ(UINT32_C(4294967286), mi64_2.v);
    EXPECT_EQ(M, mi64_2.M());
    const mod mi64_3(INT64_C(-105));
    EXPECT_EQ(UINT32_C(4294967186), mi64_3.v);
    EXPECT_EQ(M, mi64_3.M());
    const mod mi64_4(INT64_C(4294967296));
    EXPECT_EQ(UINT32_C(5), mi64_4.v);
    EXPECT_EQ(M, mi64_4.M());
    const mod mi64_5(INT64_C(1000000000000));
    EXPECT_EQ(UINT32_C(3567588488), mi64_5.v);
    EXPECT_EQ(M, mi64_5.M());
    const mod mi64_6(INT64_C(-1000000000000));
    EXPECT_EQ(UINT32_C(727378803), mi64_6.v);
    EXPECT_EQ(M, mi64_6.M());

    // from different integral type: int64_t
    // value + modulus, modulus ignored
    const mod mi64_7(INT64_C(-1000000000000), UINT32_C(12345));
    EXPECT_EQ(UINT32_C(727378803), mi64_7.v);
    EXPECT_EQ(M, mi64_7.M());

    // copy constructor
    const mod mu32_c(mu32_1);
    EXPECT_EQ(UINT32_C(10), mu32_c.v);
    EXPECT_EQ(M, mu32_c.M());
    // move constructor
    const mod mu32_m(std::move(mu32_2));
    EXPECT_EQ(UINT32_C(4294967290), mu32_m.v);
    EXPECT_EQ(M, mu32_m.M());
    // assignment
    mod mu32_a; mu32_a = mu32_1;
    EXPECT_EQ(UINT32_C(10), mu32_a.v);
    EXPECT_EQ(M, mu32_a.M());
    // move assignment
    mu32_a = std::move(mu32_3);
    EXPECT_EQ(UINT32_C(1), mu32_a.v);
    EXPECT_EQ(M, mu32_a.M());
}

TEST(modulo_uint32_test, operators_comparison) {
    const mod m1 = 10;
    const mod m2 = 20;
    ASSERT_COMPARISON_OPERATORS(0, m1, m1);
    ASSERT_COMPARISON_OPERATORS(0, m2, m2);
    ASSERT_COMPARISON_OPERATORS(-1, m1, m2);
    ASSERT_COMPARISON_OPERATORS(+1, m2, m1);
}

TEST(modulo_uint32_test, operators_arithmetic) {
    const uint32_t M = UINT32_C(4294967291);
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
    EXPECT_EQ(mod(UINT32_C(954437175)), m1 / m2);
    EXPECT_EQ(mod(1), m1 % m2);
    EXPECT_EQ(mod(2), m2 + m1);
    EXPECT_EQ(mod(16), m2 - m1);
    EXPECT_EQ(mod(-9), -m2);
    EXPECT_EQ(mod(-63), m2 * m1);
    EXPECT_EQ(mod(UINT32_C(3067833778)), m2 / m1);
    EXPECT_EQ(mod(9), m2 % m1);
    EXPECT_EQ(mod(3), m3 / m1);
    EXPECT_EQ(mod(UINT32_C(1431655764)), m1 / m3);
}

TEST(modulo_uint32_test, operators_inplace) {
    const uint32_t M = UINT32_C(4294967291);
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
    EXPECT_EQ(mod(UINT32_C(954437175)), mr);
    mr = m1; mr %= m2;
    EXPECT_EQ(mod(1), mr);
    mr = m2; mr += m1;
    EXPECT_EQ(mod(2), mr);
    mr = m2; mr -= m1;
    EXPECT_EQ(mod(16), mr);
    mr = m2; mr *= m1;
    EXPECT_EQ(mod(-63), mr);
    mr = m2; mr /= m1;
    EXPECT_EQ(mod(UINT32_C(3067833778)), mr);
    mr = m2; mr %= m1;
    EXPECT_EQ(mod(9), mr);
    mr = m3; mr /= m1;
    EXPECT_EQ(mod(3), m3 / m1);
    mr = m1; mr /= m3;
    EXPECT_EQ(mod(UINT32_C(1431655764)), m1 / m3);
}

TEST(modulo_uint32_test, operators_inplace_self) {
    const uint32_t M = UINT32_C(4294967291);
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

TEST(modulo_uint32_test, casts) {
    const uint32_t M = UINT32_C(4294967291);
    const mod m1 = -7;
    const mod e0 = zeroOf(m1);
    const mod e1 = identityOf(m1);
    EXPECT_EQ(UINT32_C(0), e0.v);
    EXPECT_EQ(M, e0.M());
    EXPECT_EQ(UINT32_C(1), e1.v);
    EXPECT_EQ(M, e1.M());
    const mod m3 = castOf<mod>(INT64_C(1000000000000));
    EXPECT_EQ(UINT32_C(3567588488), m3.v);
    EXPECT_EQ(M, m3.M());
    const mod m5 = castOf(m1, -5);
    EXPECT_EQ(UINT32_C(4294967286), m5.v);
    EXPECT_EQ(M, m5.M());
    const mod m6 = castOf(m1, m5);
    EXPECT_EQ(UINT32_C(4294967286), m6.v);
    EXPECT_EQ(M, m6.M());
    const mod m7 = castOf<mod>(m5);
    EXPECT_EQ(UINT32_C(4294967286), m7.v);
    EXPECT_EQ(M, m7.M());
    EXPECT_EQ(UINT32_C(4), modT(UINT32_C(4294967295), M));
    const mod m8 = powT(m1, 10);
    EXPECT_EQ(UINT32_C(282475249), m8.v);
    EXPECT_EQ(M, m8.M());
}

TEST(modulo_uint8_test, modulo_normalize_bruteforce) {
    for (int m = 1; m < (1 << 8); m++) {
        for (int v = -(1 << 7); v < (1 << 7); v++) {
            uint8_t vn0 = ((v % m) + m) % m;
            uint8_t vn = modulo_normalize(int8_t(v), m);
            EXPECT_EQ(vn0, vn) << int(vn) << " != " << int(v) << " % " << m;
        }
    }
    for (int m = 1; m < (1 << 8); m++) {
        for (int v = 0; v < (1 << 8); v++) {
            uint8_t vn0 = ((v % m) + m) % m;
            uint8_t vn = modulo_normalize(uint8_t(v), m);
            EXPECT_EQ(vn0, vn) << int(vn) << " != " << int(v) << " % " << m;
        }
    }
}

TEST(modulo_uint8_test, modulo_inv_int_bruteforce) {
    for (int m = 1; m < (1 << 8); m++) {
        for (int v = 1; v < m; v++) {
            if (gcd(m, v) != 1) continue;
            uint8_t vi = modulo_inv_int<uint8_t>(v, m);
            EXPECT_TRUE(vi < m);
            uint8_t e = (uint16_t(v) * vi) % uint8_t(m);
            EXPECT_EQ(1, e) << v << " * " << int(vi) << " != 1  mod " << m;
        }
    }
}

TEST(modulo_uint16_test, modulo_inv_int_bruteforce) {
    for (int m = 1; m < (1 << 16); m += 1000) { // step 1000 for speed
        for (int v = 1; v < m; v++) {
            if (gcd(m, v) != 1) continue;
            uint16_t vi = modulo_inv_int<uint16_t>(v, m);
            EXPECT_TRUE(vi < m);
            uint16_t e = (uint32_t(v) * vi) % uint16_t(m);
            EXPECT_EQ(1, e) << v << " * " << int(vi) << " != 1  mod " << m;
        }
    }
}
