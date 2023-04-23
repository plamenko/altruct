#include "altruct/structure/math/modulo.h"
#include "structure_test_util.h"

#include "gtest/gtest.h"

#include <functional>

using namespace std;
using namespace altruct::math;
using namespace altruct::test_util;

// second largest prime that fits uint32_t: 2147483629 = 2^31 - 19
typedef modulo<int32_t, 2147483629, modulo_storage::CONSTANT> mod;

TEST(modulo_int32_test, standalone_functions_1000000007) {
    const int32_t O = INT32_C(0);
    const int32_t M = INT32_C(1000000007);
    EXPECT_EQ(O + 0, modulo_normalize(INT32_C(-2000000014), M));
    EXPECT_EQ(O + 0, modulo_normalize(O + 0, M));
    EXPECT_EQ(O + 12, modulo_normalize(UINT32_C(4000000040), M));
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
    EXPECT_EQ(INT32_C(500000004), modulo_inv(O + 2, M));
    EXPECT_EQ(O + 2, modulo_inv(INT32_C(500000004), M));
    EXPECT_EQ(INT32_C(333333336), modulo_inv(O + 3, M));
    EXPECT_EQ(O + 3, modulo_inv(INT32_C(333333336), M));
    EXPECT_EQ(O + 0, modulo_div(O + 0, O + 7, M));
    EXPECT_EQ(O + 7, modulo_div(O + 7, O + 1, M));
    EXPECT_EQ(INT32_C(428571432), modulo_div(O + 3, O + 7, M));
    EXPECT_EQ(O + 7, modulo_div(O + 3, INT32_C(428571432), M));
}

TEST(modulo_int32_test, standalone_functions_2000000011) {
    EXPECT_EQ(INT32_C(0), modulo_normalize(INT32_C(0), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(2000000002), modulo_normalize(INT32_C(-2000000020), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(2000000006), modulo_add(INT32_C(2000000004), INT32_C(2), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(2000000002), modulo_add(INT32_C(2000000004), INT32_C(2000000009), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(13), modulo_add(INT32_C(13), INT32_C(2000000011), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(2000000002), modulo_sub(INT32_C(2000000004), INT32_C(2), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(2000000010), modulo_sub(INT32_C(2000000004), INT32_C(2000000005), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(14), modulo_sub(INT32_C(14), INT32_C(2000000011), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(0), modulo_neg(INT32_C(0), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(2000000009), modulo_neg(INT32_C(2), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(3), modulo_neg(INT32_C(2000000008), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(15), modulo_mul(INT32_C(3), INT32_C(5), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(2000000005), modulo_mul(INT32_C(3), INT32_C(2000000009), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(2000000005), modulo_mul(INT32_C(2000000008), INT32_C(2), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(18), modulo_mul(INT32_C(2000000008), INT32_C(2000000005), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(1), modulo_inv(INT32_C(1), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(2000000010), modulo_inv(INT32_C(2000000010), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(1000000006), modulo_inv(INT32_C(2), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(2), modulo_inv(INT32_C(1000000006), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(1333333341), modulo_inv(INT32_C(3), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(3), modulo_inv(INT32_C(1333333341), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(0), modulo_div(INT32_C(0), INT32_C(7), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(7), modulo_div(INT32_C(7), INT32_C(1), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(571428575), modulo_div(INT32_C(3), INT32_C(7), INT32_C(2000000011)));
    EXPECT_EQ(INT32_C(7), modulo_div(INT32_C(3), INT32_C(571428575), INT32_C(2000000011)));
}

TEST(modulo_int32_test, modulo_gcd_ex) {
    int32_t ni1, ni2;
    modulo_gcd_ex(1134903170, 1836311903, ni1, ni2);
    EXPECT_EQ(1134903170, ni1);
    EXPECT_EQ(433494437, ni2);
    modulo_gcd_ex(1836311903, 1134903170, ni1, ni2);
    EXPECT_EQ(433494437, ni1);
    EXPECT_EQ(1134903170, ni2);
    modulo_gcd_ex(2147450880, 1836311903, ni1, ni2);
    EXPECT_EQ(459437288, ni1);
    EXPECT_EQ(1610167967, ni2);
    modulo_gcd_ex(1836311903, 2147450880, ni1, ni2);
    EXPECT_EQ(1610167967, ni1);
    EXPECT_EQ(459437288, ni2);
}

TEST(modulo_int32_test, constructor) {
    const int32_t M = INT32_C(2147483629);
    // default
    const mod m1;
    EXPECT_EQ(INT32_C(0), m1.v);
    EXPECT_EQ(M, m1.M());
    // value only
    const mod m2(INT32_C(10));
    EXPECT_EQ(INT32_C(10), m2.v);
    EXPECT_EQ(M, m2.M());
    // value + modulus, modulus ignored
    const mod m3(INT32_C(13), INT32_C(12345));
    EXPECT_EQ(INT32_C(13), m3.v);
    EXPECT_EQ(M, m3.M());

    // from different integral type: uint32_t
    const mod mu32_0(UINT32_C(0));
    EXPECT_EQ(INT32_C(0), mu32_0.v);
    EXPECT_EQ(M, mu32_0.M());
    const mod mu32_1(UINT32_C(10));
    EXPECT_EQ(INT32_C(10), mu32_1.v);
    EXPECT_EQ(M, mu32_1.M());
    const mod mu32_2(UINT32_C(2147483628)); // -1
    EXPECT_EQ(INT32_C(2147483628), mu32_2.v);
    EXPECT_EQ(M, mu32_2.M());
    const mod mu32_3(UINT32_C(2147483630)); // +1
    EXPECT_EQ(INT32_C(1), mu32_3.v);
    EXPECT_EQ(M, mu32_3.M());

    // from same integral type: int32_t
    const mod mi32_0(INT32_C(0));
    EXPECT_EQ(INT32_C(0), mi32_0.v);
    EXPECT_EQ(M, mi32_0.M());
    const mod mi32_1(INT32_C(20));
    EXPECT_EQ(INT32_C(20), mi32_1.v);
    EXPECT_EQ(M, mi32_1.M());
    const mod mi32_2(INT32_C(-2));
    EXPECT_EQ(INT32_C(2147483627), mi32_2.v);
    EXPECT_EQ(M, mi32_2.M());
    const mod mi32_3(INT32_C(-102));
    EXPECT_EQ(INT32_C(2147483527), mi32_3.v);
    EXPECT_EQ(M, mi32_3.M());

    // from different integral type: uint64_t
    const mod mu64_0(UINT64_C(0));
    EXPECT_EQ(INT32_C(0), mu64_0.v);
    EXPECT_EQ(M, mu64_0.M());
    const mod mu64_1(UINT64_C(40));
    EXPECT_EQ(INT32_C(40), mu64_1.v);
    EXPECT_EQ(M, mu64_1.M());
    const mod mu64_2(UINT64_C(4294967254)); // -4
    EXPECT_EQ(INT32_C(2147483625), mu64_2.v);
    EXPECT_EQ(M, mu64_2.M());
    const mod mu64_3(UINT64_C(4294967154)); // -104
    EXPECT_EQ(INT32_C(2147483525), mu64_3.v);
    EXPECT_EQ(M, mu64_3.M());
    const mod mu64_4(UINT64_C(4294967262)); // 4
    EXPECT_EQ(INT32_C(4), mu64_4.v);
    EXPECT_EQ(M, mu64_4.M());
    const mod mu64_5(UINT64_C(1000000000000));
    EXPECT_EQ(INT32_C(1420112515), mu64_5.v);
    EXPECT_EQ(M, mu64_5.M());

    // from different integral type: int64_t
    const mod mi64_0(INT64_C(0));
    EXPECT_EQ(INT32_C(0), mi64_0.v);
    EXPECT_EQ(M, mi64_0.M());
    const mod mi64_1(INT64_C(50));
    EXPECT_EQ(INT32_C(50), mi64_1.v);
    EXPECT_EQ(M, mi64_1.M());
    const mod mi64_2(INT64_C(-5));
    EXPECT_EQ(INT32_C(2147483624), mi64_2.v);
    EXPECT_EQ(M, mi64_2.M());
    const mod mi64_3(INT64_C(-105));
    EXPECT_EQ(INT32_C(2147483524), mi64_3.v);
    EXPECT_EQ(M, mi64_3.M());
    const mod mi64_4(INT64_C(4294967296));
    EXPECT_EQ(INT32_C(38), mi64_4.v);
    EXPECT_EQ(M, mi64_4.M());
    const mod mi64_5(INT64_C(1000000000000));
    EXPECT_EQ(INT32_C(1420112515), mi64_5.v);
    EXPECT_EQ(M, mi64_5.M());
    const mod mi64_6(INT64_C(-1000000000000));
    EXPECT_EQ(INT32_C(727371114), mi64_6.v);
    EXPECT_EQ(M, mi64_6.M());

    // from different integral type: uint64_t
    // value + modulus, modulus ignored
    const mod mu64_7(UINT64_C(1000000000000), INT32_C(12345));
    EXPECT_EQ(INT32_C(1420112515), mu64_7.v);
    EXPECT_EQ(M, mu64_7.M());

    // copy constructor
    const mod mi32_c(mi32_1);
    EXPECT_EQ(INT32_C(20), mi32_c.v);
    EXPECT_EQ(M, mi32_c.M());
    // move constructor
    const mod mi32_m(std::move(mi32_2));
    EXPECT_EQ(INT32_C(2147483627), mi32_m.v);
    EXPECT_EQ(M, mi32_m.M());
    // assignment
    mod mi32_a; mi32_a = mi32_1;
    EXPECT_EQ(INT32_C(20), mi32_a.v);
    EXPECT_EQ(M, mi32_a.M());
    // move assignment
    mi32_a = std::move(mi32_3);
    EXPECT_EQ(INT32_C(2147483527), mi32_a.v);
    EXPECT_EQ(M, mi32_a.M());
}

TEST(modulo_int32_test, operators_comparison) {
    const mod m1 = 10;
    const mod m2 = 20;
    ASSERT_COMPARISON_OPERATORS(0, m1, m1);
    ASSERT_COMPARISON_OPERATORS(0, m2, m2);
    ASSERT_COMPARISON_OPERATORS(-1, m1, m2);
    ASSERT_COMPARISON_OPERATORS(+1, m2, m1);
}

TEST(modulo_int32_test, operators_arithmetic) {
    const int32_t M = INT32_C(2147483629);
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
    EXPECT_EQ(mod(INT32_C(1670265044)), m1 / m2);
    EXPECT_EQ(mod(3), m1 % m2);
    EXPECT_EQ(mod(2), m2 + m1);
    EXPECT_EQ(mod(16), m2 - m1);
    EXPECT_EQ(mod(-9), -m2);
    EXPECT_EQ(mod(-63), m2 * m1);
    EXPECT_EQ(mod(INT32_C(1227133501)), m2 / m1);
    EXPECT_EQ(mod(9), m2 % m1);
    EXPECT_EQ(mod(3), m3 / m1);
    EXPECT_EQ(mod(INT32_C(1431655753)), m1 / m3);
}

TEST(modulo_int32_test, operators_inplace) {
    const int32_t M = INT32_C(2147483629);
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
    EXPECT_EQ(mod(INT32_C(1670265044)), mr);
    mr = m1; mr %= m2;
    EXPECT_EQ(mod(3), mr);
    mr = m2; mr += m1;
    EXPECT_EQ(mod(2), mr);
    mr = m2; mr -= m1;
    EXPECT_EQ(mod(16), mr);
    mr = m2; mr *= m1;
    EXPECT_EQ(mod(-63), mr);
    mr = m2; mr /= m1;
    EXPECT_EQ(mod(INT32_C(1227133501)), mr);
    mr = m2; mr %= m1;
    EXPECT_EQ(mod(9), mr);
    mr = m3; mr /= m1;
    EXPECT_EQ(mod(3), m3 / m1);
    mr = m1; mr /= m3;
    EXPECT_EQ(mod(INT32_C(1431655753)), m1 / m3);
}

TEST(modulo_int32_test, operators_inplace_self) {
    const int32_t M = INT32_C(2147483629);
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

TEST(modulo_int32_test, casts) {
    const int32_t M = INT32_C(2147483629);
    const mod m1 = -7;
    const mod e0 = zeroOf(m1);
    const mod e1 = identityOf(m1);
    EXPECT_EQ(INT32_C(0), e0.v);
    EXPECT_EQ(M, e0.M());
    EXPECT_EQ(INT32_C(1), e1.v);
    EXPECT_EQ(M, e1.M());
    const mod m3 = castOf<mod>(INT64_C(1000000000000));
    EXPECT_EQ(INT32_C(1420112515), m3.v);
    EXPECT_EQ(M, m3.M());
    const mod m5 = castOf(m1, -5);
    EXPECT_EQ(INT32_C(2147483624), m5.v);
    EXPECT_EQ(M, m5.M());
    const mod m6 = castOf(m1, m5);
    EXPECT_EQ(INT32_C(2147483624), m6.v);
    EXPECT_EQ(M, m6.M());
    const mod m7 = castOf<mod>(m5);
    EXPECT_EQ(INT32_C(2147483624), m7.v);
    EXPECT_EQ(M, m7.M());
    EXPECT_EQ(INT32_C(4), modT(INT32_C(2147483633), M));
    const mod m8 = powT(m1, 100);
    EXPECT_EQ(INT32_C(681305249), m8.v);
    EXPECT_EQ(M, m8.M());
}

TEST(modulo_int8_test, modulo_normalize_bruteforce) {
    for (int m = 1; m < (1 << 7); m++) {
        for (int v = -(1 << 7); v < (1 << 7); v++) {
            int8_t vn0 = ((v % m) + m) % m;
            int8_t vn = modulo_normalize(int8_t(v), m);
            EXPECT_EQ(vn0, vn) << int(vn) << " != " << int(v) << " % " << m;
        }
    }
    for (int m = 1; m < (1 << 7); m++) {
        for (int v = 0; v < (1 << 8); v++) {
            int8_t vn0 = ((v % m) + m) % m;
            int8_t vn = modulo_normalize(uint8_t(v), m);
            EXPECT_EQ(vn0, vn) << int(vn) << " != " << int(v) << " % " << m;
        }
    }
}

TEST(modulo_int8_test, modulo_inv_int_bruteforce) {
    for (int m = 1; m < (1 << 7); m++) {
        for (int v = 1; v < m; v++) {
            if (gcd(m, v) != 1) continue;
            int8_t vi = modulo_inv_int<int8_t>(v, m);
            EXPECT_TRUE(vi < m);
            int8_t e = (int16_t(v) * vi) % int8_t(m);
            EXPECT_EQ(1, e) << v << " * " << int(vi) << " != 1  mod " << m;
        }
    }
}

TEST(modulo_int16_test, modulo_inv_int_bruteforce) {
    for (int m = 1; m < (1 << 15); m += 1000) { // step 1000 for speed
        for (int v = 1; v < m; v++) {
            if (gcd(m, v) != 1) continue;
            int16_t vi = modulo_inv_int<int16_t>(v, m);
            EXPECT_TRUE(vi < m);
            int16_t e = (int32_t(v) * vi) % int16_t(m);
            EXPECT_EQ(1, e) << v << " * " << int(vi) << " != 1  mod " << m;
        }
    }
}
