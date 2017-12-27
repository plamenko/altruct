#include "altruct/algorithm/collections/collections.h"
#include "altruct/structure/math/series.h"
#include "altruct/structure/math/modulo.h"
#include "structure_test_util.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;
using namespace altruct::collections;
using namespace altruct::test_util;

namespace {
typedef moduloX<int> modx;
typedef polynom<modx> polyx;
typedef seriesX<modx> serx;

vector<moduloX<int>> to_modx(int M, const vector<int>& v) {
    return transform(v, [&](int a) {
        return moduloX<int>(a, M);
    });
}

polyx make_polyx(int M, const vector<int>& v) {
    polyx p(to_modx(M, v));
    p.ZERO_COEFF = modx(0, M);
    return p;
}

serx make_serx(int M, const vector<int>& v) {
    return serx(make_polyx(M, v));
}
}

TEST(series_modx_test, constructor) {
    const auto p = make_polyx(1009, { 1, 2, 3, 4 });
    serx s0;
    EXPECT_EQ(make_polyx(1009, {}), s0.p);
    EXPECT_EQ(1, s0.N());
    serx s1(5);
    EXPECT_EQ(make_polyx(1009, { 5 }), s1.p);
    EXPECT_EQ(1, s1.N());
    serx s2(p);
    EXPECT_EQ(p, s2.p);
    EXPECT_EQ(4, s2.N());
    serx s3(s2);
    EXPECT_EQ(p, s3.p);
    EXPECT_EQ(4, s3.N());
    serx s4(p.c.begin(), p.c.end());
    EXPECT_EQ(p, s4.p);
    EXPECT_EQ(4, s4.N());
    serx s5(&p[0], &p[0] + p.size());
    EXPECT_EQ(p, s5.p);
    EXPECT_EQ(4, s5.N());
    serx s7{ modx(1, 1009), modx(2, 1009), modx(3, 1009), modx(4, 1009) };
    EXPECT_EQ(p, s7.p);
    EXPECT_EQ(4, s7.N());
    serx s8(to_modx(1009, { 1, 2, 3, 4 }));
    EXPECT_EQ(p, s8.p);
    EXPECT_EQ(4, s8.N());
    serx s9(make_polyx(1009, { 1, 2, 3, 4 }));
    EXPECT_EQ(p, s9.p);
    EXPECT_EQ(4, s9.N());
    serx s10(serx(p.c.begin(), p.c.end()));
    EXPECT_EQ(p, s10.p);
    EXPECT_EQ(4, s10.N());
    serx s11 = make_serx(1009, { 1, 2, 3, 4 });
    EXPECT_EQ(p, s11.p);
    EXPECT_EQ(4, s11.N());
    serx s12 = s11;
    EXPECT_EQ(p, s12.p);
    EXPECT_EQ(4, s12.N());
}

TEST(series_modx_test, swap) {
    auto s1 = make_serx(1009, { 1, 2, 3, 4 });
    auto s2 = make_serx(1009, { 5, 6, 7 });
    s1.swap(s2);
    EXPECT_EQ(make_polyx(1009, { 5, 6, 7 }), s1.p);
    EXPECT_EQ(3, s1.N());
    EXPECT_EQ(make_polyx(1009, { 1, 2, 3, 4 }), s2.p);
    EXPECT_EQ(4, s2.N());
}

TEST(series_modx_test, resize) {
    auto s = make_serx(1009, {1, 2, 3, 4});
    s.resize(6);
    EXPECT_EQ(make_polyx(1009, { 1, 2, 3, 4, 0, 0 }), s.p);
    EXPECT_EQ(6, s.N());
    s.resize(3);
    EXPECT_EQ(make_polyx(1009, { 1, 2, 3 }), s.p);
    EXPECT_EQ(3, s.N());
}

TEST(series_modx_test, size) {
    auto s = make_serx(1009, { 1, 2, 3 });
    EXPECT_EQ(3, s.size());
}

TEST(series_modx_test, at) {
    auto s = make_serx(1009, { 2, 3, 5, 7 });
    EXPECT_EQ(modx(2, 1009), s.at(0));
    EXPECT_EQ(modx(3, 1009), s.at(1));
    EXPECT_EQ(modx(5, 1009), s.at(2));
    EXPECT_EQ(modx(7, 1009), s.at(3));
    EXPECT_EQ(modx(0, 1009), s.at(4));
    EXPECT_EQ(modx(0, 1009), s.at(100));
    EXPECT_EQ(4, s.size());
}

TEST(series_modx_test, operator_const_brackets) {
    const auto s = make_serx(1009, { 2, 3, 5, 7 });
    EXPECT_EQ(modx(2, 1009), s[0]);
    EXPECT_EQ(modx(7, 1009), s[3]);
    EXPECT_EQ(modx(0, 1009), s[4]);
    EXPECT_EQ(modx(0, 1009), s[100]);
    EXPECT_EQ(modx(4, 1009), s.size());
}

TEST(series_modx_test, operator_brackets) {
    auto s = make_serx(1009, {});
    s.resize(4);
    s[3] = modx(3, 1009);
    EXPECT_EQ(modx(0, 1009), s[0]);
    EXPECT_EQ(modx(0, 1009), s[4]);
    EXPECT_EQ(modx(3, 1009), s[3]);
    EXPECT_EQ(modx(0, 1009), s[4]);
    EXPECT_EQ(modx(0, 1009), s[100]);
    EXPECT_EQ(4, s.size());
}

TEST(series_modx_test, operators_comparison) {
    const auto s1 = make_serx(1009, { 4 });
    const auto s2 = make_serx(1009, { 1, 3, 5, 7 });
    const auto s3 = make_serx(1009, { 1, 3, 5, 7, 0, 0, 0 });
    ASSERT_COMPARISON_OPERATORS(0, s1, s1);
    ASSERT_COMPARISON_OPERATORS(0, s2, s2);
    ASSERT_COMPARISON_OPERATORS(0, s3, s3);
    ASSERT_COMPARISON_OPERATORS(-1, s1, s2);
    ASSERT_COMPARISON_OPERATORS(+1, s2, s1);
    ASSERT_COMPARISON_OPERATORS(0, s2, s3);
    ASSERT_COMPARISON_OPERATORS(0, s3, s2);
}

TEST(series_modx_test, inverse) {
    const auto s = make_serx(1009, { 1, -3, 5, 7 });
    const auto si = s.inverse();
    EXPECT_EQ(make_polyx(1009, { 1, 3, 4, -10 }), si.p);
    EXPECT_EQ(4, si.N());
}

TEST(series_modx_test, operators_arithmetic) {
    const auto s1 = make_serx(1009, { 4, 1, 0, 0 });
    const auto s2 = make_serx(1009, { 1, -3, 5, 7 });
    EXPECT_EQ(make_serx(1009, { 5, -2, 5, 7 }), s1 + s2);
    EXPECT_EQ(make_serx(1009, { 3, 4, -5, -7 }), s1 - s2);
    EXPECT_EQ(make_serx(1009, { -1, 3, -5, -7 }), -s2);
    EXPECT_EQ(make_serx(1009, { 4, -11, 17, 33 }), s1 * s2);
    EXPECT_EQ(make_serx(1009, { 4, 13, 19, -36 }), s1 / s2);
    EXPECT_EQ(make_serx(1009, { 5, -2, 5, 7 }), s2 + s1);
    EXPECT_EQ(make_serx(1009, { -3, -4, 5, 7 }), s2 - s1);
    EXPECT_EQ(make_serx(1009, { -4, -1 }), -s1);
    EXPECT_EQ(make_serx(1009, { 4, -11, 17, 33 }), s2 * s1);
    EXPECT_EQ(make_serx(1009, { 757, 819, 301, 431 }), s2 / s1);
    EXPECT_EQ(make_serx(1009, { 11, -33, 55, 77 }), s2 * modx(11, 1009));
    EXPECT_EQ(make_serx(1009, { 367, 917, 826, 551 }), s2 / modx(11, 1009));
}

TEST(series_modx_test, operators_inplace) {
    const auto s1 = make_serx(1009, { 4, 1, 0, 0 });
    const auto s2 = make_serx(1009, { 1, -3, 5, 7 });
    const auto s3 = make_serx(1009, { 11, -33, 55, 77 });
    serx sr;
    sr = s1; sr += s2;
    EXPECT_EQ(make_serx(1009, { 5, -2, 5, 7 }), sr);
    sr = s1; sr -= s2;
    EXPECT_EQ(make_serx(1009, { 3, 4, -5, -7 }), sr);
    sr = s1; sr *= s2;
    EXPECT_EQ(make_serx(1009, { 4, -11, 17, 33 }), sr);
    sr = s1; sr /= s2;
    EXPECT_EQ(make_serx(1009, { 4, 13, 19, -36 }), sr);
    sr = s2; sr += s1;
    EXPECT_EQ(make_serx(1009, { 5, -2, 5, 7 }), sr);
    sr = s2; sr -= s1;
    EXPECT_EQ(make_serx(1009, { -3, -4, 5, 7 }), sr);
    sr = s2; sr *= s1;
    EXPECT_EQ(make_serx(1009, { 4, -11, 17, 33 }), sr);
    sr = s2; sr /= s1;
    EXPECT_EQ(make_serx(1009, { 757, 819, 301, 431 }), sr);
    sr = s2; sr *= modx(11, 1009);
    EXPECT_EQ(make_serx(1009, { 11, -33, 55, 77 }), sr);
    sr = s3; sr /= modx(11, 1009);
    EXPECT_EQ(make_serx(1009, { 1, -3, 5, 7 }), sr);
}

TEST(series_modx_test, operators_inplace_self) {
    const auto s1 = make_serx(1009, { 2, -3, 5, 7 });
    serx sr;
    sr = s1; sr += sr;
    EXPECT_EQ(make_serx(1009, { 4, -6, 10, 14 }), sr);
    sr = s1; sr -= sr;
    EXPECT_EQ(make_serx(1009, { 0, 0, 0, 0 }), sr);
    sr = s1; sr *= sr;
    EXPECT_EQ(make_serx(1009, { 4, -12, 29, -2 }), sr);
    sr = s1; sr /= sr;
    EXPECT_EQ(make_serx(1009, { 1 }), sr);
    EXPECT_EQ(4, sr.N());
}

TEST(series_modx_test, shift) {
    const auto s = make_serx(1009, { 7, 5, -3, 4, 2, 1, -8 });
    EXPECT_EQ(make_serx(1009, { 0, 0, 0, 7, 5, -3, 4 }), s.shift(+3));
    EXPECT_EQ(make_serx(1009, { 4, 2, 1, -8, 0, 0, 0 }), s.shift(-3));
}

TEST(series_modx_test, sub_mul) {
    const auto s = make_serx(1009, { 7, 5, -3, 4 });
    EXPECT_EQ(make_serx(1009, { 7, -15, -27, -108 }), s.sub_mul(modx(-3, 1009)));
}

TEST(series_modx_test, sub_pow) {
    const auto s = make_serx(1009, { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 });
    EXPECT_EQ(make_serx(1009, { 1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0, 0, 5 }), s.sub_pow(3));
}

TEST(series_modx_test, derivative) {
    const auto s = make_serx(1009, { 7, 5, -3, 4 });
    const auto sd = s.derivative();
    EXPECT_EQ(make_serx(1009, { 5, -6, 12 }), sd);
    EXPECT_EQ(4, sd.N());
}

TEST(series_modx_test, integral) {
    const auto s = make_serx(1009, { 7, 8, 15, -4, 20 });
    const auto si0 = s.integral();
    const auto si3 = s.integral(3);
    EXPECT_EQ(make_serx(1009, { 0, 7, 4, 5, -1 }), si0);
    EXPECT_EQ(5, si0.N());
    EXPECT_EQ(make_serx(1009, { 3, 7, 4, 5, -1 }), si3);
    EXPECT_EQ(5, si3.N());
}

TEST(series_modx_test, exp) {
    const auto s = make_serx(1009, { 0, 2, 3, 5, 7 });
    EXPECT_EQ(make_serx(1009, { 1, 2, 5, 685, 869 }), s.exp());
}

TEST(series_modx_test, ln) {
    const auto s = make_serx(1009, { 1, -36, 654, -7836, +68673 });
    EXPECT_EQ(make_serx(1009, { 0, -36, 6, 156, 399 }), s.ln());
    EXPECT_EQ(make_serx(1009, { 5, -36, 6, 156, 399 }), s.ln(5));
}

TEST(series_modx_test, pow) {
    const auto s1 = make_serx(1009, { 1, 2, 3, 5, 7, 11, 13, 17, 19, 23 });
    EXPECT_EQ(make_serx(1009, { 1, 6, 21, 59, 144, 321, 663, 1284, 2358, 4133 }), s1.pow(3, 0));
    const auto s2 = make_serx(1009, { 4, 2, 3, 5, 7, 11, 13, 17, 19, 23 });
    EXPECT_EQ(make_serx(1009, { 64, 96, 192, 392, 720, 1338, 2247, 3741, 5958, 9326 }), s2.pow(3, 0));
    const auto s3 = make_serx(1009, { 0, 0, 4, 2, 3, 5, 7, 11, 13, 17 });
    EXPECT_EQ(make_serx(1009, { 0, 0, 0, 0, 0, 0, 64, 96, 192, 392 }), s3.pow(3, 0));
    EXPECT_EQ(make_serx(1009, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }), s3.pow(5, 0));
}

TEST(series_modx_test, static_exp) {
    EXPECT_EQ(make_serx(1009, { 1, 0, 0, 0, 0 }), serx::exp(modx(0, 1009), 5));
    EXPECT_EQ(make_serx(1009, { 1, 1, 505, 841, 967 }), serx::exp(modx(1, 1009), 5));
    EXPECT_EQ(make_serx(1009, { 1, 1008, 505, 168, 967 }), serx::exp(modx(-1, 1009), 5));
    EXPECT_EQ(make_serx(1009, { 1, 2, 2, 674, 337 }), serx::exp(modx(2, 1009), 5));
    EXPECT_EQ(make_serx(1009, { 1, 30, 450, 464, 453 }), serx::exp(modx(30, 1009), 5));
}

TEST(series_modx_test, make_exp_ord) {
    const auto s1 = make_serx(1009, { 2, -3, 5, -9, 12 });
    const auto s2 = s1.make_exponential();
    EXPECT_EQ(make_serx(1009, { 2, 1006, 507, 503, 505 }), s2);
    EXPECT_EQ(make_serx(1009, { 2, 1006, 10, 955, 288 }), s1.make_ordinary());
}

TEST(series_modx_test, of) {
    EXPECT_EQ(make_serx(1009, { 0, 1, 3, 6, 10, 15, 21, 28, 36, 45 }), serx::of([](int n){ return modx(n, 1009) * (n + 1) / 2; }, 10));
}

TEST(series_modx_test, casts) {
    serx s1{ { 2, 1009 }, { 3, 1009 }, { 5, 1009 } };
    EXPECT_EQ(2, s1[0].v);
    EXPECT_EQ(1009, s1[0].M());
    EXPECT_EQ(3, s1[1].v);
    EXPECT_EQ(1009, s1[1].M());
    EXPECT_EQ(5, s1[2].v);
    EXPECT_EQ(1009, s1[2].M());
    EXPECT_EQ(0, s1[3].v);
    EXPECT_EQ(1009, s1[3].M());
    serx e0 = zeroT<serx>::of(s1);
    EXPECT_EQ(0, e0[0].v);
    EXPECT_EQ(1009, e0[0].M());
    EXPECT_EQ(0, e0.p.deg());
    serx e1 = identityT<serx>::of(s1);
    EXPECT_EQ(1, e1[0].v);
    EXPECT_EQ(1009, e1[0].M());
    EXPECT_EQ(0, e1.p.deg());

    const serx s2 = castOf<serx>(s1);
    EXPECT_EQ(2, s2[0].v);
    EXPECT_EQ(1009, s2[0].M());
    EXPECT_EQ(3, s2[1].v);
    EXPECT_EQ(1009, s2[1].M());
    EXPECT_EQ(5, s2[2].v);
    EXPECT_EQ(1009, s2[2].M());
    EXPECT_EQ(0, s2[3].v);
    EXPECT_EQ(1009, s2[3].M());
    const serx s3 = castOf(e1, s1);
    EXPECT_EQ(2, s3[0].v);
    EXPECT_EQ(1009, s3[0].M());
    EXPECT_EQ(3, s3[1].v);
    EXPECT_EQ(1009, s3[1].M());
    EXPECT_EQ(5, s3[2].v);
    EXPECT_EQ(1009, s3[2].M());
    EXPECT_EQ(0, s3[3].v);
    EXPECT_EQ(1009, s3[3].M());
    const serx s4 = castOf(e1, 4);
    EXPECT_EQ(4, s4[0].v);
    EXPECT_EQ(1009, s4[0].M());
    EXPECT_EQ(0, s4[1].v);
    EXPECT_EQ(1009, s4[1].M());
}
