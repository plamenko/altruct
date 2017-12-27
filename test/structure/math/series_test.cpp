#include "altruct/structure/math/series.h"
#include "altruct/structure/math/modulo.h"
#include "structure_test_util.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;
using namespace altruct::test_util;

namespace {
class A {
public:
    double v;
    A(double v) : v(v) {}
    A(int v) : v(v) {}
    bool operator == (const A& rhs) const { return v == rhs.v; }
    bool operator < (const A& rhs) const { return v < rhs.v; }
};
}

TEST(series_test, constructor) {
    const polynom<int> p = { 1, 2, 3, 4 };
    series<int, 4> s1;
    EXPECT_EQ((polynom<int>{}), s1.p);
    series<int, 4> s2(p);
    EXPECT_EQ(p, s2.p);
    series<int, 4> s3(s2);
    EXPECT_EQ(p, s3.p);
    series<int, 4> s4(p.c.begin(), p.c.end());
    EXPECT_EQ(p, s4.p);
    series<int, 4> s5(&p[0], &p[0] + p.size());
    EXPECT_EQ(p, s5.p);
    series<int, 4> s6(5);
    EXPECT_EQ((polynom<int>{ 5 }), s6.p);
    series<int, 4> s7{ 1, 2, 3, 4 };
    EXPECT_EQ(p, s7.p);
    series<A, 4> q1(5);
    EXPECT_EQ((polynom<A>{5.0}), q1.p);
    series<A, 4> q2(A(5.3));
    EXPECT_EQ((polynom<A>{5.3}), q2.p);
    series<int, 4> s8(vector<int>{ 1, 2, 3, 4 });
    EXPECT_EQ(p, s8.p);
    series<int, 4> s9(polynom<int>{ 1, 2, 3, 4 });
    EXPECT_EQ(p, s9.p);
    series<int, 4> s10(series<int, 4>{ 1, 2, 3, 4 });
    EXPECT_EQ(p, s10.p);
    series<int, 4> s11 = series<int, 4>{ 1, 2, 3, 4 };
    EXPECT_EQ(p, s11.p);
    series<int, 4> s12 = s11;
    EXPECT_EQ(p, s12.p);
}

TEST(series_test, swap) {
    series<int, 4> s1{ 1, 2, 3, 4 };
    series<int, 4> s2{ 5, 6, 7 };
    s1.swap(s2);
    EXPECT_EQ((polynom<int>{ 5, 6, 7 }), s1.p);
    EXPECT_EQ((polynom<int>{ 1, 2, 3, 4 }), s2.p);
}

TEST(series_test, size) {
    series<int, 5> s{ 1, 2, 3 };
    EXPECT_EQ(5, s.size());
}

TEST(series_test, at) {
    const series<int, 4> s{ 2, 3, 5, 7 };
    EXPECT_EQ(2, s.at(0));
    EXPECT_EQ(3, s.at(1));
    EXPECT_EQ(5, s.at(2));
    EXPECT_EQ(7, s.at(3));
    EXPECT_EQ(0, s.at(4));
    EXPECT_EQ(0, s.at(100));
    EXPECT_EQ(4, s.size());
}

TEST(series_test, operator_const_brackets) {
    const series<int, 4> s{ 2, 3, 5, 7 };
    EXPECT_EQ(2, s[0]);
    EXPECT_EQ(7, s[3]);
    EXPECT_EQ(0, s[4]);
    EXPECT_EQ(0, s[100]);
    EXPECT_EQ(4, s.size());
}

TEST(series_test, operator_brackets) {
    series<int, 4> s;
    s[3] = 3;
    EXPECT_EQ(0, s[0]);
    EXPECT_EQ(0, s[4]);
    EXPECT_EQ(3, s[3]);
    EXPECT_EQ(0, s[4]);
    EXPECT_EQ(0, s[100]);
    EXPECT_EQ(4, s.size());
}

TEST(series_test, operators_comparison) {
    const series<int, 4> s1{ 4 };
    const series<int, 4> s2{ 1, 3, 5, 7 };
    const series<int, 4> s3{ 1, 3, 5, 7, 0, 0, 0 };
    ASSERT_COMPARISON_OPERATORS(0, s1, s1);
    ASSERT_COMPARISON_OPERATORS(0, s2, s2);
    ASSERT_COMPARISON_OPERATORS(0, s3, s3);
    ASSERT_COMPARISON_OPERATORS(-1, s1, s2);
    ASSERT_COMPARISON_OPERATORS(+1, s2, s1);
    ASSERT_COMPARISON_OPERATORS(0, s2, s3);
    ASSERT_COMPARISON_OPERATORS(0, s3, s2);
}

TEST(series_test, inverse) {
    const series<int, 4> s{ 1, -3, 5, 7 };
    series<int, 4> si = s.inverse();
    EXPECT_EQ((polynom<int>{ 1, 3, 4, -10 }), si.p);
}

TEST(series_test, operators_arithmetic) {
    const series<double, 4> s1{ 4, 1 };
    const series<double, 4> s2{ 1, -3, 5, 7 };
    EXPECT_EQ((series<double, 4>{ 5, -2, 5, 7 }), s1 + s2);
    EXPECT_EQ((series<double, 4>{ 3, 4, -5, -7 }), s1 - s2);
    EXPECT_EQ((series<double, 4>{ -1, 3, -5, -7 }), -s2);
    EXPECT_EQ((series<double, 4>{ 4, -11, 17, 33 }), s1 * s2);
    EXPECT_EQ((series<double, 4>{ 4, 13, 19, -36 }), s1 / s2);
    EXPECT_EQ((series<double, 4>{ 5, -2, 5, 7 }), s2 + s1);
    EXPECT_EQ((series<double, 4>{ -3, -4, 5, 7 }), s2 - s1);
    EXPECT_EQ((series<double, 4>{ -4, -1 }), -s1);
    EXPECT_EQ((series<double, 4>{ 4, -11, 17, 33 }), s2 * s1);
    EXPECT_EQ((series<double, 4>{ 1 / 4., -13 / 16., 93 / 64., 355 / 256. }), s2 / s1);
    EXPECT_EQ((series<double, 4>{ 11, -33, 55, 77 }), s2 * 11);
}

TEST(series_test, operators_inplace) {
    const series<double, 4> s1{ 4, 1 };
    const series<double, 4> s2{ 1, -3, 5, 7 };
    const series<double, 4> s3{ 11, -33, 55, 77 };
    series<double, 4> sr;
    sr = s1; sr += s2;
    EXPECT_EQ((series<double, 4>{ 5, -2, 5, 7 }), sr);
    sr = s1; sr -= s2;
    EXPECT_EQ((series<double, 4>{ 3, 4, -5, -7 }), sr);
    sr = s1; sr *= s2;
    EXPECT_EQ((series<double, 4>{ 4, -11, 17, 33 }), sr);
    sr = s1; sr /= s2;
    EXPECT_EQ((series<double, 4>{ 4, 13, 19, -36 }), sr);
    sr = s2; sr += s1;
    EXPECT_EQ((series<double, 4>{ 5, -2, 5, 7 }), sr);
    sr = s2; sr -= s1;
    EXPECT_EQ((series<double, 4>{ -3, -4, 5, 7 }), sr);
    sr = s2; sr *= s1;
    EXPECT_EQ((series<double, 4>{ 4, -11, 17, 33 }), sr);
    sr = s2; sr /= s1;
    EXPECT_EQ((series<double, 4>{ 1 / 4., -13 / 16., 93 / 64., 355 / 256. }), sr);
    sr = s2; sr *= 11;
    EXPECT_EQ((series<double, 4>{ 11, -33, 55, 77 }), sr);
    sr = s3; sr /= 11;
    EXPECT_EQ((series<double, 4>{ 1, -3, 5, 7 }), sr);
}

TEST(series_test, operators_inplace_self) {
    const series<double, 4> s1{ 2, -3, 5, 7 };
    series<double, 4> sr;
    sr = s1; sr += sr;
    EXPECT_EQ((series<double, 4>{ 4, -6, 10, 14 }), sr);
    sr = s1; sr -= sr;
    EXPECT_EQ((series<double, 4>{ 0, 0, 0, 0 }), sr);
    sr = s1; sr *= sr;
    EXPECT_EQ((series<double, 4>{ 4, -12, 29, -2 }), sr);
    sr = s1; sr /= sr;
    EXPECT_EQ((series<double, 4>{ 1 }), sr);
}

TEST(series_test, shift) {
    const series<int, 7> s{ 7, 5, -3, 4, 2, 1, -8 };
    EXPECT_EQ((series<int, 7>{ 0, 0, 0, 7, 5, -3, 4 }), s.shift(+3));
    EXPECT_EQ((series<int, 7>{ 4, 2, 1, -8, 0, 0, 0 }), s.shift(-3));
}

TEST(series_test, sub_mul) {
    const series<int, 4> s{ 7, 5, -3, 4 };
    EXPECT_EQ((series<int, 4>{ 7, -15, -27, -108 }), s.sub_mul(-3));
}

TEST(series_test, sub_pow) {
    const series<int, 13> s{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 };
    EXPECT_EQ((series<int, 13>{ 1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0, 0, 5 }), s.sub_pow(3));
}

TEST(series_test, derivative) {
    const series<int, 4> s{ 7, 5, -3, 4 };
    EXPECT_EQ((series<int, 4>{ 5, -6, 12 }), s.derivative());
}

TEST(series_test, integral) {
    const series<int, 5> s{ 7, 8, 15, -4, 20 };
    EXPECT_EQ((series<int, 5>{ 0, 7, 4, 5, -1 }), s.integral());
    EXPECT_EQ((series<int, 5>{ 3, 7, 4, 5, -1 }), s.integral(3));
}

TEST(series_test, exp) {
    const series<double, 5> s{ 0, 2, 3, 5, 7 };
    EXPECT_EQ((series<double, 5>{ 1, 2, 5, 12 + 1/3., 28 + 1/6.}), s.exp());
}

TEST(series_test, ln) {
    const series<double, 5> s1{ 1, -36, 654, -7836, +68673 };
    EXPECT_EQ((series<double, 5>{ 0, -36, 6, 156, 399 }), s1.ln());
    EXPECT_EQ((series<double, 5>{ 5, -36, 6, 156, 399 }), s1.ln(5));
}

TEST(series_test, pow) {
    const series<double, 10> s1{ 1, 2, 3, 5, 7, 11, 13, 17, 19, 23 };
    EXPECT_EQ((series<double, 10>{1, 6, 21, 59, 144, 321, 663, 1284, 2358, 4133}), s1.pow(3, 0));
    const series<double, 10> s2{ 4, 2, 3, 5, 7, 11, 13, 17, 19, 23 };
    EXPECT_EQ((series<double, 10>{64, 96, 192, 392, 720, 1338, 2247, 3741, 5958, 9326}), s2.pow(3, 0));
    const series<double, 10> s3{ 0, 0, 4, 2, 3, 5, 7, 11, 13, 17 };
    EXPECT_EQ((series<double, 10>{0, 0, 0, 0, 0, 0, 64, 96, 192, 392}), s3.pow(3, 0));
    EXPECT_EQ((series<double, 10>{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}), s3.pow(5, 0));
}

TEST(series_test, static_exp) {
    EXPECT_EQ((series<double, 5>{ 1, 0, 0, 0, 0 }), (series<double, 5>::exp(0)));
    EXPECT_EQ((series<double, 5>{ 1 / 1., 1 / 1., 1 / 2., 1 / 6., 1 / 24. }), (series<double, 5>::exp(1)));
    EXPECT_EQ((series<double, 5>{ 1 / 1., -1 / 1., 1 / 2., -1 / 6., 1 / 24. }), (series<double, 5>::exp(-1)));
    EXPECT_EQ((series<double, 5>{ 1 / 1., 2 / 1., 4 / 2., 8 / 6., 16 / 24. }), (series<double, 5>::exp(2)));
    EXPECT_EQ((series<double, 5>{ 1 / 1., 30 / 1., 900 / 2., 27000 / 6., 810000 / 24. }), (series<double, 5>::exp(30)));
}

TEST(series_test, make_exp_ord) {
    const series<double, 5> s1{ 2, -3, 5, -9, 12 };
    auto s2 = s1.make_exponential();
    EXPECT_EQ((series<double, 5>{ 2 / 1., -3 / 1., 5 / 2., -9 / 6., 12 / 24. }), s2);
    EXPECT_EQ((series<double, 5>{ 2 * 1., -3 * 1., 5 * 2., -9 * 6., 12 * 24. }), s1.make_ordinary());
}

TEST(series_test, of) {
    EXPECT_EQ((series<int, 10>{ 0, 1, 3, 6, 10, 15, 21, 28, 36, 45 }), (series<int, 10>::of([](int n){ return n * (n + 1) / 2; })));
}

TEST(series_test, casts) {
    typedef modulo<int, 1009> mod;
    typedef polynom<mod> poly;
    typedef series<mod, 4> ser;
    ser s1{ { 2 }, { 3 }, { 5 } };
    EXPECT_EQ(2, s1[0].v);
    EXPECT_EQ(1009, s1[0].M());
    EXPECT_EQ(3, s1[1].v);
    EXPECT_EQ(1009, s1[1].M());
    EXPECT_EQ(5, s1[2].v);
    EXPECT_EQ(1009, s1[2].M());
    EXPECT_EQ(0, s1[3].v);
    EXPECT_EQ(1009, s1[3].M());
    ser e0 = zeroT<ser>::of(s1);
    EXPECT_EQ(0, e0[0].v);
    EXPECT_EQ(1009, e0[0].M());
    EXPECT_EQ(0, e0.p.deg());
    ser e1 = identityT<ser>::of(s1);
    EXPECT_EQ(1, e1[0].v);
    EXPECT_EQ(1009, e1[0].M());
    EXPECT_EQ(0, e1.p.deg());

    const ser s2 = castOf<ser>(s1);
    EXPECT_EQ(2, s2[0].v);
    EXPECT_EQ(1009, s2[0].M());
    EXPECT_EQ(3, s2[1].v);
    EXPECT_EQ(1009, s2[1].M());
    EXPECT_EQ(5, s2[2].v);
    EXPECT_EQ(1009, s2[2].M());
    EXPECT_EQ(0, s2[3].v);
    EXPECT_EQ(1009, s2[3].M());
    const ser s3 = castOf(e1, s1);
    EXPECT_EQ(2, s3[0].v);
    EXPECT_EQ(1009, s3[0].M());
    EXPECT_EQ(3, s3[1].v);
    EXPECT_EQ(1009, s3[1].M());
    EXPECT_EQ(5, s3[2].v);
    EXPECT_EQ(1009, s3[2].M());
    EXPECT_EQ(0, s3[3].v);
    EXPECT_EQ(1009, s3[3].M());
    const ser s4 = castOf(e1, 4);
    EXPECT_EQ(4, s4[0].v);
    EXPECT_EQ(1009, s4[0].M());
    EXPECT_EQ(0, s4[1].v);
    EXPECT_EQ(1009, s4[1].M());
    const ser s5 = castOf<ser>(5);
    EXPECT_EQ(5, s5[0].v);
    EXPECT_EQ(1009, s5[0].M());
    EXPECT_EQ(0, s5[1].v);
    EXPECT_EQ(1009, s5[1].M());
}
