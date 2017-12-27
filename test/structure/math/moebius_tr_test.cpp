#include "altruct/structure/math/moebius_tr.h"
#include "altruct/structure/math/with_infinity.h"
#include "altruct/structure/math/complex.h"
#include "altruct/structure/math/modulo.h"
#include "structure_test_util.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;
using namespace altruct::test_util;

typedef modulo<int, 1009, modulo_storage::CONSTANT> mod;
typedef altruct::math::complex<mod> cplx;
typedef with_infinity<cplx> winf;
typedef moebius_tr<winf> moeb_tr;

template<typename W>
vector<int> to_vec_w(const W& w) {
    return{ w.v.a.v, w.v.a.M(), w.v.b.v, w.v.b.M(), w.v.D().v, w.v.D().M(), w.is_inf };
}

template<typename T>
vector<vector<int>> to_vec(const T& t) {
    return{ to_vec_w(t.a), to_vec_w(t.b), to_vec_w(t.c), to_vec_w(t.d), { t.s } };
}

TEST(moebius_tr_test, constructor) {
    moeb_tr t1;
    EXPECT_EQ((vector<vector<int>> {{ 0, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 1, 1009, 0, 1009, 1008, 1009, 0 }, { 1 }}), to_vec(t1));
    moeb_tr t2(10);
    EXPECT_EQ((vector<vector<int>> {{ 10, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 1, 1009, 0, 1009, 1008, 1009, 0 }, { 1 }}), to_vec(t2));
    moeb_tr t3(cplx(+2, -5));
    EXPECT_EQ((vector<vector<int>> {{ 2, 1009, 1004, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 1, 1009, 0, 1009, 1008, 1009, 0 }, { 1 }}), to_vec(t3));
    moeb_tr t4(cplx(+2, -5), cplx(3, 7), cplx(-1, 4), cplx(8, 6));
    EXPECT_EQ((vector<vector<int>> {{ 2, 1009, 1004, 1009, 1008, 1009, 0 }, { 3, 1009, 7, 1009, 1008, 1009, 0 }, { 1008, 1009, 4, 1009, 1008, 1009, 0 }, { 8, 1009, 6, 1009, 1008, 1009, 0 }, { 1 }}), to_vec(t4));
    moeb_tr t5(cplx(+2, -5), cplx(3, 7), cplx(-1, 4), cplx(8, 6), -1);
    EXPECT_EQ((vector<vector<int>> {{ 2, 1009, 1004, 1009, 1008, 1009, 0 }, { 3, 1009, 7, 1009, 1008, 1009, 0 }, { 1008, 1009, 4, 1009, 1008, 1009, 0 }, { 8, 1009, 6, 1009, 1008, 1009, 0 }, { -1 }}), to_vec(t5));
    moeb_tr t6(t5);
    EXPECT_EQ((vector<vector<int>> {{ 2, 1009, 1004, 1009, 1008, 1009, 0 }, { 3, 1009, 7, 1009, 1008, 1009, 0 }, { 1008, 1009, 4, 1009, 1008, 1009, 0 }, { 8, 1009, 6, 1009, 1008, 1009, 0 }, { -1 }}), to_vec(t6));
}

TEST(moebius_tr_test, operators_comparison) {
    moeb_tr t0(cplx(+1, -5), cplx(3, 7), cplx(-1, 4), cplx(8, 6), +1);
    moeb_tr t1(cplx(+1, -5), cplx(3, 7), cplx(-1, 4), cplx(8, 6), +1);
    moeb_tr t2(cplx(+2, -5), cplx(3, 7), cplx(-1, 4), cplx(8, 6), +1);
    moeb_tr t3(cplx(+1, -5), cplx(4, 7), cplx(-1, 4), cplx(8, 6), +1);
    moeb_tr t4(cplx(+1, -5), cplx(3, 7), cplx(-1, 5), cplx(8, 6), +1);
    moeb_tr t5(cplx(+1, -5), cplx(3, 7), cplx(-1, 4), cplx(8, 7), +1);
    moeb_tr t6(cplx(+1, -5), cplx(3, 7), cplx(-1, 4), cplx(8, 6), -1);
    ASSERT_COMPARISON_OPERATORS(0, t0, t1);
    ASSERT_COMPARISON_OPERATORS(-1, t0, t2);
    ASSERT_COMPARISON_OPERATORS(-1, t0, t3);
    ASSERT_COMPARISON_OPERATORS(-1, t0, t4);
    ASSERT_COMPARISON_OPERATORS(-1, t0, t5);
    ASSERT_COMPARISON_OPERATORS(+1, t0, t6);
}

TEST(moebius_tr_test, operators_arithmetic) {
    moeb_tr t1(cplx(+1, -5), cplx(3, 7), cplx(-1, 4), cplx(8, 6), +1);
    moeb_tr t2(cplx(-2, +3), cplx(-4, 5), cplx(6, -7), cplx(-8, -9), +1);
    moeb_tr t3(cplx(+1, -5), cplx(3, 7), cplx(-1, 4), cplx(8, 6), -1);
    moeb_tr t4(cplx(-2, +3), cplx(-4, 5), cplx(6, -7), cplx(-8, -9), -1);
    EXPECT_EQ((vector<vector<int>>{{ 80, 1009, 34, 1009, 1008, 1009, 0 }, { 60, 1009, 951, 1009, 1008, 1009, 0 }, { 80, 1009, 978, 1009, 1008, 1009, 0 }, { 983, 1009, 868, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(t1 * t2));
    EXPECT_EQ((vector<vector<int>>{{ 80, 1009, 34, 1009, 1008, 1009, 0 }, { 60, 1009, 951, 1009, 1008, 1009, 0 }, { 80, 1009, 978, 1009, 1008, 1009, 0 }, { 983, 1009, 868, 1009, 1008, 1009, 0 }, { -1 }}), to_vec(t1 * t4));
    EXPECT_EQ((vector<vector<int>>{{ 961, 1009, 70, 1009, 1008, 1009, 0 }, { 893, 1009, 995, 1009, 1008, 1009, 0 }, { 20, 1009, 87, 1009, 1008, 1009, 0 }, { 915, 1009, 13, 1009, 1008, 1009, 0 }, { -1 }}), to_vec(t3 * t2));
    EXPECT_EQ((vector<vector<int>>{{ 961, 1009, 70, 1009, 1008, 1009, 0 }, { 893, 1009, 995, 1009, 1008, 1009, 0 }, { 20, 1009, 87, 1009, 1008, 1009, 0 }, { 915, 1009, 13, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(t3 * t4));
    EXPECT_EQ((vector<vector<int>>{{ 882, 1009, 440, 1009, 1008, 1009, 0 }, { 328, 1009, 592, 1009, 1008, 1009, 0 }, { 252, 1009, 94, 1009, 1008, 1009, 0 }, { 731, 1009, 927, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(t1 / t2));
    EXPECT_EQ((vector<vector<int>>{{ 606, 1009, 644, 1009, 1008, 1009, 0 }, { 253, 1009, 694, 1009, 1008, 1009, 0 }, { 960, 1009, 560, 1009, 1008, 1009, 0 }, { 429, 1009, 132, 1009, 1008, 1009, 0 }, { -1 }}), to_vec(t1 / t4));
    EXPECT_EQ((vector<vector<int>>{{ 757, 1009, 945, 1009, 1008, 1009, 0 }, { 807, 1009, 592, 1009, 1008, 1009, 0 }, { 881, 1009, 914, 1009, 1008, 1009, 0 }, { 50, 1009, 183, 1009, 1008, 1009, 0 }, { -1 }}), to_vec(t3 / t2));
    EXPECT_EQ((vector<vector<int>> {{ 630, 1009, 443, 1009, 1008, 1009, 0 }, { 75, 1009, 795, 1009, 1008, 1009, 0 }, { 933, 1009, 751, 1009, 1008, 1009, 0 }, { 303, 1009, 877, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(t3 / t4));
}

TEST(moebius_tr_test, operators_inplace) {
    moeb_tr t1(cplx(+1, -5), cplx(3, 7), cplx(-1, 4), cplx(8, 6), +1);
    moeb_tr t2(cplx(-2, +3), cplx(-4, 5), cplx(6, -7), cplx(-8, -9), +1);
    moeb_tr t3(cplx(+1, -5), cplx(3, 7), cplx(-1, 4), cplx(8, 6), -1);
    moeb_tr t4(cplx(-2, +3), cplx(-4, 5), cplx(6, -7), cplx(-8, -9), -1);
    moeb_tr tr;
    tr = t1; tr *= t2;
    EXPECT_EQ((vector<vector<int>>{{ 80, 1009, 34, 1009, 1008, 1009, 0 }, { 60, 1009, 951, 1009, 1008, 1009, 0 }, { 80, 1009, 978, 1009, 1008, 1009, 0 }, { 983, 1009, 868, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(tr));
    tr = t1; tr *= t4;
    EXPECT_EQ((vector<vector<int>>{{ 80, 1009, 34, 1009, 1008, 1009, 0 }, { 60, 1009, 951, 1009, 1008, 1009, 0 }, { 80, 1009, 978, 1009, 1008, 1009, 0 }, { 983, 1009, 868, 1009, 1008, 1009, 0 }, { -1 }}), to_vec(tr));
    tr = t3; tr *= t2;
    EXPECT_EQ((vector<vector<int>>{{ 961, 1009, 70, 1009, 1008, 1009, 0 }, { 893, 1009, 995, 1009, 1008, 1009, 0 }, { 20, 1009, 87, 1009, 1008, 1009, 0 }, { 915, 1009, 13, 1009, 1008, 1009, 0 }, { -1 }}), to_vec(tr));
    tr = t3; tr *= t4;
    EXPECT_EQ((vector<vector<int>>{{ 961, 1009, 70, 1009, 1008, 1009, 0 }, { 893, 1009, 995, 1009, 1008, 1009, 0 }, { 20, 1009, 87, 1009, 1008, 1009, 0 }, { 915, 1009, 13, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(tr));
    tr = t1; tr /= t2;
    EXPECT_EQ((vector<vector<int>>{{ 882, 1009, 440, 1009, 1008, 1009, 0 }, { 328, 1009, 592, 1009, 1008, 1009, 0 }, { 252, 1009, 94, 1009, 1008, 1009, 0 }, { 731, 1009, 927, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(tr));
    tr = t1; tr /= t4;
    EXPECT_EQ((vector<vector<int>>{{ 606, 1009, 644, 1009, 1008, 1009, 0 }, { 253, 1009, 694, 1009, 1008, 1009, 0 }, { 960, 1009, 560, 1009, 1008, 1009, 0 }, { 429, 1009, 132, 1009, 1008, 1009, 0 }, { -1 }}), to_vec(tr));
    tr = t3; tr /= t2;
    EXPECT_EQ((vector<vector<int>>{{ 757, 1009, 945, 1009, 1008, 1009, 0 }, { 807, 1009, 592, 1009, 1008, 1009, 0 }, { 881, 1009, 914, 1009, 1008, 1009, 0 }, { 50, 1009, 183, 1009, 1008, 1009, 0 }, { -1 }}), to_vec(tr));
    tr = t3; tr /= t4;
    EXPECT_EQ((vector<vector<int>> {{ 630, 1009, 443, 1009, 1008, 1009, 0 }, { 75, 1009, 795, 1009, 1008, 1009, 0 }, { 933, 1009, 751, 1009, 1008, 1009, 0 }, { 303, 1009, 877, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(tr));
}

TEST(moebius_tr_test, operators_inplace_self) {
    moeb_tr t1(cplx(+1, -5), cplx(3, 7), cplx(-1, 4), cplx(8, 6), +1);
    moeb_tr t3(cplx(+1, -5), cplx(3, 7), cplx(-1, 4), cplx(8, 6), -1);
    moeb_tr tr;
    tr = t1; tr *= tr;
    EXPECT_EQ((vector<vector<int>>{{ 954, 1009, 1004, 1009, 1008, 1009, 0 }, { 20, 1009, 66, 1009, 1008, 1009, 0 }, { 996, 1009, 35, 1009, 1008, 1009, 0 }, { 1006, 1009, 101, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(tr));
    tr = t3; tr *= tr;
    EXPECT_EQ((vector<vector<int>>{{ 51, 1009, 990, 1009, 1008, 1009, 0 }, { 34, 1009, 16, 1009, 1008, 1009, 0 }, { 1004, 1009, 970, 1009, 1008, 1009, 0 }, { 125, 1009, 19, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(tr));
    tr = t1; tr /= tr;
    EXPECT_EQ((vector<vector<int>>{{ 1, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 1, 1009, 0, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(tr));
    tr = t3; tr /= tr; tr.normalize();
    EXPECT_EQ((vector<vector<int>>{{ 1, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 1, 1009, 0, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(tr));
}

TEST(moebius_tr_test, inverse) {
    moeb_tr t1(cplx(+1, -5), cplx(3, 7), cplx(-1, 4), cplx(8, 6), +1);
    moeb_tr t3(cplx(+1, -5), cplx(3, 7), cplx(-1, 4), cplx(8, 6), -1);
    EXPECT_EQ((vector<vector<int>>{{ 559, 1009, 667, 1009, 1008, 1009, 0 }, { 611, 1009, 316, 1009, 1008, 1009, 0 }, { 386, 1009, 145, 1009, 1008, 1009, 0 }, { 426, 1009, 928, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(t1.inverse()));
    EXPECT_EQ((vector<vector<int>> {{ 30, 1009, 982, 1009, 1008, 1009, 0 }, { 51, 1009, 453, 1009, 1008, 1009, 0 }, { 66, 1009, 944, 1009, 1008, 1009, 0 }, { 26, 1009, 161, 1009, 1008, 1009, 0 }, { -1 }}), to_vec(t3.inverse()));
}

TEST(moebius_tr_test, normalize) {
    moeb_tr t1(cplx(+1, -5), cplx(3, 7), cplx(-1, 4), cplx(8, 6), +1);
    moeb_tr t3(cplx(+1, -5), cplx(3, 7), cplx(-1, 4), cplx(8, 6), -1);
    EXPECT_EQ((vector<vector<int>>{{ 1, 1009, 0, 1009, 1008, 1009, 0 }, { 154, 1009, 777, 1009, 1008, 1009, 0 }, { 38, 1009, 194, 1009, 1008, 1009, 0 }, { 232, 1009, 157, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(t1.normalize()));
    EXPECT_EQ((vector<vector<int>>{{ 1, 1009, 0, 1009, 1008, 1009, 0 }, { 154, 1009, 777, 1009, 1008, 1009, 0 }, { 38, 1009, 194, 1009, 1008, 1009, 0 }, { 232, 1009, 157, 1009, 1008, 1009, 0 }, { -1 }}), to_vec(t3.normalize()));
    moeb_tr t5(cplx(0), cplx(3, 7), cplx(-1, 4), cplx(8, 6), +1);
    moeb_tr t7(cplx(0), cplx(3, 7), cplx(-1, 4), cplx(8, 6), -1);
    EXPECT_EQ((vector<vector<int>>{{ 0, 1009, 0, 1009, 1008, 1009, 0 }, { 1, 1009, 0, 1009, 1008, 1009, 0 }, { 157, 1009, 644, 1009, 1008, 1009, 0 }, { 697, 1009, 730, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(t5.normalize()));
    EXPECT_EQ((vector<vector<int>>{{ 0, 1009, 0, 1009, 1008, 1009, 0 }, { 1, 1009, 0, 1009, 1008, 1009, 0 }, { 157, 1009, 644, 1009, 1008, 1009, 0 }, { 697, 1009, 730, 1009, 1008, 1009, 0 }, { -1 }}), to_vec(t7.normalize()));
    moeb_tr t9(cplx(0), cplx(0), cplx(-1, 4), cplx(8, 6), +1);
    moeb_tr t11(cplx(0), cplx(0), cplx(-1, 4), cplx(8, 6), -1);
    EXPECT_EQ((vector<vector<int>>{{ 0, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 1, 1009, 0, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(t9.normalize()));
    EXPECT_EQ((vector<vector<int>>{{ 0, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 1, 1009, 0, 1009, 1008, 1009, 0 }, { -1 }}), to_vec(t11.normalize()));
}

TEST(moebius_tr_test, apply) {
    moeb_tr t1(cplx(+1, -5), cplx(3, 7), cplx(-1, 4), cplx(8, 6), +1);
    moeb_tr t2(cplx(-2, +3), cplx(-4, 5), cplx(6, -7), cplx(-8, -9), +1);
    moeb_tr t3(cplx(+1, -5), cplx(3, 7), cplx(-1, 4), cplx(8, 6), -1);
    moeb_tr t4(cplx(-2, +3), cplx(-4, 5), cplx(6, -7), cplx(-8, -9), -1);
    moeb_tr t5(cplx(0), cplx(3, 7), cplx(-1, 4), cplx(8, 6), +1);
    moeb_tr t7(cplx(0), cplx(3, 7), cplx(-1, 4), cplx(8, 6), -1);
    moeb_tr t9(cplx(0), cplx(0), cplx(-1, 4), cplx(8, 6), +1);
    moeb_tr t11(cplx(0), cplx(0), cplx(-1, 4), cplx(8, 6), -1);
    winf w(cplx(11, 23));
    winf wi(cplx(0, 0), 1);
    winf ws(-cplx(-8, -9) / cplx(6, -7));
    
    EXPECT_EQ((vector<int>{ 976, 1009, 639, 1009, 1008, 1009, 0 }), to_vec_w(t2(w)));
    EXPECT_EQ((vector<int>{ 214, 1009, 4, 1009, 1008, 1009, 0 }), to_vec_w(t4(w)));

    // anti-moebius(z) == moebius(z*)
    EXPECT_EQ((vector<int>{ 214, 1009, 4, 1009, 1008, 1009, 0 }), to_vec_w(t2(conjugateT<winf>::of(w))));
    EXPECT_EQ((vector<int>{ 976, 1009, 639, 1009, 1008, 1009, 0 }), to_vec_w(t4(conjugateT<winf>::of(w))));

    // singularity point maps to infinity
    EXPECT_EQ((vector<int>{ 586, 1009, 358, 1009, 1008, 1009, 1}), to_vec_w(t2(ws)));
    EXPECT_EQ((vector<int>{ 586, 1009, 358, 1009, 1008, 1009, 1 }), to_vec_w(t4(conjugateT<winf>::of(ws))));

    // product of transformations equals their composition
    EXPECT_EQ((vector<int>{ 163, 1009, 52, 1009, 1008, 1009, 0 }), to_vec_w(t1(t2(w))));
    EXPECT_EQ((vector<int>{ 163, 1009, 52, 1009, 1008, 1009, 0 }), to_vec_w((t1*t2)(w)));
    EXPECT_EQ((vector<int>{ 365, 1009, 311, 1009, 1008, 1009, 0 }), to_vec_w(t1(t4(w))));
    EXPECT_EQ((vector<int>{ 365, 1009, 311, 1009, 1008, 1009, 0 }), to_vec_w((t1*t4)(w)));
    EXPECT_EQ((vector<int>{ 718, 1009, 203, 1009, 1008, 1009, 0 }), to_vec_w(t3(t2(w))));
    EXPECT_EQ((vector<int>{ 718, 1009, 203, 1009, 1008, 1009, 0 }), to_vec_w((t3*t2)(w)));
    EXPECT_EQ((vector<int>{ 227, 1009, 557, 1009, 1008, 1009, 0 }), to_vec_w(t3(t4(w))));
    EXPECT_EQ((vector<int>{ 227, 1009, 557, 1009, 1008, 1009, 0 }), to_vec_w((t3*t4)(w)));

    // inverse
    EXPECT_EQ((vector<int>{ 11, 1009, 23, 1009, 1008, 1009, 0 }), to_vec_w(t2.inverse()(t2(w))));
    EXPECT_EQ((vector<int>{ 11, 1009, 23, 1009, 1008, 1009, 0 }), to_vec_w((t4.inverse()*t4)(w)));

    // normalize
    EXPECT_EQ((vector<int>{ 561, 1009, 712, 1009, 1008, 1009, 0 }), to_vec_w(t1(w)));
    EXPECT_EQ((vector<int>{ 561, 1009, 712, 1009, 1008, 1009, 0 }), to_vec_w(t1.normalize()(w)));
    EXPECT_EQ((vector<int>{ 30, 1009, 190, 1009, 1008, 1009, 0 }), to_vec_w(t3(w)));
    EXPECT_EQ((vector<int>{ 30, 1009, 190, 1009, 1008, 1009, 0 }), to_vec_w(t3.normalize()(w)));
    EXPECT_EQ((vector<int>{ 721, 1009, 789, 1009, 1008, 1009, 0 }), to_vec_w(t5(w)));
    EXPECT_EQ((vector<int>{ 721, 1009, 789, 1009, 1008, 1009, 0 }), to_vec_w(t5.normalize()(w)));
    EXPECT_EQ((vector<int>{ 370, 1009, 921, 1009, 1008, 1009, 0 }), to_vec_w(t7(w)));
    EXPECT_EQ((vector<int>{ 370, 1009, 921, 1009, 1008, 1009, 0 }), to_vec_w(t7.normalize()(w)));
    EXPECT_EQ((vector<int>{ 0, 1009, 0, 1009, 1008, 1009, 0 }), to_vec_w(t9(w)));
    EXPECT_EQ((vector<int>{ 0, 1009, 0, 1009, 1008, 1009, 0 }), to_vec_w(t9.normalize()(w)));
    EXPECT_EQ((vector<int>{ 0, 1009, 0, 1009, 1008, 1009, 0 }), to_vec_w(t11(w)));
    EXPECT_EQ((vector<int>{ 0, 1009, 0, 1009, 1008, 1009, 0 }), to_vec_w(t11.normalize()(w)));

    // at infinity
    EXPECT_EQ((vector<int>{ 973, 1009, 463, 1009, 1008, 1009, 0 }), to_vec_w(t2(wi)));
    EXPECT_EQ((vector<int>{ 973, 1009, 463, 1009, 1008, 1009, 0 }), to_vec_w(t4(wi)));
}

TEST(moebius_tr_test, builtin_transformations) {
    EXPECT_EQ((vector<int>{ 15, 1009, 999, 1009, 1008, 1009, 0 }), to_vec_w(moeb_tr::scaling(cplx(5))(cplx(3, -2))));
    EXPECT_EQ((vector<int>{ 8, 1009, 1007, 1009, 1008, 1009, 0 }), to_vec_w(moeb_tr::translation(cplx(5, -6))(cplx(3, 4))));
    EXPECT_EQ((vector<int>{ 404, 1009, 204, 1009, 1008, 1009, 0 }), to_vec_w(moeb_tr::rotation(cplx(3, 4) / cplx(5))(cplx(2, 1))));
    EXPECT_EQ((vector<int>{ 1007, 1009, 3, 1009, 1008, 1009, 0 }), to_vec_w(moeb_tr::flip_x()(cplx(2, 3))));
    EXPECT_EQ((vector<int>{ 2, 1009, 1006, 1009, 1008, 1009, 0 }), to_vec_w(moeb_tr::flip_y()(cplx(2, 3))));
    EXPECT_EQ((vector<int>{ 6, 1009, 8, 1009, 1008, 1009, 0 }), to_vec_w(moeb_tr::inversion()(cplx(3, 4) / cplx(50))));
}

TEST(moebius_tr_test, casts) {
    moeb_tr t0;
    moeb_tr t1(cplx(+1, -5), cplx(3, 7), cplx(-1, 4), cplx(8, 6), +1);
    moeb_tr t3(cplx(+1, -5), cplx(3, 7), cplx(-1, 4), cplx(8, 6), -1);
    EXPECT_EQ((vector<vector<int>>{{ 1, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 1, 1009, 0, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(identityOf(t1)));
    EXPECT_EQ((vector<vector<int>>{{ 1, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 1, 1009, 0, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(identityOf(t3)));
    EXPECT_EQ((vector<vector<int>>{{ 0, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 1, 1009, 0, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(zeroOf(t1)));
    EXPECT_EQ((vector<vector<int>>{{ 0, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 1, 1009, 0, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(zeroOf(t3)));
    EXPECT_EQ((vector<vector<int>>{{ 3, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 1, 1009, 0, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(castOf<moeb_tr>(3)));
    EXPECT_EQ((vector<vector<int>>{{ 4, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 1, 1009, 0, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(castOf<moeb_tr>(t1, 4)));
    EXPECT_EQ((vector<vector<int>>{{ 4, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 0, 1009, 0, 1009, 1008, 1009, 0 }, { 1, 1009, 0, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(castOf<moeb_tr>(t3, 4)));
    EXPECT_EQ((vector<vector<int>>{{ 1, 1009, 1004, 1009, 1008, 1009, 0 }, { 3, 1009, 7, 1009, 1008, 1009, 0 }, { 1008, 1009, 4, 1009, 1008, 1009, 0 }, { 8, 1009, 6, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(castOf<moeb_tr>(t0, t1)));
    EXPECT_EQ((vector<vector<int>>{{ 1, 1009, 1004, 1009, 1008, 1009, 0 }, { 3, 1009, 7, 1009, 1008, 1009, 0 }, { 1008, 1009, 4, 1009, 1008, 1009, 0 }, { 8, 1009, 6, 1009, 1008, 1009, 0 }, { -1 }}), to_vec(castOf<moeb_tr>(t0, t3)));
    EXPECT_EQ((vector<vector<int>>{{ 1, 1009, 1004, 1009, 1008, 1009, 0 }, { 3, 1009, 7, 1009, 1008, 1009, 0 }, { 1008, 1009, 4, 1009, 1008, 1009, 0 }, { 8, 1009, 6, 1009, 1008, 1009, 0 }, { +1 }}), to_vec(castOf<moeb_tr>(t1)));
    EXPECT_EQ((vector<vector<int>>{{ 1, 1009, 1004, 1009, 1008, 1009, 0 }, { 3, 1009, 7, 1009, 1008, 1009, 0 }, { 1008, 1009, 4, 1009, 1008, 1009, 0 }, { 8, 1009, 6, 1009, 1008, 1009, 0 }, { -1 }}), to_vec(castOf<moeb_tr>(t3)));
}
