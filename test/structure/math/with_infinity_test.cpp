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

template<typename W>
vector<int> to_vec(const W& w) {
    return{ w.v.a.v, w.v.a.M(), w.v.b.v, w.v.b.M(), w.v.D().v, w.v.D().M(), w.is_inf };
}

TEST(with_infinity_test, constructor) {
    winf w1;
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 0}), to_vec(w1));
    winf w2(10);
    EXPECT_EQ((vector<int>{10, 1009, 0, 1009, 1008, 1009, 0}), to_vec(w2));
    winf w3(cplx(+2, -5));
    EXPECT_EQ((vector<int>{2, 1009, 1004, 1009, 1008, 1009, 0}), to_vec(w3));
    winf w4(cplx(+2, -5), 1);
    EXPECT_EQ((vector<int>{2, 1009, 1004, 1009, 1008, 1009, 1}), to_vec(w4));
    winf w5(w4);
    EXPECT_EQ((vector<int>{2, 1009, 1004, 1009, 1008, 1009, 1}), to_vec(w5));
}

TEST(with_infinity_test, operators_comparison) {
    ASSERT_COMPARISON_OPERATORS(0, winf(cplx(2, 5)), winf(cplx(2, 5)));
    ASSERT_COMPARISON_OPERATORS(+1, winf(cplx(2, 5)), winf(cplx(2, 3)));
    ASSERT_COMPARISON_OPERATORS(-1, winf(cplx(1, 5)), winf(cplx(2, 3)));
    ASSERT_COMPARISON_OPERATORS(-1, winf(cplx(2, 5)), winf(cplx(2, 5), 1));
    ASSERT_COMPARISON_OPERATORS(-1, winf(cplx(2, 5)), winf(cplx(2, 3), 1));
    ASSERT_COMPARISON_OPERATORS(-1, winf(cplx(1, 5)), winf(cplx(2, 3), 1));
    ASSERT_COMPARISON_OPERATORS(+1, winf(cplx(2, 5), 1), winf(cplx(2, 5)));
    ASSERT_COMPARISON_OPERATORS(+1, winf(cplx(2, 5), 1), winf(cplx(2, 3)));
    ASSERT_COMPARISON_OPERATORS(+1, winf(cplx(1, 5), 1), winf(cplx(2, 3)));
    ASSERT_COMPARISON_OPERATORS(0, winf(cplx(2, 5), 1), winf(cplx(2, 5), 1));
    ASSERT_COMPARISON_OPERATORS(+1, winf(cplx(2, 5), 1), winf(cplx(2, 3), 1));
    ASSERT_COMPARISON_OPERATORS(-1, winf(cplx(1, 5), 1), winf(cplx(2, 3), 1));
}

TEST(with_infinity_test, operators_arithmetic) {
    const winf wi(cplx(4, 6), 1);
    const winf w0(cplx(0));
    const winf w1(cplx(2, -5));
    const winf w2(cplx(3, 4));
    EXPECT_EQ((vector<int>{5, 1009, 1008, 1009, 1008, 1009, 0}), to_vec(w1 + w2));
    EXPECT_EQ((vector<int>{1008, 1009, 1000, 1009, 1008, 1009, 0}), to_vec(w1 - w2));
    EXPECT_EQ((vector<int>{1007, 1009, 5, 1009, 1008, 1009, 0}), to_vec(-w1));
    EXPECT_EQ((vector<int>{26, 1009, 1002, 1009, 1008, 1009, 0}), to_vec(w1 * w2));
    EXPECT_EQ((vector<int>{847, 1009, 887, 1009, 1008, 1009, 0}), to_vec(w1 / w2));
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 0}), to_vec(w1 % w2));
    EXPECT_EQ((vector<int>{6, 1009, 994, 1009, 1008, 1009, 0}), to_vec(w1 * cplx(3)));
    EXPECT_EQ((vector<int>{1, 1009, 502, 1009, 1008, 1009, 0}), to_vec(w1 / cplx(2)));

    EXPECT_EQ((vector<int>{2, 1009, 1004, 1009, 1008, 1009, 1}), to_vec(w1 + wi));  // make infinity
    EXPECT_EQ((vector<int>{6, 1009, 1, 1009, 1008, 1009, 1}), to_vec(wi + w1));     // leave infinity & add values
    EXPECT_EQ((vector<int>{4, 1009, 6, 1009, 1008, 1009, 1}), to_vec(wi + wi));     // leave infinity

    EXPECT_EQ((vector<int>{2, 1009, 1004, 1009, 1008, 1009, 1}), to_vec(w1 - wi));  // make infinity
    EXPECT_EQ((vector<int>{2, 1009, 11, 1009, 1008, 1009, 1}), to_vec(wi - w1));    // leave infinity & subtract values
    EXPECT_EQ((vector<int>{4, 1009, 6, 1009, 1008, 1009, 1}), to_vec(wi - wi));     // leave infinity (undefined by contract)

    EXPECT_EQ((vector<int>{2, 1009, 1004, 1009, 1008, 1009, 1}), to_vec(w1 * wi));  // make infinity
    EXPECT_EQ((vector<int>{38, 1009, 1001, 1009, 1008, 1009, 1}), to_vec(wi * w1)); // leave infinity & multiply values
    EXPECT_EQ((vector<int>{4, 1009, 6, 1009, 1008, 1009, 1}), to_vec(wi * wi));     // leave infinity
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 1}), to_vec(w0 * wi));     // make infinity & zero values (undefined by contract)
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 1}), to_vec(wi * w0));     // make infinity & zero values (undefined by contract)

    EXPECT_EQ((vector<int>{2, 1009, 1004, 1009, 1008, 1009, 1}), to_vec(w1 / w0));  // make infinity
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 0}), to_vec(w1 / wi));     // make zero
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 1}), to_vec(w0 / w0));     // make infinity & zero values (undefined by contract)
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 0}), to_vec(w0 / wi));     // leave zero
    EXPECT_EQ((vector<int>{208, 1009, 523, 1009, 1008, 1009, 1}), to_vec(wi / w1)); // leave infinity & divide values
    EXPECT_EQ((vector<int>{4, 1009, 6, 1009, 1008, 1009, 1}), to_vec(wi / w0));     // leave infinity
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 1}), to_vec(wi / wi));     // make infinity & zero values (undefined by contract)
}

TEST(with_infinity_test, operators_inplace) {
    const winf wi(cplx(4, 6), 1);
    const winf w0(cplx(0));
    const winf w1(cplx(2, -5));
    const winf w2(cplx(3, 4));
    winf wr;
    wr = w1; wr += w2;
    EXPECT_EQ((vector<int>{5, 1009, 1008, 1009, 1008, 1009, 0}), to_vec(wr));
    wr = w1; wr -= w2;
    EXPECT_EQ((vector<int>{1008, 1009, 1000, 1009, 1008, 1009, 0}), to_vec(wr));
    wr = w1; wr *= w2;
    EXPECT_EQ((vector<int>{26, 1009, 1002, 1009, 1008, 1009, 0}), to_vec(wr));
    wr = w1; wr /= w2;
    EXPECT_EQ((vector<int>{847, 1009, 887, 1009, 1008, 1009, 0}), to_vec(wr));
    wr = w1; wr %= w2;
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 0}), to_vec(wr));
    wr = w1; wr *= cplx(3);
    EXPECT_EQ((vector<int>{6, 1009, 994, 1009, 1008, 1009, 0}), to_vec(wr));
    wr = w1; wr /= cplx(2);
    EXPECT_EQ((vector<int>{1, 1009, 502, 1009, 1008, 1009, 0}), to_vec(wr));

    wr = w1; wr += wi;
    EXPECT_EQ((vector<int>{2, 1009, 1004, 1009, 1008, 1009, 1}), to_vec(wr));  // make infinity
    wr = wi; wr += w1;
    EXPECT_EQ((vector<int>{6, 1009, 1, 1009, 1008, 1009, 1}), to_vec(wr));     // leave infinity & add values
    wr = wi; wr += wi;
    EXPECT_EQ((vector<int>{4, 1009, 6, 1009, 1008, 1009, 1}), to_vec(wr));     // leave infinity

    wr = w1; wr -= wi;
    EXPECT_EQ((vector<int>{2, 1009, 1004, 1009, 1008, 1009, 1}), to_vec(wr));  // make infinity
    wr = wi; wr -= w1;
    EXPECT_EQ((vector<int>{2, 1009, 11, 1009, 1008, 1009, 1}), to_vec(wr));    // leave infinity & subtract values
    wr = wi; wr -= wi;
    EXPECT_EQ((vector<int>{4, 1009, 6, 1009, 1008, 1009, 1}), to_vec(wr));     // leave infinity (undefined by contract)

    wr = w1; wr *= wi;
    EXPECT_EQ((vector<int>{2, 1009, 1004, 1009, 1008, 1009, 1}), to_vec(wr));  // make infinity
    wr = wi; wr *= w1;
    EXPECT_EQ((vector<int>{38, 1009, 1001, 1009, 1008, 1009, 1}), to_vec(wr)); // leave infinity & multiply values
    wr = wi; wr *= wi;
    EXPECT_EQ((vector<int>{4, 1009, 6, 1009, 1008, 1009, 1}), to_vec(wr));     // leave infinity
    wr = w0; wr *= wi;
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 1}), to_vec(wr));     // make infinity & zero values (undefined by contract)
    wr = wi; wr *= w0;
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 1}), to_vec(wr));     // make infinity & zero values (undefined by contract)

    wr = w1; wr /= w0;
    EXPECT_EQ((vector<int>{2, 1009, 1004, 1009, 1008, 1009, 1}), to_vec(wr));  // make infinity
    wr = w1; wr /= wi;
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 0}), to_vec(wr));     // make zero
    wr = w0; wr /= w0;
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 1}), to_vec(wr));     // make infinity & zero values (undefined by contract)
    wr = w0; wr /= wi;
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 0}), to_vec(wr));     // leave zero
    wr = wi; wr /= w1;
    EXPECT_EQ((vector<int>{208, 1009, 523, 1009, 1008, 1009, 1}), to_vec(wr)); // leave infinity & divide values
    wr = wi; wr /= w0;
    EXPECT_EQ((vector<int>{4, 1009, 6, 1009, 1008, 1009, 1}), to_vec(wr));     // leave infinity
    wr = wi; wr /= wi;
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 1}), to_vec(wr));     // make infinity & zero values (undefined by contract)
}

TEST(with_infinity_test, operators_inplace_self) {
    const winf wi(cplx(4, 6), 1);
    const winf w1(cplx(2, -5));
    winf wr;
    wr = w1; wr += wr;
    EXPECT_EQ((vector<int>{4, 1009, 999, 1009, 1008, 1009, 0}), to_vec(wr));
    wr = w1; wr -= wr;
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 0}), to_vec(wr));
    wr = w1; wr *= wr;
    EXPECT_EQ((vector<int>{988, 1009, 989, 1009, 1008, 1009, 0}), to_vec(wr));
    wr = w1; wr /= wr;
    EXPECT_EQ((vector<int>{1, 1009, 0, 1009, 1008, 1009, 0}), to_vec(wr));
    wr = w1; wr %= wr;
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 0}), to_vec(wr));

    wr = wi; wr += wr;
    EXPECT_EQ((vector<int>{4, 1009, 6, 1009, 1008, 1009, 1}), to_vec(wr));      // leave infinity
    wr = wi; wr -= wr;
    EXPECT_EQ((vector<int>{4, 1009, 6, 1009, 1008, 1009, 1}), to_vec(wr));      // leave infinity (undefined by contract)
    wr = wi; wr *= wr;
    EXPECT_EQ((vector<int>{4, 1009, 6, 1009, 1008, 1009, 1}), to_vec(wr));      // leave infinity
    wr = wi; wr /= wr;
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 1}), to_vec(wr));      // make infinity & zero values (undefined by contract)
    wr = wi; wr %= wr;
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 1}), to_vec(wr));      // make infinity & zero values (undefined by contract)
}

TEST(with_infinity_test, inverse) {
    const winf wi(cplx(4, 6), 1);
    const winf w0(cplx(0));
    const winf w1(cplx(2, -5));
    EXPECT_EQ((vector<int>{348, 1009, 870, 1009, 1008, 1009, 0}), to_vec(w1.inverse()));
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 1}), to_vec(w0.inverse()));
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 0}), to_vec(wi.inverse()));
}

TEST(with_infinity_test, conjugate) {
    const winf wi(cplx(4, 6), 1);
    const winf w0(cplx(0));
    const winf w1(cplx(2, -5));
    EXPECT_EQ((vector<int>{2, 1009, 5, 1009, 1008, 1009, 0}), to_vec(conjugateT<winf>::of(w1)));
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 0}), to_vec(conjugateT<winf>::of(w0)));
    EXPECT_EQ((vector<int>{4, 1009, 1003, 1009, 1008, 1009, 1}), to_vec(conjugateT<winf>::of(wi)));
}

TEST(with_infinity_test, casts) {
    const winf wi(cplx(4, 6), 1);
    const winf w2(cplx(2, -5));
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 0}), to_vec(zeroT<winf>::of(w2)));
    EXPECT_EQ((vector<int>{1, 1009, 0, 1009, 1008, 1009, 0}), to_vec(identityT<winf>::of(w2)));
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 1008, 1009, 0}), to_vec(zeroT<winf>::of(wi)));
    EXPECT_EQ((vector<int>{1, 1009, 0, 1009, 1008, 1009, 0}), to_vec(identityT<winf>::of(wi)));
    EXPECT_EQ((vector<int>{3, 1009, 0, 1009, 1008, 1009, 0}), to_vec(castOf<winf>(3)));
    EXPECT_EQ((vector<int>{4, 1009, 0, 1009, 1008, 1009, 0}), to_vec(castOf<winf>(w2, 4)));
    EXPECT_EQ((vector<int>{4, 1009, 0, 1009, 1008, 1009, 0}), to_vec(castOf<winf>(wi, 4)));
    EXPECT_EQ((vector<int>{2, 1009, 1004, 1009, 1008, 1009, 0}), to_vec(castOf<winf>(wi, w2)));
    EXPECT_EQ((vector<int>{4, 1009, 6, 1009, 1008, 1009, 1}), to_vec(castOf<winf>(w2, wi)));
    EXPECT_EQ((vector<int>{2, 1009, 1004, 1009, 1008, 1009, 0}), to_vec(castOf<winf>(w2)));
    EXPECT_EQ((vector<int>{4, 1009, 6, 1009, 1008, 1009, 1}), to_vec(castOf<winf>(wi)));
}
