#include "structure/math/galois_field_2.h"

#include "structure_test_util.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;
using namespace altruct::test_util;

typedef galois_field_2<> gf2;

TEST(galois_field_2, constructor) {
	gf2 x1;
	EXPECT_EQ(0, x1.v);
	gf2 x2(123);
	EXPECT_EQ(123, x2.v);
}

TEST(galois_field_2, operators_comparison) {
	ASSERT_COMPARISON_OPERATORS(0, gf2(123), gf2(123));
	ASSERT_COMPARISON_OPERATORS(-1, gf2(123), gf2(241));
	ASSERT_COMPARISON_OPERATORS(+1, gf2(241), gf2(123));
}

TEST(galois_field_2, operators_arithmetic) {
	const gf2 x1(123);
	const gf2 x2(241);
	EXPECT_EQ(gf2(123 ^ 241), x1 + x2);
	EXPECT_EQ(gf2(123 ^ 241), x1 - x2);
	EXPECT_EQ(gf2(123), -x1);
	EXPECT_EQ(gf2(123 & 241), x1 * x2);
	EXPECT_EQ(gf2(123 | ~241), x1 / x2);
	EXPECT_EQ(gf2(241 ^ 123), x2 + x1);
	EXPECT_EQ(gf2(241 ^ 123), x2 - x1);
	EXPECT_EQ(gf2(241), -x2);
	EXPECT_EQ(gf2(241 & 123), x2 * x1);
	EXPECT_EQ(gf2(241 | ~123), x2 / x1);
}

TEST(galois_field_2, operators_inplace) {
	const gf2 x1(123);
	const gf2 x2(241);
	gf2 xr;
	xr = x1; xr += x2;
	EXPECT_EQ(gf2(123 ^ 241), xr);
	xr = x1; xr -= x2;
	EXPECT_EQ(gf2(123 ^ 241), xr);
	xr = x1; xr *= x2;
	EXPECT_EQ(gf2(123 & 241), xr);
	xr = x1; xr /= x2;
	EXPECT_EQ(gf2(123 | ~241), xr);
	xr = x2; xr += x1;
	EXPECT_EQ(gf2(241 ^ 123), xr);
	xr = x2; xr -= x1;
	EXPECT_EQ(gf2(241 ^ 123), xr);
	xr = x2; xr *= x1;
	EXPECT_EQ(gf2(241 & 123), xr);
	xr = x2; xr /= x1;
	EXPECT_EQ(gf2(241 | ~123), xr);
}

TEST(galois_field_2, operators_inplace_self) {
	const gf2 x1(123);
	gf2 xr;
	xr = x1; xr += xr;
	EXPECT_EQ(gf2(0), xr);
	xr = x1; xr -= xr;
	EXPECT_EQ(gf2(0), xr);
	xr = x1; xr *= xr;
	EXPECT_EQ(gf2(123), xr);
	xr = x1; xr /= xr;
	EXPECT_EQ(gf2(~0), xr);
}

TEST(galois_field_2, identity) {
	const gf2 x(123);
	const gf2 e0 = zeroT<gf2>::of(x);
	const gf2 e1 = identityT<gf2>::of(x);
	EXPECT_EQ(0, e0.v);
	EXPECT_EQ(~0, e1.v);
}
