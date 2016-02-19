#include "structure/math/fraction.h"

#include "gtest/gtest.h"

//#include <functional>

using namespace std;
using namespace altruct::math;

typedef fraction<int> frac;

TEST(fraction_test, constructor) {
	frac f1;
	EXPECT_EQ(0, f1.p);
	EXPECT_EQ(1, f1.q);
	frac f2(10);
	EXPECT_EQ(10, f2.p);
	EXPECT_EQ(1, f2.q);
	frac f3(10, 6);
	EXPECT_EQ(5, f3.p);
	EXPECT_EQ(3, f3.q);
	frac f4(f3);
	EXPECT_EQ(5, f4.p);
	EXPECT_EQ(3, f4.q);
	frac f5(10, -6);
	EXPECT_EQ(-5, f5.p);
	EXPECT_EQ(3, f5.q);
}

TEST(fraction_test, operators_comparison) {
	const frac f1(20, 31);
	const frac f2(3, 4);
	EXPECT_EQ(false, f1 == f2);
	EXPECT_EQ(true, f1 != f2);
	EXPECT_EQ(true, f1 < f2);
	EXPECT_EQ(false, f1 > f2);
	EXPECT_EQ(true, f1 <= f2);
	EXPECT_EQ(false, f1 >= f2);
	EXPECT_EQ(false, f2 == f1);
	EXPECT_EQ(true, f2 != f1);
	EXPECT_EQ(false, f2 < f1);
	EXPECT_EQ(true, f2 > f1);
	EXPECT_EQ(false, f2 <= f1);
	EXPECT_EQ(true, f2 >= f1);
	EXPECT_EQ(true, f2 == f2);
	EXPECT_EQ(false, f2 != f2);
	EXPECT_EQ(false, f2 < f2);
	EXPECT_EQ(false, f2 > f2);
	EXPECT_EQ(true, f2 <= f2);
	EXPECT_EQ(true, f2 >= f2);
}

TEST(fraction_test, operators_arithmetic) {
	const frac f1(5, 6);
	const frac f2(3, 10);
	EXPECT_EQ(frac(17, 15), f1 + f2);
	EXPECT_EQ(frac(8, 15), f1 - f2);
	EXPECT_EQ(frac(-5, 6), -f1);
	EXPECT_EQ(frac(1, 4), f1 * f2);
	EXPECT_EQ(frac(25, 9), f1 / f2);
	EXPECT_EQ(frac(0), f1 % f2);
	EXPECT_EQ(frac(17, 15), f2 + f1);
	EXPECT_EQ(frac(-8, 15), f2 - f1);
	EXPECT_EQ(frac(-3, 10), -f2);
	EXPECT_EQ(frac(1, 4), f2 * f1);
	EXPECT_EQ(frac(9, 25), f2 / f1);
	EXPECT_EQ(frac(0), f2 % f1);
}

TEST(fraction_test, operators_inplace) {
	const frac f1(5, 6);
	const frac f2(3, 10);
	frac fr;
	fr = f1; fr += f2;
	EXPECT_EQ(frac(17, 15), fr);
	fr = f1; fr -= f2;
	EXPECT_EQ(frac(8, 15), fr);
	fr = f1; fr *= f2;
	EXPECT_EQ(frac(1, 4), fr);
	fr = f1; fr /= f2;
	EXPECT_EQ(frac(25, 9), fr);
	fr = f1; fr %= f2;
	EXPECT_EQ(frac(0), fr);
	fr = f2; fr += f1;
	EXPECT_EQ(frac(17, 15), fr);
	fr = f2; fr -= f1;
	EXPECT_EQ(frac(-8, 15), fr);
	fr = f2; fr *= f1;
	EXPECT_EQ(frac(1, 4), fr);
	fr = f2; fr /= f1;
	EXPECT_EQ(frac(9, 25), fr);
	fr = f2; fr %= f1;
	EXPECT_EQ(frac(0), fr);
}

TEST(fraction_test, operators_inplace_self) {
	const frac f1(3, 7);
	frac fr;
	fr = f1; fr += fr;
	EXPECT_EQ(frac(6, 7), fr);
	fr = f1; fr -= fr;
	EXPECT_EQ(frac(0), fr);
	fr = f1; fr *= fr;
	EXPECT_EQ(frac(9, 49), fr);
	fr = f1;
	fr /= fr;
	EXPECT_EQ(frac(1), fr);
	fr = f1; fr %= fr;
	EXPECT_EQ(frac(0), fr);
}
