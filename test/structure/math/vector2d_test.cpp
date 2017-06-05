#include "altruct/structure/math/vector2d.h"

#include "structure_test_util.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;
using namespace altruct::test_util;

typedef vector2d<int> pnti;
typedef vector2d<double, int> pntd;

TEST(vector2d_test, constructor) {
	pnti p1;
	EXPECT_EQ(0, p1.x);
	EXPECT_EQ(0, p1.y);
	pnti p2(+3, -5);
	EXPECT_EQ(+3, p2.x);
	EXPECT_EQ(-5, p2.y);
	pnti p3(p2);
	EXPECT_EQ(+3, p3.x);
	EXPECT_EQ(-5, p3.y);
	pntd p4;
	EXPECT_EQ(0.0, p4.x);
	EXPECT_EQ(0.0, p4.y);
	EXPECT_EQ(0, p4.data);
	pntd p5(+3.5, -5.2);
	EXPECT_EQ(+3.5, p5.x);
	EXPECT_EQ(-5.2, p5.y);
	EXPECT_EQ(0, p5.data);
	pntd p6(+3.5, -5.2, 4);
	EXPECT_EQ(+3.5, p6.x);
	EXPECT_EQ(-5.2, p6.y);
	EXPECT_EQ(4, p6.data);
	pntd p7{ +3.5, -5.2, 4 };
	EXPECT_EQ(+3.5, p7.x);
	EXPECT_EQ(-5.2, p7.y);
	EXPECT_EQ(4, p7.data);
	pntd p8(p7);
	EXPECT_EQ(+3.5, p8.x);
	EXPECT_EQ(-5.2, p8.y);
	EXPECT_EQ(4, p8.data);
}

TEST(vector2d_test, operators_comparison) {
	const pntd p1{ 1, 2 };
	const pntd p2{ 0, 5 };
	const pntd p3{ 1, 0 };
	const pntd p4{ 5, 0 };
	const pntd p5{ 1, 5 };
	ASSERT_COMPARISON_OPERATORS(0, p1, p1);
	ASSERT_COMPARISON_OPERATORS(+1, p1, p2);
	ASSERT_COMPARISON_OPERATORS(-1, p2, p1);
	ASSERT_COMPARISON_OPERATORS(+1, p1, p3);
	ASSERT_COMPARISON_OPERATORS(-1, p3, p1);
	ASSERT_COMPARISON_OPERATORS(-1, p1, p4);
	ASSERT_COMPARISON_OPERATORS(+1, p4, p1);
	ASSERT_COMPARISON_OPERATORS(-1, p1, p5);
	ASSERT_COMPARISON_OPERATORS(+1, p5, p1);
}

TEST(vector2d_test, operators_arithmetic) {
	const pntd p1{ 1, 2 };
	const pntd p2{ -3, 5 };
	const pntd p3{ 8, -3 };
	EXPECT_EQ((pntd{ -2, 7 }), p1 + p2);
	EXPECT_EQ((pntd{ 4, -3 }), p1 - p2);
	EXPECT_EQ((pntd{ -1, -2 }), -p1);
	EXPECT_EQ((pntd{ -3, 10 }), p1 * p2);
	EXPECT_EQ((pntd{ -3, 2.5}), p2 / p1);
	EXPECT_EQ((pntd{ -3, -6 }), p1 * -3);
	EXPECT_EQ((pntd{ 0.5, 1 }), p1 / 2);
	EXPECT_EQ(-39, p2 & p3);
	EXPECT_EQ(-31, p2 ^ p3);
	EXPECT_EQ(-43, p1.dot(p2, p3));
	EXPECT_EQ(-1, p1.cross(p2, p3));
}

TEST(vector2d_test, operators_inplace) {
	const pntd p1{ 1, 2 };
	const pntd p2{ -3, 5 };
	const pntd p3{ 8, -3 };
	pntd pr;
	pr = p1; pr += p2;
	EXPECT_EQ((pntd{ -2, 7 }), pr);
	pr = p1; pr -= p2;
	EXPECT_EQ((pntd{ 4, -3 }), pr);
	pr = p1; pr *= p2;
	EXPECT_EQ((pntd{ -3, 10 }), pr);
	pr = p2; pr /= p1;
	EXPECT_EQ((pntd{ -3, 2.5 }), pr);
	pr = p1; pr *= -3;
	EXPECT_EQ((pntd{ -3, -6 }), pr);
	pr = p1; pr /= 2;
	EXPECT_EQ((pntd{ 0.5, 1 }), pr);
}

TEST(vector2d_test, operators_inplace_self) {
	const pntd p1{ -3, 5 };
	pntd pr;
	pr = p1; pr += pr;
	EXPECT_EQ((pntd{ -6, 10 }), pr);
	pr = p1; pr -= pr;
	EXPECT_EQ((pntd{ 0, 0 }), pr);
	pr = p1; pr *= pr;
	EXPECT_EQ((pntd{ 9, 25 }), pr);
	pr = p1; pr /= pr;
	EXPECT_EQ((pntd{ 1, 1 }), pr);
}

TEST(vector2d_test, other) {
	const pntd pE{ 888, 888 };
	const pntd p0{ 0, 0 };
	const pntd p1{ -3, 4 };
	const pntd p2{ 8, 3 };
	EXPECT_EQ(p0, p0.unit());
	EXPECT_EQ(pE, p0.unit(pE));
	EXPECT_EQ((pntd{ -3, 4 } / 5), p1.unit());
	EXPECT_EQ((pntd{ -3, 4 } / 5), p1.unit(pE));
	EXPECT_EQ((pntd{ -36, 23 }), p1.rot(p2));
	EXPECT_EQ((pntd{ -12, 41 }), p1.irot(p2));
	EXPECT_EQ((pntd{ -4, -3 }), p1.perp());
	EXPECT_EQ(5, p1.abs1());
	EXPECT_EQ(25, p1.abs2());
	EXPECT_EQ(55, p2.diff2());
}
