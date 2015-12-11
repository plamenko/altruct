#include "structure/math/polynom.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

TEST(polynom_test, constructor) {
	vector<int> c = { 1, 2, 3, 4 };
	polynom<int> p1;
	EXPECT_EQ((vector<int>{}), p1.c);
	polynom<int> p2(c);
	EXPECT_EQ(c, p2.c);
	polynom<int> p3(p2);
	EXPECT_EQ(c, p3.c);
	polynom<int> p4(c.begin(), c.end());
	EXPECT_EQ(c, p4.c);
	polynom<int> p5(&c[0], &c[0] + c.size());
	EXPECT_EQ(c, p5.c);
	polynom<int> p6(5);
	EXPECT_EQ((vector<int>{ 5 }), p6.c);
	polynom<int> p7{ 1, 2, 3, 4 };
	EXPECT_EQ(c, p7.c);
}

TEST(polynom_test, swap) {
	polynom<int> p1{ 1, 2, 3, 4 };
	polynom<int> p2{ 5, 6, 7 };
	p1.swap(p2);
	EXPECT_EQ((vector<int>{ 5, 6, 7 }), p1.c);
	EXPECT_EQ((vector<int>{ 1, 2, 3, 4 }), p2.c);
}

TEST(polynom_test, shrink_to_fit) {
	polynom<int> p{ 1, 2, 3, 4, 0, 0 };
	EXPECT_EQ(6, p.c.size());
	p.shrink_to_fit();
	EXPECT_EQ(4, p.c.size());
	EXPECT_EQ((vector<int>{ 1, 2, 3, 4 }), p.c);
}

TEST(polynom_test, reserve) {
	polynom<int> p{ 1, 2, 3, 4 };
	EXPECT_EQ(4, p.c.size());
	p.reserve(6);
	EXPECT_EQ(6, p.c.size());
	EXPECT_EQ((vector<int>{ 1, 2, 3, 4, 0, 0 }), p.c);
}

TEST(polynom_test, size) {
	polynom<int> p{ 1, 2, 3, 4 };
	EXPECT_EQ(4, p.size());
}

TEST(polynom_test, at) {
	const polynom<int> p{ 2, 3, 5, 7 };
	EXPECT_EQ(2, p.at(0));
	EXPECT_EQ(7, p.at(3));
	EXPECT_EQ(0, p.at(4));
	EXPECT_EQ(0, p.at(100));
	EXPECT_EQ(4, p.size());
}

TEST(polynom_test, operator_const_brackets) {
	const polynom<int> p{ 2, 3, 5, 7 };
	EXPECT_EQ(2, p[0]);
	EXPECT_EQ(7, p[3]);
	EXPECT_EQ(0, p[4]);
	EXPECT_EQ(0, p[100]);
	EXPECT_EQ(4, p.size());
}

TEST(polynom_test, operator_brackets) {
	polynom<int> p;
	p[3] = 3;
	EXPECT_EQ(0, p[0]);
	EXPECT_EQ(0, p[4]);
	EXPECT_EQ(3, p[3]);
	EXPECT_EQ(0, p[4]);
	EXPECT_EQ(0, p[100]);
	EXPECT_EQ(101, p.size());
}

TEST(polynom_test, degree) {
	const polynom<int> p1;
	EXPECT_EQ(0, p1.deg());
	const polynom<int> p2{ 4 };
	EXPECT_EQ(0, p2.deg());
	const polynom<int> p3{ 0, 3 };
	EXPECT_EQ(1, p3.deg());
	const polynom<int> p4{ 2, 3, 5, 7 };
	EXPECT_EQ(3, p4.deg());
	const polynom<int> p5{ 2, 3, 5, 7, 0, 0 };
	EXPECT_EQ(3, p5.deg());
}

TEST(polynom_test, leading_coefficient) {
	const polynom<int> p1;
	EXPECT_EQ(0, p1.leading_coeff());
	const polynom<int> p2{ 4 };
	EXPECT_EQ(4, p2.leading_coeff());
	const polynom<int> p3{ 0, 3 };
	EXPECT_EQ(3, p3.leading_coeff());
	const polynom<int> p4{ 2, 3, 5, 7 };
	EXPECT_EQ(7, p4.leading_coeff());
	const polynom<int> p5{ 2, 3, 5, 7, 0, 0 };
	EXPECT_EQ(7, p5.leading_coeff());
}

TEST(polynom_test, is_power) {
	const polynom<int> p1;
	EXPECT_EQ(false, p1.is_power());
	const polynom<int> p2{ 4 };
	EXPECT_EQ(false, p2.is_power());
	const polynom<int> p3{ 1 };
	EXPECT_EQ(true, p3.is_power());
	const polynom<int> p4{ 0, 0, 0, 3 };
	EXPECT_EQ(false, p4.is_power());
	const polynom<int> p5{ 0, 0, 0, 1 };
	EXPECT_EQ(true, p5.is_power());
	const polynom<int> p6{ 0, 0, 0, 1, 0 };
	EXPECT_EQ(true, p6.is_power());
}
