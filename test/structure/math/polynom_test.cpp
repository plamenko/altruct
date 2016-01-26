#include "structure/math/polynom.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

TEST(polynom_test, constructor) {
	const vector<int> c = { 1, 2, 3, 4 };
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

TEST(polynom_test, cmp) {
	const polynom<int> p0{ };
	const polynom<int> p1{ 4 };
	const polynom<int> p2{ 1, 3, 5, 7 };
	const polynom<int> p3{ 2, 3, 5, 7 };
	const polynom<int> p4{ 2, 3, 5, 7, 0, 0 };
	const polynom<int> p5{ 2, 3, 6, 7 };
	EXPECT_EQ(0, polynom<int>::cmp(p0, p0));
	EXPECT_EQ(-1, polynom<int>::cmp(p0, p1));
	EXPECT_EQ(+1, polynom<int>::cmp(p1, p0));
	EXPECT_EQ(0, polynom<int>::cmp(p1, p1));
	EXPECT_EQ(-1, polynom<int>::cmp(p0, p2));
	EXPECT_EQ(+1, polynom<int>::cmp(p2, p0));
	EXPECT_EQ(-1, polynom<int>::cmp(p1, p2));
	EXPECT_EQ(+1, polynom<int>::cmp(p2, p1));
	EXPECT_EQ(0, polynom<int>::cmp(p2, p2));
	EXPECT_EQ(-1, polynom<int>::cmp(p2, p3));
	EXPECT_EQ(+1, polynom<int>::cmp(p3, p2));
	EXPECT_EQ(0, polynom<int>::cmp(p3, p4));
	EXPECT_EQ(0, polynom<int>::cmp(p4, p3));
	EXPECT_EQ(-1, polynom<int>::cmp(p4, p5));
	EXPECT_EQ(+1, polynom<int>::cmp(p5, p4));
}

TEST(polynom_test, neg) {
	const polynom<int> p0{};
	const polynom<int> p1{ 4 };
	const polynom<int> p2{ 1, -3, 5, 7 };
	const polynom<int> p3{ 2, 3, 5, -7, 0, 0 };
	polynom<int> pr;
	polynom<int>::neg(pr, p0);
	EXPECT_EQ((polynom<int>{}), pr);
	polynom<int>::neg(pr, p1);
	EXPECT_EQ((polynom<int>{ -4 }), pr);
	polynom<int>::neg(pr, p2);
	EXPECT_EQ((polynom<int>{ -1, 3, -5, -7 }), pr);
	polynom<int>::neg(pr, p3);
	EXPECT_EQ((polynom<int>{ -2, -3, -5, 7 }), pr);
	// inplace
	polynom<int>::neg(pr, pr);
	EXPECT_EQ((polynom<int>{ 2, 3, 5, -7 }), pr);
}

TEST(polynom_test, add) {
	const polynom<int> p0{};
	const polynom<int> p1{ 4 };
	const polynom<int> p2{ 1, -3, 5, 7 };
	const polynom<int> p3{ 2, 3, 5, -7, 0, 0 };
	polynom<int> pr;
	polynom<int>::add(pr, p0, p0);
	EXPECT_EQ((polynom<int>{}), pr);
	polynom<int>::add(pr, p0, p1);
	EXPECT_EQ((polynom<int>{ 4 }), pr);
	polynom<int>::add(pr, p1, p0);
	EXPECT_EQ((polynom<int>{ 4 }), pr);
	polynom<int>::add(pr, p1, p1);
	EXPECT_EQ((polynom<int>{ 8 }), pr);
	polynom<int>::add(pr, p0, p2);
	EXPECT_EQ((polynom<int>{ 1, -3, 5, 7 }), pr);
	polynom<int>::add(pr, p2, p0);
	EXPECT_EQ((polynom<int>{ 1, -3, 5, 7 }), pr);
	polynom<int>::add(pr, p1, p2);
	EXPECT_EQ((polynom<int>{ 5, -3, 5, 7 }), pr);
	polynom<int>::add(pr, p2, p1);
	EXPECT_EQ((polynom<int>{ 5, -3, 5, 7 }), pr);
	polynom<int>::add(pr, p2, p3);
	EXPECT_EQ((polynom<int>{ 3, 0, 10 }), pr);
	polynom<int>::add(pr, p3, p2);
	EXPECT_EQ((polynom<int>{ 3, 0, 10 }), pr);
	polynom<int>::add(pr, p3, p3);
	EXPECT_EQ((polynom<int>{ 4, 6, 10, -14 }), pr);
	// inplace
	pr = { 4, 6, 10, -14 };
	polynom<int>::add(pr, pr, p1);
	EXPECT_EQ((polynom<int>{ 8, 6, 10, -14 }), pr);
	polynom<int>::add(pr, p1, pr);
	EXPECT_EQ((polynom<int>{ 12, 6, 10, -14 }), pr);
	polynom<int>::add(pr, pr, pr);
	EXPECT_EQ((polynom<int>{ 24, 12, 20, -28 }), pr);
}

TEST(polynom_test, sub) {
	const polynom<int> p0{};
	const polynom<int> p1{ 4 };
	const polynom<int> p2{ 1, -3, 5, 7 };
	const polynom<int> p3{ 2, 3, 5, -7, 0, 0 };
	polynom<int> pr;
	polynom<int>::sub(pr, p0, p0);
	EXPECT_EQ((polynom<int>{}), pr);
	polynom<int>::sub(pr, p0, p1);
	EXPECT_EQ((polynom<int>{ -4 }), pr);
	polynom<int>::sub(pr, p1, p0);
	EXPECT_EQ((polynom<int>{ 4 }), pr);
	polynom<int>::sub(pr, p1, p1);
	EXPECT_EQ((polynom<int>{ }), pr);
	polynom<int>::sub(pr, p0, p2);
	EXPECT_EQ((polynom<int>{ -1, 3, -5, -7 }), pr);
	polynom<int>::sub(pr, p2, p0);
	EXPECT_EQ((polynom<int>{ 1, -3, 5, 7 }), pr);
	polynom<int>::sub(pr, p1, p2);
	EXPECT_EQ((polynom<int>{ 3, 3, -5, -7 }), pr);
	polynom<int>::sub(pr, p2, p1);
	EXPECT_EQ((polynom<int>{ -3, -3, 5, 7 }), pr);
	polynom<int>::sub(pr, p2, p3);
	EXPECT_EQ((polynom<int>{ -1, -6, 0, 14 }), pr);
	polynom<int>::sub(pr, p3, p2);
	EXPECT_EQ((polynom<int>{ 1, 6, 0, -14 }), pr);
	polynom<int>::sub(pr, p3, p3);
	EXPECT_EQ((polynom<int>{ }), pr);
	// inplace
	pr = { 1, 6, 10, -14 };
	polynom<int>::sub(pr, pr, p1);
	EXPECT_EQ((polynom<int>{ -3, 6, 10, -14 }), pr);
	polynom<int>::sub(pr, p1, pr);
	EXPECT_EQ((polynom<int>{ 7, -6, -10, +14 }), pr);
	polynom<int>::sub(pr, pr, pr);
	EXPECT_EQ((polynom<int>{}), pr);
}

TEST(polynom_test, mul) {
	const polynom<int> p0{};
	const polynom<int> p1{ 4 };
	const polynom<int> p2{ 1, -3, 5, 7 };
	const polynom<int> p3{ 2, 3, 5, -7, 0, 0 };
	polynom<int> pr;
	polynom<int>::mul(pr, p0, p0);
	EXPECT_EQ((polynom<int>{}), pr);
	polynom<int>::mul(pr, p0, p1);
	EXPECT_EQ((polynom<int>{}), pr);
	polynom<int>::mul(pr, p1, p0);
	EXPECT_EQ((polynom<int>{}), pr);
	polynom<int>::mul(pr, p1, p1);
	EXPECT_EQ((polynom<int>{ 16 }), pr);
	polynom<int>::mul(pr, p0, p2);
	EXPECT_EQ((polynom<int>{}), pr);
	polynom<int>::mul(pr, p2, p0);
	EXPECT_EQ((polynom<int>{}), pr);
	polynom<int>::mul(pr, p1, p2);
	EXPECT_EQ((polynom<int>{ 4, -12, 20, 28 }), pr);
	polynom<int>::mul(pr, p2, p1);
	EXPECT_EQ((polynom<int>{ 4, -12, 20, 28 }), pr);
	polynom<int>::mul(pr, p2, p3);
	EXPECT_EQ((polynom<int>{ 2, -3, 6, 7, 67, 0, -49 }), pr);
	polynom<int>::mul(pr, p3, p2);
	EXPECT_EQ((polynom<int>{ 2, -3, 6, 7, 67, 0, -49 }), pr);
	polynom<int>::mul(pr, p3, p3);
	EXPECT_EQ((polynom<int>{ 4, 12, 29, 2, -17, -70, 49 }), pr);
	// inplace
	pr = { 2, 3, 5, -7 };
	polynom<int>::mul(pr, pr, p1);
	EXPECT_EQ((polynom<int>{ 8, 12, 20, -28 }), pr);
	polynom<int>::mul(pr, p1, pr);
	EXPECT_EQ((polynom<int>{ 32, 48, 80, -112 }), pr);
	polynom<int>::mul(pr, pr, pr);
	EXPECT_EQ((polynom<int>{ 1024, 3072, 7424, 512, -4352, -17920, 12544 }), pr);
}

TEST(polynom_test, quot_rem) {
	const polynom<int> p0{};
	const polynom<int> p1{ 6 };
	const polynom<int> p2{ 1, -3, 0, -2, 0, 0 };
	const polynom<int> p3{ 12, 18, 30, -42, 36, 0, 24, 0 };
	polynom<int> pr;
	polynom<int>::quot_rem(pr, p0, p1);
	EXPECT_EQ((polynom<int>{}), pr);
	polynom<int>::quot_rem(pr, p0, p2);
	EXPECT_EQ((polynom<int>{}), pr);
	polynom<int>::quot_rem(pr, p1, p1);
	EXPECT_EQ((polynom<int>{ 1 }), pr);
	polynom<int>::quot_rem(pr, p1, p2);
	EXPECT_EQ((polynom<int>{ 6 }), pr);
	polynom<int>::quot_rem(pr, p3, p1);
	EXPECT_EQ((polynom<int>{ 2, 3, 5, -7, 6, 0, 4 }), pr);
	polynom<int>::quot_rem(pr, p3, p2);
	EXPECT_EQ((polynom<int>{-3, 63, 30, 15, 0, 0, -12}), pr);
	polynom<int>::quot_rem(pr, p2, p3);
	EXPECT_EQ((polynom<int>{ 1, -3, 0, -2 }), pr);
	// inplace
	pr = p3;
	polynom<int>::quot_rem(pr, pr, p2);
	EXPECT_EQ((polynom<int>{-3, 63, 30, 15, 0, 0, -12}), pr);
	pr = p2;
	polynom<int>::quot_rem(pr, pr, p3);
	EXPECT_EQ((polynom<int>{ 1, -3, 0, -2 }), pr);
	pr = p3;
	polynom<int>::quot_rem(pr, pr, p3);
	EXPECT_EQ((polynom<int>{ 0, 0, 0, 0, 0, 0, 1 }), pr);
}

TEST(polynom_test, div) {
	const polynom<int> p0{};
	const polynom<int> p1{ 6 };
	const polynom<int> p2{ 1, -3, 0, -2, 0, 0 };
	const polynom<int> p3{ 12, 18, 30, -42, 36, 0, 24, 0 };
	polynom<int> pr;
	polynom<int>::div(pr, p0, p1);
	EXPECT_EQ((polynom<int>{}), pr);
	polynom<int>::div(pr, p0, p2);
	EXPECT_EQ((polynom<int>{}), pr);
	polynom<int>::div(pr, p1, p1);
	EXPECT_EQ((polynom<int>{ 1 }), pr);
	polynom<int>::div(pr, p1, p2);
	EXPECT_EQ((polynom<int>{ }), pr);
	polynom<int>::div(pr, p3, p1);
	EXPECT_EQ((polynom<int>{ 2, 3, 5, -7, 6, 0, 4 }), pr);
	polynom<int>::div(pr, p3, p2);
	EXPECT_EQ((polynom<int>{ 15, 0, 0, -12}), pr);
	polynom<int>::div(pr, p2, p3);
	EXPECT_EQ((polynom<int>{}), pr);
	// inplace
	pr = p3;
	polynom<int>::div(pr, pr, p2);
	EXPECT_EQ((polynom<int>{ 15, 0, 0, -12}), pr);
	pr = p2;
	polynom<int>::div(pr, pr, p3);
	EXPECT_EQ((polynom<int>{}), pr);
	pr = p3;
	polynom<int>::div(pr, pr, p3);
	EXPECT_EQ((polynom<int>{ 1 }), pr);
}

TEST(polynom_test, mod) {
	const polynom<int> p0{};
	const polynom<int> p1{ 6 };
	const polynom<int> p2{ 1, -3, 0, -2, 0, 0 };
	const polynom<int> p3{ 12, 18, 30, -42, 36, 0, 24, 0 };
	polynom<int> pr;
	polynom<int>::mod(pr, p0, p1);
	EXPECT_EQ((polynom<int>{}), pr);
	polynom<int>::mod(pr, p0, p2);
	EXPECT_EQ((polynom<int>{}), pr);
	polynom<int>::mod(pr, p1, p1);
	EXPECT_EQ((polynom<int>{}), pr);
	polynom<int>::mod(pr, p1, p2);
	EXPECT_EQ((polynom<int>{ 6 }), pr);
	polynom<int>::mod(pr, p3, p1);
	EXPECT_EQ((polynom<int>{}), pr);
	polynom<int>::mod(pr, p3, p2);
	EXPECT_EQ((polynom<int>{-3, 63, 30 }), pr);
	polynom<int>::mod(pr, p2, p3);
	EXPECT_EQ((polynom<int>{ 1, -3, 0, -2 }), pr);
	// inplace
	pr = p3;
	polynom<int>::mod(pr, pr, p2);
	EXPECT_EQ((polynom<int>{-3, 63, 30 }), pr);
	pr = p2;
	polynom<int>::mod(pr, pr, p3);
	EXPECT_EQ((polynom<int>{ 1, -3, 0, -2 }), pr);
	pr = p3;
	polynom<int>::mod(pr, pr, p3);
	EXPECT_EQ((polynom<int>{ 0, 0, 0, 0, 0, 0 }), pr);
}

TEST(polynom_test, muls) {
	const polynom<int> p0{};
	const polynom<int> p1{ 4 };
	const polynom<int> p2{ 1, -3, 5, 7 };
	const polynom<int> p3{ 2, 3, 5, -7, 0, 0 };
	polynom<int> pr;
	polynom<int>::mul(pr, p0, 11);
	EXPECT_EQ((polynom<int>{}), pr);
	polynom<int>::mul(pr, p1, 11);
	EXPECT_EQ((polynom<int>{ 44 }), pr);
	polynom<int>::mul(pr, p2, 11);
	EXPECT_EQ((polynom<int>{ 11, -33, 55, 77 }), pr);
	polynom<int>::mul(pr, p3, 11);
	EXPECT_EQ((polynom<int>{ 22, 33, 55, -77 }), pr);
	// inplace
	pr = { 2, 3, 5, -7 };
	polynom<int>::mul(pr, pr, 11);
	EXPECT_EQ((polynom<int>{ 22, 33, 55, -77 }), pr);
}

TEST(polynom_test, divs) {
	const polynom<int> p0{};
	const polynom<int> p1{ 44 };
	const polynom<int> p2{ 11, -33, 55, 77 };
	const polynom<int> p3{ 22, 33, 55, -77, 0, 0 };
	polynom<int> pr;
	polynom<int>::div(pr, p0, 11);
	EXPECT_EQ((polynom<int>{}), pr);
	polynom<int>::div(pr, p1, 11);
	EXPECT_EQ((polynom<int>{ 4 }), pr);
	polynom<int>::div(pr, p2, 11);
	EXPECT_EQ((polynom<int>{ 1, -3, 5, 7 }), pr);
	polynom<int>::div(pr, p3, 11);
	EXPECT_EQ((polynom<int>{ 2, 3, 5, -7 }), pr);
	// inplace
	pr = { 22, 33, 55, -77 };
	polynom<int>::div(pr, pr, 11);
	EXPECT_EQ((polynom<int>{ 2, 3, 5, -7 }), pr);
}

TEST(polynom_test, operators_comparison) {
	const polynom<int> p1{ 4 };
	const polynom<int> p2{ 1, 3, 5, 7 };
	EXPECT_EQ(false, p1 == p2);
	EXPECT_EQ(true, p1 != p2);
	EXPECT_EQ(true, p1 < p2);
	EXPECT_EQ(false, p1 > p2);
	EXPECT_EQ(true, p1 <= p2);
	EXPECT_EQ(false, p1 >= p2);
	EXPECT_EQ(false, p2 == p1);
	EXPECT_EQ(true, p2 != p1);
	EXPECT_EQ(false, p2 < p1);
	EXPECT_EQ(true, p2 > p1);
	EXPECT_EQ(false, p2 <= p1);
	EXPECT_EQ(true, p2 >= p1);
	EXPECT_EQ(true, p2 == p2);
	EXPECT_EQ(false, p2 != p2);
	EXPECT_EQ(false, p2 < p2);
	EXPECT_EQ(false, p2 > p2);
	EXPECT_EQ(true, p2 <= p2);
	EXPECT_EQ(true, p2 >= p2);
}

TEST(polynom_test, operators_arithmetic) {
	const polynom<int> p1{ 4, 1 };
	const polynom<int> p2{ 1, -3, 5, 7 };
	const polynom<int> p3{ 11, -33, 55, 77 };
	EXPECT_EQ((polynom<int>{ 5, -2, 5, 7 }), p1 + p2);
	EXPECT_EQ((polynom<int>{ 3, 4, -5, -7 }), p1 - p2);
	EXPECT_EQ((polynom<int>{ -1, 3, -5, -7 }), -p2);
	EXPECT_EQ((polynom<int>{ 4, -11, 17, 33, 7 }), p1 * p2);
	EXPECT_EQ((polynom<int>{ }), p1 / p2);
	EXPECT_EQ((polynom<int>{ 4, 1 }), p1 % p2);
	EXPECT_EQ((polynom<int>{ 5, -2, 5, 7 }), p2 + p1);
	EXPECT_EQ((polynom<int>{ -3, -4, 5, 7 }), p2 - p1);
	EXPECT_EQ((polynom<int>{ -4, -1 }), -p1);
	EXPECT_EQ((polynom<int>{ 4, -11, 17, 33, 7 }), p2 * p1);
	EXPECT_EQ((polynom<int>{ 89, -23, 7 }), p2 / p1);
	EXPECT_EQ((polynom<int>{ -355 }), p2 % p1);
	EXPECT_EQ((polynom<int>{ 11, -33, 55, 77 }), p2 * 11);
	EXPECT_EQ((polynom<int>{ 1, -3, 5, 7 }), p3 / 11);
}

TEST(polynom_test, operators_inplace) {
	const polynom<int> p1{ 4, 1 };
	const polynom<int> p2{ 1, -3, 5, 7 };
	const polynom<int> p3{ 11, -33, 55, 77 };
	polynom<int> pr;
	pr = p1; pr += p2;
	EXPECT_EQ((polynom<int>{ 5, -2, 5, 7 }), pr);
	pr = p1; pr -= p2;
	EXPECT_EQ((polynom<int>{ 3, 4, -5, -7 }), pr);
	pr = p1; pr *= p2;
	EXPECT_EQ((polynom<int>{ 4, -11, 17, 33, 7 }), pr);
	pr = p1; pr /= p2;
	EXPECT_EQ((polynom<int>{ }), pr);
	pr = p1; pr %= p2;
	EXPECT_EQ((polynom<int>{ 4, 1 }), pr);
	pr = p2; pr += p1;
	EXPECT_EQ((polynom<int>{ 5, -2, 5, 7 }), pr);
	pr = p2; pr -= p1;
	EXPECT_EQ((polynom<int>{ -3, -4, 5, 7 }), pr);
	pr = p2; pr *= p1;
	EXPECT_EQ((polynom<int>{ 4, -11, 17, 33, 7 }), pr);
	pr = p2; pr /= p1;
	EXPECT_EQ((polynom<int>{ 89, -23, 7 }), pr);
	pr = p2; pr %= p1;
	EXPECT_EQ((polynom<int>{ -355 }), pr);
	pr = p2; pr *= 11;
	EXPECT_EQ((polynom<int>{ 11, -33, 55, 77 }), pr);
	pr = p3; pr /= 11;
	EXPECT_EQ((polynom<int>{ 1, -3, 5, 7 }), pr);
}

TEST(polynom_test, operators_inplace_self) {
	const polynom<int> p1{ 2, -3, 5, 7 };
	polynom<int> pr;
	pr = p1; pr += pr;
	EXPECT_EQ((polynom<int>{ 4, -6, 10, 14 }), pr);
	pr = p1; pr -= pr;
	EXPECT_EQ((polynom<int>{ 0, 0, 0, 0 }), pr);
	pr = p1; pr *= pr;
	EXPECT_EQ((polynom<int>{ 4, -12, 29, -2, -17, 70, 49 }), pr);
	pr = p1; pr %= pr;
	EXPECT_EQ((polynom<int>{ 0, 0, 0 }), pr);
	pr = p1; pr /= pr;
	EXPECT_EQ((polynom<int>{ 1 }), pr);
}

TEST(polynom_test, eval) {
	const polynom<int> p1{ 7, 5, -3, 2 };
	vector<int> ve1;
	for (int x = -3; x <= 4; x++) {
		ve1.push_back(p1.eval(x));
	}
	EXPECT_EQ((vector<int>{-89, -31, -3, 7, 11, 21, 49, 107}), ve1);
	vector<double> ve2;
	for (int x = -3; x <= 4; x++) {
		ve2.push_back(p1.eval(x + 0.5));
	}
	EXPECT_EQ((vector<double>{-55.5, -14, 3.5, 9, 14.5, 32, 73.5, 151}), ve2);
}

TEST(polynom_test, derivative) {
	const polynom<int> p1{ 7, 5, -3, 4 };
	const polynom<int> pd = p1.derivative();
	EXPECT_EQ((polynom<int>{ 5, -6, 12 }), pd);
}

TEST(polynom_test, search) {
	const polynom<int> p1{ 7, -5, -13, 4 };
	EXPECT_DOUBLE_EQ(-0.8262501959871101, p1.search(-1e100, -0.177, 0.0));
	EXPECT_DOUBLE_EQ(+0.6112574125565371, p1.search(-0.177, +2.344, 0.0));
	EXPECT_DOUBLE_EQ(+3.4649927834305730, p1.search(+2.344, +1e100, 0.0));
}

TEST(polynom_test, zeros) {
	const polynom<int> p1{ 7, -5, -13, 4 };
	vector<double> vz = p1.find_zeros(1e100);
	for (auto &z : vz) { z = round(z * 1e12) / 1e12; }
	EXPECT_EQ((vector<double>{-0.826250195987, 0.611257412557, 3.464992783431}), vz);
}
