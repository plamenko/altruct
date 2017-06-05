#include "altruct/structure/math/complex.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

typedef complex<double> cplx;

TEST(complex_test, constructor) {
	cplx z1;
	EXPECT_EQ(0.0, z1.a);
	EXPECT_EQ(0.0, z1.b);
	cplx z2(5.0);
	EXPECT_EQ(5.0, z2.a);
	EXPECT_EQ(0.0, z2.b);
	cplx z3(+2, -5);
	EXPECT_EQ(+2.0, z3.a);
	EXPECT_EQ(-5.0, z3.b);
	cplx z4(z3);
	EXPECT_EQ(+2.0, z4.a);
	EXPECT_EQ(-5.0, z4.b);
}

template<typename T>
void test_comparison(bool eq, bool lt, const T& lhs, const T& rhs) {
	ASSERT_FALSE(eq && lt);
	EXPECT_EQ(eq, lhs == rhs);
	EXPECT_EQ(!eq, lhs != rhs);
	EXPECT_EQ(lt, lhs < rhs);
	EXPECT_EQ(!(lt || eq), lhs > rhs);
	EXPECT_EQ((lt || eq), lhs <= rhs);
	EXPECT_EQ(!lt, lhs >= rhs);
}

TEST(complex_test, operators_comparison) {
	test_comparison(true, false, cplx(2, 5), cplx(2, 5));
	test_comparison(false, false, cplx(2, 5), cplx(2, 3));
	test_comparison(false, true, cplx(2, 5), cplx(2, 7));
	test_comparison(false, true, cplx(2, 5), cplx(4, 5));
	test_comparison(false, true, cplx(2, 5), cplx(4, 3));
	test_comparison(false, true, cplx(2, 5), cplx(4, 7));
	test_comparison(false, false, cplx(2, 5), cplx(1, 5));
	test_comparison(false, false, cplx(2, 5), cplx(1, 3));
	test_comparison(false, false, cplx(2, 5), cplx(1, 7));
}

TEST(complex_test, operators_arithmetic) {
	const cplx z1(2.0, -5.0);
	const cplx z2(3.0, +4.0);
	const cplx z3(3.0, -2.0);
	EXPECT_EQ(cplx(5, -1), z1 + z2);
	EXPECT_EQ(cplx(-1, -9), z1 - z2);
	EXPECT_EQ(cplx(-2, 5), -z1);
	EXPECT_EQ(cplx(26, -7), z1 * z2);
	EXPECT_EQ(cplx(16, -11) / 13.0, z1 / z3);
	EXPECT_EQ(cplx(0, 0), z1 % z2);
	EXPECT_EQ(cplx(5, -1), z2 + z1);
	EXPECT_EQ(cplx(1, 9), z2 - z1);
	EXPECT_EQ(cplx(-3, -4), -z2);
	EXPECT_EQ(cplx(26, -7), z2 * z1);
	EXPECT_EQ(cplx(-7.0, 17.5), z1 * -3.5);
	EXPECT_EQ(cplx(1.0, -2.5), z1 / 2);
}

TEST(complex_test, operators_inplace) {
	const cplx z1(2, -5);
	const cplx z2(3, 4);
	const cplx z3(3, -2);
	cplx zr;
	zr = z1; zr += z2;
	EXPECT_EQ(cplx(5, -1), zr);
	zr = z1; zr -= z2;
	EXPECT_EQ(cplx(-1, -9), zr);
	zr = z1; zr *= z2;
	EXPECT_EQ(cplx(26, -7), zr);
	zr = z1; zr /= z3;
	EXPECT_EQ(cplx(16, -11) / 13.0, zr);
	zr = z1; zr %= z2;
	EXPECT_EQ(cplx(0, 0), zr);
	zr = z2; zr += z1;
	EXPECT_EQ(cplx(5, -1), zr);
	zr = z2; zr -= z1;
	EXPECT_EQ(cplx(1, 9), zr);
	zr = z2; zr *= z1;
	EXPECT_EQ(cplx(26, -7), zr);
	zr = z1; zr *= -3.5;
	EXPECT_EQ(cplx(-7, 17.5), zr);
	zr = z1; zr /= 2;
	EXPECT_EQ(cplx(1, -2.5), zr);
}

TEST(complex_test, operators_inplace_self) {
	const cplx z1(2, -5);
	cplx zr;
	zr = z1; zr += zr;
	EXPECT_EQ(cplx(4, -10), zr);
	zr = z1; zr -= zr;
	EXPECT_EQ(cplx(0, 0), zr);
	zr = z1; zr *= zr;
	EXPECT_EQ(cplx(-21, -20), zr);
	zr = z1; zr /= zr;
	EXPECT_EQ(cplx(1, 0), zr);
	zr = z1; zr %= zr;
	EXPECT_EQ(cplx(0, 0), zr);
}

TEST(complex_test, conjugate) {
	const cplx z1(2, -5);
	const cplx z2(2, 3);
	EXPECT_EQ(cplx(2, 5), z1.conjugate());
	EXPECT_EQ(cplx(2, -3), z2.conjugate());
}

TEST(complex_test, norm) {
	const cplx z1(2, -5);
	const cplx z2(3, 4);
	EXPECT_EQ(29.0, z1.norm());
	EXPECT_EQ(25.0, z2.norm());
}

TEST(complex_test, identity) {
	const cplx z1(2, -5);
	const cplx e0 = zeroT<cplx>::of(z1);
	const cplx e1 = identityT<cplx>::of(z1);
	EXPECT_EQ(0, e0.a);
	EXPECT_EQ(0, e0.b);
	EXPECT_EQ(-1, e0.D());
	EXPECT_EQ(1, e1.a);
	EXPECT_EQ(0, e1.b);
	EXPECT_EQ(-1, e1.D());
}
