#include "altruct/structure/math/quadratic.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

typedef quadratic<int, 5> quad;
typedef quadratic<int, -1> gaussian;
typedef quadraticX<int> quadx;

TEST(quadratic_test, constructor) {
	quad q1;
	EXPECT_EQ(0, q1.a);
	EXPECT_EQ(0, q1.b);
	EXPECT_EQ(5, q1.D());
	quad q2(10);
	EXPECT_EQ(10, q2.a);
	EXPECT_EQ(0, q2.b);
	EXPECT_EQ(5, q2.D());
	quad q3(+2, -5);
	EXPECT_EQ(+2, q3.a);
	EXPECT_EQ(-5, q3.b);
	EXPECT_EQ(5, q3.D());
	quad q4(+2, -5, 7); // ignore D if static
	EXPECT_EQ(+2, q4.a);
	EXPECT_EQ(-5, q4.b);
	EXPECT_EQ(5, q4.D());
	quad q5(q4);
	EXPECT_EQ(+2, q5.a);
	EXPECT_EQ(-5, q5.b);
	EXPECT_EQ(5, q5.D());
}

TEST(quadratic_test, constructor_x) {
	quadx q1;
	EXPECT_EQ(0, q1.a);
	EXPECT_EQ(0, q1.b);
	EXPECT_EQ(0, q1.D());
	quadx q2(10);
	EXPECT_EQ(10, q2.a);
	EXPECT_EQ(0, q2.b);
	EXPECT_EQ(0, q2.D());
	quadx q3(+2, -5);
	EXPECT_EQ(+2, q3.a);
	EXPECT_EQ(-5, q3.b);
	EXPECT_EQ(0, q3.D());
	quadx q4(+2, -5, 7); // ignore D if static
	EXPECT_EQ(+2, q4.a);
	EXPECT_EQ(-5, q4.b);
	EXPECT_EQ(7, q4.D());
	quadx q5(q4);
	EXPECT_EQ(+2, q5.a);
	EXPECT_EQ(-5, q5.b);
	EXPECT_EQ(7, q5.D());
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

TEST(quadratic_test, operators_comparison) {
	test_comparison(true, false, quad(2, 5), quad(2, 5));
	test_comparison(false, false, quad(2, 5), quad(2, 3));
	test_comparison(false, true, quad(2, 5), quad(2, 7));
	test_comparison(false, true, quad(2, 5), quad(4, 5));
	test_comparison(false, true, quad(2, 5), quad(4, 3));
	test_comparison(false, true, quad(2, 5), quad(4, 7));
	test_comparison(false, false, quad(2, 5), quad(1, 5));
	test_comparison(false, false, quad(2, 5), quad(1, 3));
	test_comparison(false, false, quad(2, 5), quad(1, 7));
}

TEST(quadratic_test, operators_arithmetic) {
	const quad q1(2, -5);
	const quad q2(3, 4);
	const quad q3(3, -2);
	EXPECT_EQ(quad(5, -1), q1 + q2);
	EXPECT_EQ(quad(-1, -9), q1 - q2);
	EXPECT_EQ(quad(-2, 5), -q1);
	EXPECT_EQ(quad(-94, -7), q1 * q2);
	EXPECT_EQ(quad(4, 1), q1 / q3);
	EXPECT_EQ(quad(5, -1), q1 % q2);
	EXPECT_EQ(quad(5, -1), q2 + q1);
	EXPECT_EQ(quad(1, 9), q2 - q1);
	EXPECT_EQ(quad(-3, -4), -q2);
	EXPECT_EQ(quad(-94, -7), q2 * q1);
	EXPECT_EQ(quad(-6, 15), q1 * -3);
	EXPECT_EQ(quad(1, -2), q1 / 2);
}

TEST(quadratic_test, operators_inplace) {
	const quad q1(2, -5);
	const quad q2(3, 4);
	const quad q3(3, -2);
	quad qr;
	qr = q1; qr += q2;
	EXPECT_EQ(quad(5, -1), qr);
	qr = q1; qr -= q2;
	EXPECT_EQ(quad(-1, -9), qr);
	qr = q1; qr *= q2;
	EXPECT_EQ(quad(-94, -7), qr);
	qr = q1; qr /= q3;
	EXPECT_EQ(quad(4, 1), qr);
	qr = q1; qr %= q2;
	EXPECT_EQ(quad(5, -1), qr);
	qr = q2; qr += q1;
	EXPECT_EQ(quad(5, -1), qr);
	qr = q2; qr -= q1;
	EXPECT_EQ(quad(1, 9), qr);
	qr = q2; qr *= q1;
	EXPECT_EQ(quad(-94, -7), qr);
	qr = q1; qr *= -3;
	EXPECT_EQ(quad(-6, 15), qr);
	qr = q1; qr /= 2;
	EXPECT_EQ(quad(1, -2), qr);
}

TEST(quadratic_test, operators_inplace_self) {
	const quad q1(2, -5);
	quad qr;
	qr = q1; qr += qr;
	EXPECT_EQ(quad(4, -10), qr);
	qr = q1; qr -= qr;
	EXPECT_EQ(quad(0, 0), qr);
	qr = q1; qr *= qr;
	EXPECT_EQ(quad(129, -20), qr);
	qr = q1; qr /= qr;
	EXPECT_EQ(quad(1, 0), qr);
	qr = q1; qr %= qr;
	EXPECT_EQ(quad(0, 0), qr);
}

TEST(quadratic_test, conjugate) {
	const quad q1(2, -5);
	const quad q2(2, 3);
	EXPECT_EQ(quad(2, 5), q1.conjugate());
	EXPECT_EQ(quad(2, -3), q2.conjugate());
}

TEST(quadratic_test, norm) {
	const quad q1(2, -5);
	const quad q2(3, 4);
	EXPECT_EQ(-121, q1.norm());
	EXPECT_EQ(-71, q2.norm());
	const gaussian g1(2, -5);
	const gaussian g2(3, 4);
	EXPECT_EQ(29, g1.norm());
	EXPECT_EQ(25, g2.norm());
}

TEST(quadratic_test, casts) {
	const quad q(2, -5);
	const quad e0 = zeroT<quad>::of(q);
	const quad e1 = identityT<quad>::of(q);
	EXPECT_EQ(0, e0.a);
	EXPECT_EQ(0, e0.b);
	EXPECT_EQ(5, e0.D());
	EXPECT_EQ(1, e1.a);
	EXPECT_EQ(0, e1.b);
	EXPECT_EQ(5, e1.D());
    const quad q3 = castOf<quad>(3);
    EXPECT_EQ(3, q3.a);
    EXPECT_EQ(0, q3.b);
    EXPECT_EQ(5, q3.D());
    const quad q4 = castOf(q, 4);
    EXPECT_EQ(4, q4.a);
    EXPECT_EQ(0, q4.b);
    EXPECT_EQ(5, q4.D());
    const quad q6 = castOf(q, q4);
    EXPECT_EQ(4, q6.a);
    EXPECT_EQ(0, q6.b);
    EXPECT_EQ(5, q6.D());
    const quad q7 = castOf<quad>(q4);
    EXPECT_EQ(4, q7.a);
    EXPECT_EQ(0, q7.b);
    EXPECT_EQ(5, q7.D());
}

TEST(quadratic_test, casts_x) {
	const quadx z(2, -5, -1);
	const quadx z0 = zeroT<quadx>::of(z);
	const quadx z1 = identityT<quadx>::of(z);
	EXPECT_EQ(0, z0.a);
	EXPECT_EQ(0, z0.b);
	EXPECT_EQ(-1, z0.D());
	EXPECT_EQ(1, z1.a);
	EXPECT_EQ(0, z1.b);
	EXPECT_EQ(-1, z1.D());
    const quadx z5 = castOf(z, 5);
    EXPECT_EQ(5, z5.a);
    EXPECT_EQ(0, z5.b);
    EXPECT_EQ(-1, z5.D());
    const quadx z6 = castOf(z, z5);
    EXPECT_EQ(5, z6.a);
    EXPECT_EQ(0, z6.b);
    EXPECT_EQ(-1, z6.D());
    const quadx z7 = castOf<quadx>(z5);
    EXPECT_EQ(5, z7.a);
    EXPECT_EQ(0, z7.b);
    EXPECT_EQ(-1, z7.D());
}
