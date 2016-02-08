#include "algorithm/math/base.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

TEST(base_test, absT) {
	EXPECT_EQ(0, absT(0));
	EXPECT_EQ(10, absT(10));
	EXPECT_EQ(10, absT(-10));
	EXPECT_EQ(0.0, absT(0.0));
	EXPECT_EQ(10.0, absT(10.0));
	EXPECT_EQ(10.0, absT(-10.0));
}

TEST(base_test, powT) {
	EXPECT_EQ(1, powT(0, 0));
	EXPECT_EQ(1, powT(+1, 0));
	EXPECT_EQ(1, powT(-1, 0));
	EXPECT_EQ(1, powT(+2, 0));
	EXPECT_EQ(1, powT(-2, 0));
	EXPECT_EQ(0, powT(0, 1));
	EXPECT_EQ(+1, powT(+1, 1));
	EXPECT_EQ(-1, powT(-1, 1));
	EXPECT_EQ(+2, powT(+2, 1));
	EXPECT_EQ(-2, powT(-2, 1));
	EXPECT_EQ(0, powT(0, 2));
	EXPECT_EQ(1, powT(+1, 2));
	EXPECT_EQ(1, powT(-1, 2));
	EXPECT_EQ(4, powT(+2, 2));
	EXPECT_EQ(4, powT(-2, 2));
	EXPECT_EQ(0, powT(0, 3));
	EXPECT_EQ(+1, powT(+1, 3));
	EXPECT_EQ(-1, powT(-1, 3));
	EXPECT_EQ(+8, powT(+2, 3));
	EXPECT_EQ(-8, powT(-2, 3));
	EXPECT_EQ(+8000000000000000000LL, powT(+2000000LL, 3));
	EXPECT_EQ(-8000000000000000000LL, powT(-2000000LL, 3));
	EXPECT_EQ(+8.0f, powT(+2.0f, 3));
	EXPECT_EQ(-8.0f, powT(-2.0f, 3));
	EXPECT_EQ(+8.0, powT(+2.0, 3));
	EXPECT_EQ(-8.0, powT(-2.0, 3));
}

TEST(base_test, sqT) {
	EXPECT_EQ(0, sqT(0));
	EXPECT_EQ(1, sqT(+1));
	EXPECT_EQ(1, sqT(-1));
	EXPECT_EQ(4, sqT(+2));
	EXPECT_EQ(4, sqT(-2));
	EXPECT_EQ(400000000, sqT(+20000));
	EXPECT_EQ(400000000, sqT(-20000));
	EXPECT_EQ(4000000000000000000LL, sqT(+2000000000LL));
	EXPECT_EQ(4000000000000000000LL, sqT(-2000000000LL));
	EXPECT_EQ(4.0f, sqT(+2.0f));
	EXPECT_EQ(4.0f, sqT(-2.0f));
	EXPECT_EQ(4.0, sqT(+2.0));
	EXPECT_EQ(4.0, sqT(-2.0));
}

TEST(base_test, isq) {
	EXPECT_EQ(0LL, isq(0));
	EXPECT_EQ(1LL, isq(+1));
	EXPECT_EQ(1LL, isq(-1));
	EXPECT_EQ(4LL, isq(+2));
	EXPECT_EQ(4LL, isq(-2));
	EXPECT_EQ(4000000000000000000LL, isq(+2000000000));
	EXPECT_EQ(4000000000000000000LL, isq(-2000000000));
}

TEST(base_test, isqrt) {
	EXPECT_EQ(0, isqrt(0));
	EXPECT_EQ(+1, isqrt(+1));
	EXPECT_EQ(-1, isqrt(-1));
	EXPECT_EQ(+1, isqrt(+2));
	EXPECT_EQ(-2, isqrt(-2));
	EXPECT_EQ(+3, isqrt(+9));
	EXPECT_EQ(-3, isqrt(-9));
	EXPECT_EQ(+3, isqrt(+10));
	EXPECT_EQ(-4, isqrt(-10));
	EXPECT_EQ(+3, isqrt(+15));
	EXPECT_EQ(-4, isqrt(-15));
	EXPECT_EQ(+4, isqrt(+16));
	EXPECT_EQ(-4, isqrt(-16));
	EXPECT_EQ(+1414213562, isqrt(+2000000000000000000LL));
	EXPECT_EQ(-1414213563, isqrt(-2000000000000000000LL));
	EXPECT_EQ(+2000000000, isqrt(+4000000000000000000LL));
	EXPECT_EQ(-2000000000, isqrt(-4000000000000000000LL));
}

TEST(base_test, isqrtc) {
	EXPECT_EQ(0, isqrtc(0));
	EXPECT_EQ(+1, isqrtc(+1));
	EXPECT_EQ(-1, isqrtc(-1));
	EXPECT_EQ(+2, isqrtc(+2));
	EXPECT_EQ(-1, isqrtc(-2));
	EXPECT_EQ(+3, isqrtc(+9));
	EXPECT_EQ(-3, isqrtc(-9));
	EXPECT_EQ(+4, isqrtc(+10));
	EXPECT_EQ(-3, isqrtc(-10));
	EXPECT_EQ(+4, isqrtc(+15));
	EXPECT_EQ(-3, isqrtc(-15));
	EXPECT_EQ(+4, isqrtc(+16));
	EXPECT_EQ(-4, isqrtc(-16));
	EXPECT_EQ(+1414213563, isqrtc(+2000000000000000000LL));
	EXPECT_EQ(-1414213562, isqrtc(-2000000000000000000LL));
	EXPECT_EQ(+2000000000, isqrtc(+4000000000000000000LL));
	EXPECT_EQ(-2000000000, isqrtc(-4000000000000000000LL));
}

TEST(base_test, is_square) {
	EXPECT_EQ(true, is_square(0));
	EXPECT_EQ(true, is_square(+1));
	EXPECT_EQ(false, is_square(-1));
	EXPECT_EQ(false, is_square(+2));
	EXPECT_EQ(false, is_square(-2));
	EXPECT_EQ(true, is_square(+9));
	EXPECT_EQ(false, is_square(-9));
	EXPECT_EQ(false, is_square(+2000000000000000000LL));
	EXPECT_EQ(false, is_square(-2000000000000000000LL));
	EXPECT_EQ(true, is_square(+4000000000000000000LL));
	EXPECT_EQ(false, is_square(-4000000000000000000LL));
}

TEST(base_test, icb) {
	EXPECT_EQ(0LL, icb(0));
	EXPECT_EQ(+1LL, icb(+1));
	EXPECT_EQ(-1LL, icb(-1));
	EXPECT_EQ(+8LL, icb(+2));
	EXPECT_EQ(-8LL, icb(-2));
	EXPECT_EQ(+8000000000000000000LL, icb(+2000000));
	EXPECT_EQ(-8000000000000000000LL, icb(-2000000));
}

TEST(base_test, icbrt) {
	EXPECT_EQ(0, icbrt(0));
	EXPECT_EQ(+1, icbrt(+1));
	EXPECT_EQ(-1, icbrt(-1));
	EXPECT_EQ(+1, icbrt(+2));
	EXPECT_EQ(-2, icbrt(-2));
	EXPECT_EQ(+3, icbrt(+27));
	EXPECT_EQ(-3, icbrt(-27));
	EXPECT_EQ(+3, icbrt(+28));
	EXPECT_EQ(-4, icbrt(-28));
	EXPECT_EQ(+3, icbrt(+63));
	EXPECT_EQ(-4, icbrt(-63));
	EXPECT_EQ(+4, icbrt(+64));
	EXPECT_EQ(-4, icbrt(-64));
	EXPECT_EQ(+1259921, icbrt(+2000000000000000000LL));
	EXPECT_EQ(-1259922, icbrt(-2000000000000000000LL));
	EXPECT_EQ(+2000000, icbrt(+8000000000000000000LL));
	EXPECT_EQ(-2000000, icbrt(-8000000000000000000LL));
}

TEST(base_test, icbrtc) {
	EXPECT_EQ(0, icbrtc(0));
	EXPECT_EQ(+1, icbrtc(+1));
	EXPECT_EQ(-1, icbrtc(-1));
	EXPECT_EQ(+2, icbrtc(+2));
	EXPECT_EQ(-1, icbrtc(-2));
	EXPECT_EQ(+3, icbrtc(+27));
	EXPECT_EQ(-3, icbrtc(-27));
	EXPECT_EQ(+4, icbrtc(+28));
	EXPECT_EQ(-3, icbrtc(-28));
	EXPECT_EQ(+4, icbrtc(+63));
	EXPECT_EQ(-3, icbrtc(-63));
	EXPECT_EQ(+4, icbrtc(+64));
	EXPECT_EQ(-4, icbrtc(-64));
	EXPECT_EQ(+1259922, icbrtc(+2000000000000000000LL));
	EXPECT_EQ(-1259921, icbrtc(-2000000000000000000LL));
	EXPECT_EQ(+2000000, icbrtc(+8000000000000000000LL));
	EXPECT_EQ(-2000000, icbrtc(-8000000000000000000LL));
}

TEST(base_test, is_cube) {
	EXPECT_EQ(true, is_cube(0));
	EXPECT_EQ(true, is_cube(+1));
	EXPECT_EQ(true, is_cube(-1));
	EXPECT_EQ(false, is_cube(+2));
	EXPECT_EQ(false, is_cube(-2));
	EXPECT_EQ(true, is_cube(+27));
	EXPECT_EQ(true, is_cube(-27));
	EXPECT_EQ(false, is_cube(+2000000000000000000LL));
	EXPECT_EQ(false, is_cube(-2000000000000000000LL));
	EXPECT_EQ(true, is_cube(+8000000000000000000LL));
	EXPECT_EQ(true, is_cube(-8000000000000000000LL));
}

TEST(base_test, div_floor) {
	EXPECT_EQ(+6LL, div_floor(+20, +3));
	EXPECT_EQ(-7LL, div_floor(+20, -3));
	EXPECT_EQ(-7LL, div_floor(-20, +3));
	EXPECT_EQ(+6LL, div_floor(-20, -3));
}

TEST(base_test, div_ceil) {
	EXPECT_EQ(+7LL, div_ceil(+20, +3));
	EXPECT_EQ(-6LL, div_ceil(+20, -3));
	EXPECT_EQ(-6LL, div_ceil(-20, +3));
	EXPECT_EQ(+7LL, div_ceil(-20, -3));
}

TEST(base_test, gcd) {
	EXPECT_EQ(15, gcd(15, 0));
	EXPECT_EQ(15, gcd(0, 15));
	EXPECT_EQ(15, gcd(15 * 7, 15 * 11 * 13));
	EXPECT_EQ(15, gcd(15 * 11 * 13, 15 * 7));
	EXPECT_EQ(1, gcd(7, 11 * 13));
	EXPECT_EQ(1, gcd(11 * 13, 7));
}

void gcd_ex_test_impl(int expected_g, int a, int b) {
	int g, x, y;
	g = gcd_ex(a, b, &x, &y);
	EXPECT_EQ(expected_g, g);
	EXPECT_EQ(expected_g, a * x + b * y);
}

TEST(base_test, gcd_ex) {
	gcd_ex_test_impl(15, 15, 0);
	gcd_ex_test_impl(15, 0, 15);
	gcd_ex_test_impl(15, 15 * 7, 15 * 13);
	gcd_ex_test_impl(15, 15 * 13, 15 * 7);
	gcd_ex_test_impl(1, 7, 11 * 13);
	gcd_ex_test_impl(1, 11 * 13, 7);
}

TEST(base_test, lcm) {
	EXPECT_EQ(0, lcm(15, 0));
	EXPECT_EQ(0, lcm(0, 15));
	EXPECT_EQ(15, lcm(15, 15));
	EXPECT_EQ(30, lcm(15, 10));
	EXPECT_EQ(210, lcm(14, 15));
}

void crt_test_impl(int a1, int n1, int a2, int n2) {
	int a, n; chinese_remainder(a, n, a1, n1, a2, n2);
	EXPECT_EQ(lcm(n1, n2), n);
	EXPECT_EQ(a1, a % n1);
	EXPECT_EQ(a2, a % n2);
	EXPECT_GE(a, 0);
	EXPECT_LT(a, n);
}

TEST(base_test, crt) {
	crt_test_impl(0, 10, 5, 13);
	crt_test_impl(5, 10, 3, 13);
	crt_test_impl(5, 10, 3, 14);
	crt_test_impl(4, 10, 6, 14);
	crt_test_impl(6, 14, 6, 14);
	crt_test_impl(102, 65535, 12345, 55557);
}
