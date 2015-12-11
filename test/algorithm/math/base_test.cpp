#include "algorithm/math/base.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

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
