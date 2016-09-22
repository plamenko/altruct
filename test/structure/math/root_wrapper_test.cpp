#include "structure/math/root_wrapper.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

typedef complex<double> cplx;

TEST(root_wrapper_test, complex_root_wrapper) {
	auto w = complex_root_wrapper<double>(6);
	EXPECT_EQ(8, w.size);
	auto w0 = identityOf(w); cplx r0 = w0;
	EXPECT_NEAR(1.0, r0.a, 1e-9);
	EXPECT_NEAR(0.0, r0.b, 1e-9);
	auto w1 = w0 * w; cplx r1 = w1;
	EXPECT_NEAR(0.707106781186548, r1.a, 1e-9);
	EXPECT_NEAR(0.707106781186548, r1.b, 1e-9);
	auto w2 = w1; w2 *= w1; cplx r2 = w2;
	EXPECT_NEAR(0.0, r2.a, 1e-9);
	EXPECT_NEAR(1.0, r2.b, 1e-9);
	auto w3 = w1 * w1 * w1; cplx r3 = w3;
	EXPECT_NEAR(-0.707106781186548, r3.a, 1e-9);
	EXPECT_NEAR(+0.707106781186548, r3.b, 1e-9);
	auto w4 = w2; w4 *= w2; cplx r4 = w4;
	EXPECT_NEAR(-1.0, r4.a, 1e-9);
	EXPECT_NEAR( 0.0, r4.b, 1e-9);
	auto w5 = powT(w1, 5); cplx r5 = w5;
	EXPECT_NEAR(-0.707106781186548, r5.a, 1e-9);
	EXPECT_NEAR(-0.707106781186548, r5.b, 1e-9);
	auto w6 = w1 * w2 * w3; cplx r6 = w6;
	EXPECT_NEAR( 0.0, r6.a, 1e-9);
	EXPECT_NEAR(-1.0, r6.b, 1e-9);
	auto w7 = w5 * w2; cplx r7 = w7;
	EXPECT_NEAR(+0.707106781186548, r7.a, 1e-9);
	EXPECT_NEAR(-0.707106781186548, r7.b, 1e-9);
	auto w8 = w1 * w7; cplx r8 = w8;
	EXPECT_NEAR(1.0, r0.a, 1e-9);
	EXPECT_NEAR(0.0, r0.b, 1e-9);
	auto w9 = w6 * w3; cplx r9 = w9;
	EXPECT_NEAR(0.707106781186548, r9.a, 1e-9);
	EXPECT_NEAR(0.707106781186548, r9.b, 1e-9);
}
