#include "algorithm/math/base.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

const int num_iterations = 3000;
#define EXPECT_NEAR2(expected, actual, eps) EXPECT_TRUE(absT((actual) - (expected)) <= absT(eps))

template<typename T>
class wrapped {
public:
	T v;
	wrapped(const T& v = 0) : v(v) {}

	bool operator == (const wrapped& rhs) const { return v == rhs.v; }
	bool operator != (const wrapped& rhs) const { return v != rhs.v; }
	bool operator <  (const wrapped& rhs) const { return v < rhs.v; }
	bool operator >  (const wrapped& rhs) const { return v > rhs.v; }
	bool operator <= (const wrapped& rhs) const { return v <= rhs.v; }
	bool operator >= (const wrapped& rhs) const { return v >= rhs.v; }

	wrapped  operator +  (const wrapped &rhs) const { return wrapped(v + rhs.v); }
	wrapped  operator -  (const wrapped &rhs) const { return wrapped(v - rhs.v); }
	wrapped  operator +  ()                   const { return wrapped(+v); }
	wrapped  operator -  ()                   const { return wrapped(-v); }
	wrapped  operator *  (const wrapped &rhs) const { return wrapped(v * rhs.v); }
	wrapped  operator /  (const wrapped &rhs) const { return wrapped(v / rhs.v); }
	wrapped  operator %  (const wrapped &rhs) const { return wrapped(v % rhs.v); }

	wrapped& operator += (const wrapped &rhs) { v += rhs.v; return *this; }
	wrapped& operator -= (const wrapped &rhs) { v -= rhs.v; return *this; }
	wrapped& operator *= (const wrapped &rhs) { v *= rhs.v; return *this; }
	wrapped& operator /= (const wrapped &rhs) { v /= rhs.v; return *this; }
	wrapped& operator %= (const wrapped &rhs) { v %= rhs.v; return *this; }
};

TEST(base_test, cast) {
    EXPECT_EQ(123, (castT<int, int64_t>::of(123LL)));
    EXPECT_EQ(123.0, (castT<double, int64_t>::of(123LL)));
    EXPECT_EQ(123LL, (castT<int64_t, int>::of(123)));
    EXPECT_EQ(123LL, (castT<int64_t, double>::of(123.5)));
    EXPECT_EQ(123, (castOf<int>(123LL)));
    EXPECT_EQ(123.0, (castOf<double>(123LL)));
    EXPECT_EQ(123LL, (castOf<int64_t>(123)));
    EXPECT_EQ(123LL, (castOf<int64_t>(123.5)));
    EXPECT_EQ(123LL, (castOf(7LL, 123.5)));
}

TEST(base_test, identity) {
	EXPECT_EQ(1, identityT<int>::of(0));
	EXPECT_EQ(1, identityT<int>::of(1));
	EXPECT_EQ(1, identityT<int>::of(5));
	EXPECT_EQ(1.0, identityT<double>::of(0.0));
	EXPECT_EQ(1.0, identityT<double>::of(1.0));
	EXPECT_EQ(1.0, identityT<double>::of(5.0));
    EXPECT_EQ(1, identityOf(0));
    EXPECT_EQ(1, identityOf(1));
    EXPECT_EQ(1, identityOf(5));
    EXPECT_EQ(1.0, identityOf(0.0));
    EXPECT_EQ(1.0, identityOf(1.0));
    EXPECT_EQ(1.0, identityOf(5.0));
}

TEST(base_test, zero) {
	EXPECT_EQ(0, zeroT<int>::of(0));
	EXPECT_EQ(0, zeroT<int>::of(1));
	EXPECT_EQ(0, zeroT<int>::of(5));
	EXPECT_EQ(0.0, zeroT<double>::of(0.0));
	EXPECT_EQ(0.0, zeroT<double>::of(1.0));
	EXPECT_EQ(0.0, zeroT<double>::of(5.0));
    EXPECT_EQ(0, zeroOf(0));
    EXPECT_EQ(0, zeroOf(1));
    EXPECT_EQ(0, zeroOf(5));
    EXPECT_EQ(0.0, zeroOf(0.0));
    EXPECT_EQ(0.0, zeroOf(1.0));
    EXPECT_EQ(0.0, zeroOf(5.0));
}

TEST(base_test, absT) {
	EXPECT_EQ(0, absT(0));
	EXPECT_EQ(10, absT(10));
	EXPECT_EQ(10, absT(-10));
	EXPECT_EQ(0.0, absT(0.0));
	EXPECT_EQ(10.0, absT(10.0));
	EXPECT_EQ(10.0, absT(-10.0));
}

TEST(base_test, minT) {
	EXPECT_EQ(2, minT(2, 5));
	EXPECT_EQ(-5, minT(2, -5));
	EXPECT_EQ(-2, minT(-2, 5));
	EXPECT_EQ(-5, minT(-2, -5));
	EXPECT_EQ(2, minT(5, 2));
	EXPECT_EQ(-5, minT(-5, 2));
	EXPECT_EQ(-2, minT(5, -2));
	EXPECT_EQ(-5, minT(-5, -2));
	EXPECT_EQ(2.0, minT(2.0, 5.0));
	EXPECT_EQ(-5.0, minT(2.0, -5.0));
	EXPECT_EQ(-2.0, minT(-2.0, 5.0));
	EXPECT_EQ(-5.0, minT(-2.0, -5.0));
}

TEST(base_test, maxT) {
	EXPECT_EQ(5, maxT(2, 5));
	EXPECT_EQ(2, maxT(2, -5));
	EXPECT_EQ(5, maxT(-2, 5));
	EXPECT_EQ(-2, maxT(-2, -5));
	EXPECT_EQ(5, maxT(5, 2));
	EXPECT_EQ(2, maxT(-5, 2));
	EXPECT_EQ(5, maxT(5, -2));
	EXPECT_EQ(-2, maxT(-5, -2));
	EXPECT_EQ(5.0, maxT(2.0, 5.0));
	EXPECT_EQ(2.0, maxT(2.0, -5.0));
	EXPECT_EQ(5.0, maxT(-2.0, 5.0));
	EXPECT_EQ(-2.0, maxT(-2.0, -5.0));
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
	int x = 1;
	for (int i = 0; i < 20; i++) {
		EXPECT_EQ(x, powT(3, i)) << "powT(3, " << i << ")";
		x *= 3;
	}
	double y = 1;
	for (int i = 0; i < 20; i++) {
		EXPECT_NEAR(y, powT(2.5, i), 1e-6) << "powT(2.5, " << i << ")";
		y *= 2.5;
	}
}

template<typename I>
void test_sqT_uint() {
	EXPECT_EQ(I(0), sqT<I>(0));
	EXPECT_EQ(I(1), sqT<I>(+1));
	EXPECT_EQ(I(4), sqT<I>(+2));
	EXPECT_EQ(I(100), sqT<I>(+10));
	for (int i = 0; i < num_iterations; i++) {
		EXPECT_EQ(I(i * i), sqT(I(i))) << "sqT(" << i << ")";
	}
}

template<typename I>
void test_sqT_int() {
	EXPECT_EQ(I(1), sqT<I>(-1));
	EXPECT_EQ(I(4), sqT<I>(-2));
	EXPECT_EQ(I(100), sqT<I>(-10));
	for (int i = -num_iterations; i < num_iterations; i++) {
		EXPECT_EQ(I(i * i), sqT(I(-i))) << "sqT(" << i << ")";
	}
}

template<typename F>
void test_sqT_float() {
	float eps = 1e-6f;
	EXPECT_EQ(F(6.25), sqT<F>(+2.5));
	EXPECT_EQ(F(6.25), sqT<F>(-2.5));
	for (int i = -num_iterations; i < num_iterations; i++) {
		F epsi = eps * (abs(i) + 1);
		EXPECT_NEAR2(F(i + 0.0f) * F(i + 0.0f), sqT(F(i + 0.0f)), epsi) << "sqT(" << (i + 0.0f) << ")";
		EXPECT_NEAR2(F(i + 0.5f) * F(i + 0.5f), sqT(F(i + 0.5f)), epsi) << "sqT(" << (i + 0.5f) << ")";
	}
}

TEST(base_test, sqT) {
	test_sqT_int<int8_t>();
	test_sqT_int<int16_t>();
	test_sqT_int<int32_t>();
	test_sqT_int<int64_t>();
	test_sqT_int<wrapped<int>>();
	test_sqT_uint<uint8_t>();
	test_sqT_uint<uint16_t>();
	test_sqT_uint<uint32_t>();
	test_sqT_uint<uint64_t>();
	test_sqT_float<float>();
	test_sqT_float<double>();
	test_sqT_float<wrapped<float>>();
	EXPECT_EQ(400000000, sqT(+20000));
	EXPECT_EQ(400000000, sqT(-20000));
	EXPECT_EQ(4000000000000000000LL, sqT(+2000000000LL));
	EXPECT_EQ(4000000000000000000LL, sqT(-2000000000LL));
}

template<typename I>
void test_sqrtT_uint(int max_val = num_iterations) {
	EXPECT_EQ(I(0), sqrtT<I>(0));
	EXPECT_EQ(I(+1), sqrtT<I>(+1));
	EXPECT_EQ(I(+1), sqrtT<I>(+2));
	EXPECT_EQ(I(+1), sqrtT<I>(+3));
	EXPECT_EQ(I(+2), sqrtT<I>(+4));
	EXPECT_EQ(I(+3), sqrtT<I>(+9));
	EXPECT_EQ(I(+3), sqrtT<I>(+10));
	EXPECT_EQ(I(+3), sqrtT<I>(+15));
	EXPECT_EQ(I(+4), sqrtT<I>(+16));
	for (int i = 0; i <= max_val; i++) {
		EXPECT_EQ(I(int(floor(sqrt(double(i))))), sqrtT(I(i))) << "sqrtT(" << i << ")";
	}
}

template<typename I>
void test_sqrtT_int(int max_val = num_iterations) {
	EXPECT_EQ(I(-1), sqrtT<I>(-1));
	EXPECT_EQ(I(-1), sqrtT<I>(-2));
	EXPECT_EQ(I(-1), sqrtT<I>(-3));
	EXPECT_EQ(I(-2), sqrtT<I>(-4));
	EXPECT_EQ(I(-3), sqrtT<I>(-9));
	EXPECT_EQ(I(-3), sqrtT<I>(-10));
	EXPECT_EQ(I(-3), sqrtT<I>(-15));
	EXPECT_EQ(I(-4), sqrtT<I>(-16));
	for (int i = 0; i <= max_val; i++) {
		EXPECT_EQ(I(+int(floor(sqrt(double(i))))), sqrtT(I(+i))) << "sqrtT(" << (+i) << ")";
		EXPECT_EQ(I(-int(floor(sqrt(double(i))))), sqrtT(I(-i))) << "sqrtT(" << (-i) << ")";
	}
}

template<typename F>
void test_sqrtT_float(int max_val = num_iterations) {
	F eps(1e-6f);
	EXPECT_EQ(F(0), sqrtT<F>(0, eps));
	EXPECT_EQ(F(1), sqrtT<F>(1, eps));
	EXPECT_NEAR2(F(1.414213562373095), sqrtT(F(2), eps), eps);
	EXPECT_NEAR2(F(1.732050807568877), sqrtT(F(3), eps), eps);
	EXPECT_EQ(F(2), sqrtT<F>(4, eps));
	EXPECT_NEAR2(F(2.5), sqrtT(F(6.25), eps), eps);
	for (int i = 0; i <= max_val; i++) {
		F epsi = eps * (abs(i) + 1);
		EXPECT_NEAR2(F(+sqrt(i + 0.0)), sqrtT(F(+(i + 0.0)), eps), epsi) << "sqrtT(" << (+(i + 0.0)) << ")";
		EXPECT_NEAR2(F(+sqrt(i + 0.5)), sqrtT(F(+(i + 0.5)), eps), epsi) << "sqrtT(" << (+(i + 0.5)) << ")";
	}
}

TEST(base_test, sqrtT) {
	test_sqrtT_int<int8_t>(127);
	test_sqrtT_int<int16_t>();
	test_sqrtT_int<int32_t>();
	test_sqrtT_int<int64_t>();
	test_sqrtT_int<wrapped<int>>();
	test_sqrtT_uint<uint8_t>(127);
	test_sqrtT_uint<uint16_t>();
	test_sqrtT_uint<uint32_t>();
	test_sqrtT_uint<uint64_t>();
	test_sqrtT_float<float>();
	test_sqrtT_float<double>();
	test_sqrtT_float<wrapped<double>>();
	EXPECT_EQ(+1414213562LL, sqrtT(+2000000000000000000LL));
	EXPECT_EQ(-1414213562LL, sqrtT(-2000000000000000000LL));
	EXPECT_EQ(+2000000000LL, sqrtT(+4000000000000000000LL));
	EXPECT_EQ(-2000000000LL, sqrtT(-4000000000000000000LL));
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
	EXPECT_EQ(-1, isqrt(-2));
	EXPECT_EQ(+3, isqrt(+9));
	EXPECT_EQ(-3, isqrt(-9));
	EXPECT_EQ(+3, isqrt(+10));
	EXPECT_EQ(-3, isqrt(-10));
	EXPECT_EQ(+3, isqrt(+15));
	EXPECT_EQ(-3, isqrt(-15));
	EXPECT_EQ(+4, isqrt(+16));
	EXPECT_EQ(-4, isqrt(-16));
	EXPECT_EQ(+1414213562, isqrt(+2000000000000000000LL));
	EXPECT_EQ(-1414213562, isqrt(-2000000000000000000LL));
	EXPECT_EQ(+2000000000, isqrt(+4000000000000000000LL));
	EXPECT_EQ(-2000000000, isqrt(-4000000000000000000LL));
}

TEST(base_test, isqrtc) {
	EXPECT_EQ(0, isqrtc(0));
	EXPECT_EQ(+1, isqrtc(+1));
	EXPECT_EQ(-1, isqrtc(-1));
	EXPECT_EQ(+2, isqrtc(+2));
	EXPECT_EQ(-2, isqrtc(-2));
	EXPECT_EQ(+3, isqrtc(+9));
	EXPECT_EQ(-3, isqrtc(-9));
	EXPECT_EQ(+4, isqrtc(+10));
	EXPECT_EQ(-4, isqrtc(-10));
	EXPECT_EQ(+4, isqrtc(+15));
	EXPECT_EQ(-4, isqrtc(-15));
	EXPECT_EQ(+4, isqrtc(+16));
	EXPECT_EQ(-4, isqrtc(-16));
	EXPECT_EQ(+1414213563, isqrtc(+2000000000000000000LL));
	EXPECT_EQ(-1414213563, isqrtc(-2000000000000000000LL));
	EXPECT_EQ(+2000000000, isqrtc(+4000000000000000000LL));
	EXPECT_EQ(-2000000000, isqrtc(-4000000000000000000LL));
}

TEST(base_test, is_square) {
	EXPECT_TRUE(is_square(0));
	EXPECT_TRUE(is_square(+1));
	EXPECT_FALSE(is_square(-1));
	EXPECT_FALSE(is_square(+2));
	EXPECT_FALSE(is_square(-2));
	EXPECT_TRUE(is_square(+9));
	EXPECT_FALSE(is_square(-9));
	EXPECT_FALSE(is_square(+2000000000000000000LL));
	EXPECT_FALSE(is_square(-2000000000000000000LL));
	EXPECT_TRUE(is_square(+4000000000000000000LL));
	EXPECT_FALSE(is_square(-4000000000000000000LL));
}

template<typename I>
void test_cbT_uint() {
	EXPECT_EQ(I(0), cbT<I>(0));
	EXPECT_EQ(I(1), cbT<I>(+1));
	EXPECT_EQ(I(8), cbT<I>(+2));
	EXPECT_EQ(I(125), cbT<I>(+5));
}

template<typename I>
void test_cbT_int() {
	test_cbT_uint<I>();
	EXPECT_EQ(I(-1), cbT<I>(-1));
	EXPECT_EQ(I(-8), cbT<I>(-2));
	EXPECT_EQ(I(-125), cbT<I>(-5));
}

template<typename F>
void test_cbT_float() {
	test_cbT_int<F>();
	EXPECT_EQ(F(+15.625), cbT<F>(+2.5));
	EXPECT_EQ(F(-15.625), cbT<F>(-2.5));
}

TEST(base_test, cbT) {
	test_cbT_int<int8_t>();
	test_cbT_int<int16_t>();
	test_cbT_int<int32_t>();
	test_cbT_int<int64_t>();
	test_cbT_int<wrapped<int>>();
	test_cbT_uint<uint8_t>();
	test_cbT_uint<uint16_t>();
	test_cbT_uint<uint32_t>();
	test_cbT_uint<uint64_t>();
	test_cbT_float<float>();
	test_cbT_float<double>();
	test_cbT_float<wrapped<float>>();
	EXPECT_EQ(+8000000000000000000LL, cbT(+2000000LL));
	EXPECT_EQ(-8000000000000000000LL, cbT(-2000000LL));
	EXPECT_EQ(+15.625, cbT(+2.5));
	EXPECT_EQ(-15.625, cbT(-2.5));
}

template<typename I>
void test_cbrtT_uint(int max_val = num_iterations) {
	EXPECT_EQ(I(0), cbrtT<I>(0));
	EXPECT_EQ(I(+1), cbrtT<I>(+1));
	EXPECT_EQ(I(+1), cbrtT<I>(+2));
	EXPECT_EQ(I(+1), cbrtT<I>(+3));
	EXPECT_EQ(I(+1), cbrtT<I>(+4));
	EXPECT_EQ(I(+3), cbrtT<I>(+27));
	EXPECT_EQ(I(+3), cbrtT<I>(+28));
	EXPECT_EQ(I(+3), cbrtT<I>(+63));
	EXPECT_EQ(I(+4), cbrtT<I>(+64));
	for (int i = 0; i <= max_val; i++) {
		EXPECT_EQ(I(int(cbrt(double(i)))), cbrtT(I(i))) << "cbrtT(" << i << ")";
	}
}

template<typename I>
void test_cbrtT_int(int max_val = num_iterations) {
	EXPECT_EQ(I(-1), cbrtT<I>(-1));
	EXPECT_EQ(I(-1), cbrtT<I>(-2));
	EXPECT_EQ(I(-1), cbrtT<I>(-3));
	EXPECT_EQ(I(-1), cbrtT<I>(-4));
	EXPECT_EQ(I(-3), cbrtT<I>(-27));
	EXPECT_EQ(I(-3), cbrtT<I>(-28));
	EXPECT_EQ(I(-3), cbrtT<I>(-63));
	EXPECT_EQ(I(-4), cbrtT<I>(-64));
	for (int i = 0; i <= max_val; i++) {
		EXPECT_EQ(I(+int(cbrt(double(i)))), cbrtT(I(+i))) << "cbrtT(" << (+i) << ")";
		EXPECT_EQ(I(-int(cbrt(double(i)))), cbrtT(I(-i))) << "cbrtT(" << (-i) << ")";
	}
}

template<typename F>
void test_cbrtT_float(int max_val = num_iterations) {
	F eps(1e-6f);
	EXPECT_EQ(F(0), cbrtT<F>(0, eps));
	EXPECT_EQ(F(1), cbrtT<F>(1, eps));
	EXPECT_NEAR2(F(1.25992104989487), cbrtT(F(2), eps), eps);
	EXPECT_NEAR2(F(2), cbrtT(F(8), eps), eps);
	EXPECT_NEAR2(F(2.5), cbrtT(F(15.625f), eps), eps);
	for (int i = 0; i <= max_val; i++) {
		F epsi = eps * (abs(i) + 1);
		EXPECT_NEAR2(F(+cbrt(i + 0.0)), cbrtT(F(+(i + 0.0)), eps), epsi) << "cbrtT(" << (+(i + 0.0)) << ")";
		EXPECT_NEAR2(F(-cbrt(i + 0.0)), cbrtT(F(-(i + 0.0)), eps), epsi) << "cbrtT(" << (-(i + 0.0)) << ")";
		EXPECT_NEAR2(F(+cbrt(i + 0.5)), cbrtT(F(+(i + 0.5)), eps), epsi) << "cbrtT(" << (+(i + 0.5)) << ")";
		EXPECT_NEAR2(F(-cbrt(i + 0.5)), cbrtT(F(-(i + 0.5)), eps), epsi) << "cbrtT(" << (-(i + 0.5)) << ")";
	}
}

TEST(base_test, cbrtT) {
	test_cbrtT_int<int8_t>(127);
	test_cbrtT_int<int16_t>();
	test_cbrtT_int<int32_t>();
	test_cbrtT_int<int64_t>();
	test_cbrtT_int<wrapped<int>>();
	test_cbrtT_uint<uint8_t>(127);
	test_cbrtT_uint<uint16_t>();
	test_cbrtT_uint<uint32_t>();
	test_cbrtT_uint<uint64_t>();
	test_cbrtT_float<float>();
	test_cbrtT_float<double>();
	test_cbrtT_float<wrapped<double>>();
	EXPECT_EQ(+1259921LL, cbrtT(+2000000000000000000LL));
	EXPECT_EQ(-1259921LL, cbrtT(-2000000000000000000LL));
	EXPECT_EQ(+2000000LL, cbrtT(+8000000000000000000LL));
	EXPECT_EQ(-2000000LL, cbrtT(-8000000000000000000LL));
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
	EXPECT_EQ(-1, icbrt(-2));
	EXPECT_EQ(+3, icbrt(+27));
	EXPECT_EQ(-3, icbrt(-27));
	EXPECT_EQ(+3, icbrt(+28));
	EXPECT_EQ(-3, icbrt(-28));
	EXPECT_EQ(+3, icbrt(+63));
	EXPECT_EQ(-3, icbrt(-63));
	EXPECT_EQ(+4, icbrt(+64));
	EXPECT_EQ(-4, icbrt(-64));
	EXPECT_EQ(+1259921, icbrt(+2000000000000000000LL));
	EXPECT_EQ(-1259921, icbrt(-2000000000000000000LL));
	EXPECT_EQ(+2000000, icbrt(+8000000000000000000LL));
	EXPECT_EQ(-2000000, icbrt(-8000000000000000000LL));
}

TEST(base_test, icbrtc) {
	EXPECT_EQ(0, icbrtc(0));
	EXPECT_EQ(+1, icbrtc(+1));
	EXPECT_EQ(-1, icbrtc(-1));
	EXPECT_EQ(+2, icbrtc(+2));
	EXPECT_EQ(-2, icbrtc(-2));
	EXPECT_EQ(+3, icbrtc(+27));
	EXPECT_EQ(-3, icbrtc(-27));
	EXPECT_EQ(+4, icbrtc(+28));
	EXPECT_EQ(-4, icbrtc(-28));
	EXPECT_EQ(+4, icbrtc(+63));
	EXPECT_EQ(-4, icbrtc(-63));
	EXPECT_EQ(+4, icbrtc(+64));
	EXPECT_EQ(-4, icbrtc(-64));
	EXPECT_EQ(+1259922, icbrtc(+2000000000000000000LL));
	EXPECT_EQ(-1259922, icbrtc(-2000000000000000000LL));
	EXPECT_EQ(+2000000, icbrtc(+8000000000000000000LL));
	EXPECT_EQ(-2000000, icbrtc(-8000000000000000000LL));
}

TEST(base_test, is_cube) {
	EXPECT_TRUE(is_cube(0));
	EXPECT_TRUE(is_cube(+1));
	EXPECT_TRUE(is_cube(-1));
	EXPECT_FALSE(is_cube(+2));
	EXPECT_FALSE(is_cube(-2));
	EXPECT_TRUE(is_cube(+27));
	EXPECT_TRUE(is_cube(-27));
	EXPECT_FALSE(is_cube(+2000000000000000000LL));
	EXPECT_FALSE(is_cube(-2000000000000000000LL));
	EXPECT_TRUE(is_cube(+8000000000000000000LL));
	EXPECT_TRUE(is_cube(-8000000000000000000LL));
}

TEST(base_test, div_floor) {
	EXPECT_EQ(+4, div_floor(+20, +5));
	EXPECT_EQ(-4, div_floor(+20, -5));
	EXPECT_EQ(-4, div_floor(-20, +5));
	EXPECT_EQ(+4, div_floor(-20, -5));

	EXPECT_EQ(+6, div_floor(+20, +3));
	EXPECT_EQ(-7, div_floor(+20, -3));
	EXPECT_EQ(-7, div_floor(-20, +3));
	EXPECT_EQ(+6, div_floor(-20, -3));
}

TEST(base_test, div_ceil) {
	EXPECT_EQ(+4, div_ceil(+20, +5));
	EXPECT_EQ(-4, div_ceil(+20, -5));
	EXPECT_EQ(-4, div_ceil(-20, +5));
	EXPECT_EQ(+4, div_ceil(-20, -5));

	EXPECT_EQ(+7, div_ceil(+20, +3));
	EXPECT_EQ(-6, div_ceil(+20, -3));
	EXPECT_EQ(-6, div_ceil(-20, +3));
	EXPECT_EQ(+7, div_ceil(-20, -3));
}

TEST(base_test, div_round) {
	EXPECT_EQ(+2, div_round(+100, +50));
	EXPECT_EQ(-2, div_round(+100, -50));
	EXPECT_EQ(-2, div_round(-100, +50));
	EXPECT_EQ(+2, div_round(-100, -50));

	EXPECT_EQ(+3, div_round(+100, +39));
	EXPECT_EQ(-3, div_round(+100, -39));
	EXPECT_EQ(-3, div_round(-100, +39));
	EXPECT_EQ(+3, div_round(-100, -39));

	EXPECT_EQ(+3, div_round(+100, +40));
	EXPECT_EQ(-3, div_round(+100, -40));
	EXPECT_EQ(-3, div_round(-100, +40));
	EXPECT_EQ(+3, div_round(-100, -40));

	EXPECT_EQ(+2, div_round(+100, +41));
	EXPECT_EQ(-2, div_round(+100, -41));
	EXPECT_EQ(-2, div_round(-100, +41));
	EXPECT_EQ(+2, div_round(-100, -41));
}

TEST(base_test, multiple) {
	EXPECT_EQ(+20, multiple(10, +13));
	EXPECT_EQ(-10, multiple(10, -13));
	EXPECT_EQ(+20, multiple(10, +20));
	EXPECT_EQ(-20, multiple(10, -20));
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

TEST(base_test, gcd_max) {
	EXPECT_EQ(0, gcd_max(0, 0));
	EXPECT_EQ(15, gcd_max(15, 0));
	EXPECT_EQ(15, gcd_max(0, 15));
	EXPECT_EQ(25, gcd_max(5 * 5 * 5 * 7, 5 * 5 * 11 * 13));
	EXPECT_EQ(625, gcd_max(5 * 5 * 11, 5 * 5 * 5 * 5 * 7));
	EXPECT_EQ(1, gcd_max(7, 11 * 13));
	EXPECT_EQ(1, gcd_max(11 * 13, 7));
}

TEST(base_test, lcm) {
	EXPECT_EQ(0, lcm(15, 0));
	EXPECT_EQ(0, lcm(0, 15));
	EXPECT_EQ(15, lcm(15, 15));
	EXPECT_EQ(30, lcm(15, 10));
	EXPECT_EQ(210, lcm(14, 15));
}
