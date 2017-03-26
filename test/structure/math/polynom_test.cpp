#include "algorithm/math/fft.h"
#include "structure/math/polynom.h"
#include "structure/math/modulo.h"
#include "structure_test_util.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;
using namespace altruct::test_util;

namespace {
class A {
public:
	double v;
	A(double v) : v(v) {}
	A(int v) : v(v) {}
	bool operator == (const A& rhs) const { return v == rhs.v; }
};

typedef modulo<int, 1012924417> mod;
}

namespace altruct {
namespace math {
template<>
struct polynom_mul<mod> {
	static void _mul_fft(mod* pr, int lr, const mod* p1, int l1, const mod* p2, int l2) {
		// We can do FFT because of a suitable modulus; 198 ^ (1 << 21) == 1 (mod 1012924417)
		// For a general modulus, we would need to compute several convolutions,
		// each with a suitable modulus, and then combine the results with CRT.
		// Alternatively, one can use complex numbers and break down the input
		// coefficients into 16bit or 11bit words for precission to suffice.
		auto r = convolution<mod>(p1, p1 + l1 + 1, p2, p2 + l2 + 1, mod(198), 1 << 21);
		r.resize(lr + 1, mod(0)); copy(r.begin(), r.end(), pr);
	}

	static void impl(mod* pr, int lr, const mod* p1, int l1, const mod* p2, int l2) {
		if (l2 < 16) {
			polynom<mod>::_mul_long(pr, lr, p1, l1, p2, l2);
		} else if (int64_t(l1) * l2 < 300000) {
			polynom<mod>::_mul_karatsuba(pr, lr, p1, l1, p2, l2);
		} else {
			_mul_fft(pr, lr, p1, l1, p2, l2);
		}
	}
};
}
}

TEST(polynom_test, constructor) {
	const vector<int> c = { 1, 2, 3, 4 };
	polynom<int> p0;
	EXPECT_EQ((vector<int>{ 0 }), p0.c);
	EXPECT_EQ(0, p0.ZERO_COEFF);
	polynom<int> p1(5);
	EXPECT_EQ((vector<int>{ 5 }), p1.c);
	EXPECT_EQ(0, p1.ZERO_COEFF);
	polynom<int> p2(c);
	EXPECT_EQ(c, p2.c);
	EXPECT_EQ(0, p2.ZERO_COEFF);
	polynom<int> p3(p2);
	EXPECT_EQ(c, p3.c);
	EXPECT_EQ(0, p3.ZERO_COEFF);
	polynom<int> p4(c.begin(), c.end());
	EXPECT_EQ(c, p4.c);
	EXPECT_EQ(0, p4.ZERO_COEFF);
	polynom<int> p5(&c[0], &c[0] + c.size());
	EXPECT_EQ(c, p5.c);
	EXPECT_EQ(0, p5.ZERO_COEFF);
	polynom<int> p6(c.end(), c.end());
	EXPECT_EQ((vector<int> { }), p6.c);
	EXPECT_EQ(0, p6.ZERO_COEFF);
	polynom<int> p7{ 1, 2, 3, 4 };
	EXPECT_EQ(c, p7.c);
	EXPECT_EQ(0, p7.ZERO_COEFF);
	polynom<A> q1(5);
	EXPECT_EQ((vector<A>{5.0}), q1.c);
	EXPECT_EQ(A(0), q1.ZERO_COEFF);
	polynom<A> q2(A(5.3));
	EXPECT_EQ((vector<A>{5.3}), q2.c);
	EXPECT_EQ(A(0), q2.ZERO_COEFF);
	polynom<int> p8(vector<int>{ 1, 2, 3, 4 });
	EXPECT_EQ(c, p8.c);
	EXPECT_EQ(0, p8.ZERO_COEFF);
	polynom<int> p9(polynom<int>{ 1, 2, 3, 4 });
	EXPECT_EQ(c, p9.c);
	EXPECT_EQ(0, p9.ZERO_COEFF);
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

TEST(polynom_test, resize) {
	polynom<int> p{ 1, 2, 3, 4, 5 };
	EXPECT_EQ(5, p.c.size());
	p.resize(3);
	EXPECT_EQ(3, p.c.size());
	EXPECT_EQ((vector<int>{ 1, 2, 3 }), p.c);
	p.resize(6);
	EXPECT_EQ(6, p.c.size());
	EXPECT_EQ((vector<int>{ 1, 2, 3, 0, 0, 0 }), p.c);
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
	EXPECT_FALSE(p1.is_power());
	const polynom<int> p2{ 4 };
	EXPECT_FALSE(p2.is_power());
	const polynom<int> p3{ 1 };
	EXPECT_TRUE(p3.is_power());
	const polynom<int> p4{ 0, 0, 0, 3 };
	EXPECT_FALSE(p4.is_power());
	const polynom<int> p5{ 0, 0, 0, 1 };
	EXPECT_TRUE(p5.is_power());
	const polynom<int> p6{ 0, 0, 0, 1, 0 };
	EXPECT_TRUE(p6.is_power());
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

template<typename P, typename F>
void do_mul(F mul, P& pr, const P& p1, const P& p2, int lr = -1) {
	int l1 = p1.deg(), l2 = p2.deg(); if (lr < 0) lr = l1 + l2;
	pr.c.resize(lr + 1, pr.ZERO_COEFF);
	mul(pr.c.data(), lr, p1.c.data(), l1, p2.c.data(), l2);
}

TEST(polynom_test, mul_size) {
	polynom<mod> p1; for (int l = 100; l >= 0; l--) p1[l] = mod(l) * (l - 1) / 2;
	polynom<mod> p2; for (int l = 80; l >= 0; l--) p2[l] = mod(l) * 3 + 5;
	polynom<mod> q11 = p1 * p1;
	polynom<mod> q12 = p1 * p2;
	polynom<mod> q_long; do_mul(polynom<mod>::_mul_long, q_long, p1, p2);
	EXPECT_EQ(q12, q_long);
	polynom<mod> q_long_inplace = p1; do_mul(polynom<mod>::_mul_long, q_long_inplace, q_long_inplace, q_long_inplace);
	EXPECT_EQ(q11, q_long_inplace);
	polynom<mod> q_kar; do_mul(polynom<mod>::_mul_karatsuba, q_kar, p1, p2);
	EXPECT_EQ(q12, q_kar);
	polynom<mod> q_kar_inplace = p1; do_mul(polynom<mod>::_mul_karatsuba, q_kar_inplace, q_kar_inplace, q_kar_inplace);
	EXPECT_EQ(q11, q_kar_inplace);
	polynom<mod> q_fft; do_mul(polynom_mul<mod>::_mul_fft, q_fft, p1, p2);
	EXPECT_EQ(q12, q_fft);
	polynom<mod> q_fft_inplace = p1; do_mul(polynom_mul<mod>::_mul_fft, q_fft_inplace, q_fft_inplace, q_fft_inplace);
	EXPECT_EQ(q11, q_fft_inplace);
	polynom<mod> q11_150(q11.c.begin(), q11.c.begin() + 150 + 1);
	polynom<mod> q12_150(q12.c.begin(), q12.c.begin() + 150 + 1);
	polynom<mod> q_long_150; do_mul(polynom<mod>::_mul_long, q_long_150, p1, p2, 150);
	EXPECT_EQ(q12_150, q_long_150);
	polynom<mod> q_long_inplace_150 = p1; do_mul(polynom<mod>::_mul_long, q_long_inplace_150, q_long_inplace_150, q_long_inplace_150, 150);
	EXPECT_EQ(q11_150, q_long_inplace_150);
	polynom<mod> q_kar_150; do_mul(polynom<mod>::_mul_karatsuba, q_kar_150, p1, p2, 150);
	EXPECT_EQ(q12_150, q_kar_150);
	polynom<mod> q_kar_inplace_150 = p1; do_mul(polynom<mod>::_mul_karatsuba, q_kar_inplace_150, q_kar_inplace_150, q_kar_inplace_150, 150);
	EXPECT_EQ(q11_150, q_kar_inplace_150);
	polynom<mod> q_fft_150; do_mul(polynom_mul<mod>::_mul_fft, q_fft_150, p1, p2, 150);
	EXPECT_EQ(q12_150, q_fft_150);
	polynom<mod> q_fft_inplace_150 = p1; do_mul(polynom_mul<mod>::_mul_fft, q_fft_inplace_150, q_fft_inplace_150, q_fft_inplace_150, 150);
	EXPECT_EQ(q11_150, q_fft_inplace_150);
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
	const polynom<int> p3{ 1, 3, 5, 7, 0, 0, 0 };
	ASSERT_COMPARISON_OPERATORS(0, p1, p1);
	ASSERT_COMPARISON_OPERATORS(0, p2, p2);
	ASSERT_COMPARISON_OPERATORS(0, p3, p3);
	ASSERT_COMPARISON_OPERATORS(-1, p1, p2);
	ASSERT_COMPARISON_OPERATORS(+1, p2, p1);
	ASSERT_COMPARISON_OPERATORS(0, p2, p3);
	ASSERT_COMPARISON_OPERATORS(0, p3, p2);
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
	EXPECT_EQ(3.5, p1(-0.5));
}

TEST(polynom_test, derivative) {
	const polynom<int> p1{ 7, 5, -3, 4 };
	const polynom<int> pd = p1.derivative();
	EXPECT_EQ((polynom<int>{ 5, -6, 12 }), pd);
}

TEST(polynom_test, integral) {
	const polynom<int> p{ 7, 8, 15, -4, 20 };
	const polynom<int> pi0 = p.integral();
	const polynom<int> pi3 = p.integral(3);
	EXPECT_EQ((polynom<int>{ 0, 7, 4, 5, -1, 4 }), pi0);
	EXPECT_EQ((polynom<int>{ 3, 7, 4, 5, -1, 4 }), pi3);
}

TEST(polynom_test, casts) {
	typedef polynom<int> poly;
	poly p1{ 2, 3, 5 };
	EXPECT_EQ(2, p1[0]);
	EXPECT_EQ(3, p1[1]);
	EXPECT_EQ(5, p1[2]);
	EXPECT_EQ(0, p1[3]);
    EXPECT_EQ(2, p1.deg());
	const poly e0 = zeroT<poly>::of(p1);
	EXPECT_EQ(0, e0[0]);
	EXPECT_EQ(0, e0.deg());
	const poly e1 = identityT<poly>::of(p1);
	EXPECT_EQ(1, e1[0]);
	EXPECT_EQ(0, e1.deg());
    const poly p2 = castOf<poly>(p1);
    EXPECT_EQ(2, p2[0]);
    EXPECT_EQ(3, p2[1]);
    EXPECT_EQ(5, p2[2]);
    EXPECT_EQ(0, p2[3]);
    EXPECT_EQ(2, p2.deg());
    const poly p3 = castOf(e1, p1);
    EXPECT_EQ(2, p3[0]);
    EXPECT_EQ(3, p3[1]);
    EXPECT_EQ(5, p3[2]);
    EXPECT_EQ(0, p3[3]);
    EXPECT_EQ(2, p3.deg());
    const poly p4 = castOf(e1, 4);
    EXPECT_EQ(4, p4[0]);
    EXPECT_EQ(0, p4[1]);
    EXPECT_EQ(0, p4.deg());
    const poly p5 = castOf<poly>(5);
    EXPECT_EQ(5, p5[0]);
    EXPECT_EQ(0, p5[1]);
    EXPECT_EQ(0, p5.deg());
}
