#include "altruct/algorithm/collections/collections.h"
#include "altruct/algorithm/math/fft.h"
#include "altruct/structure/math/polynom.h"
#include "altruct/structure/math/modulo.h"
#include "structure_test_util.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;
using namespace altruct::collections;
using namespace altruct::test_util;

namespace {
typedef moduloX<int> modx;
typedef polynom<modx> polyx;

const int M = 1012924417;

vector<moduloX<int>> to_modx(int M, const vector<int>& v) {
	return transform(v, [&](int a) {
		return moduloX<int>(a, M);
	});
}

polyx make_polyx(int M, const vector<int>& v) {
	polyx p(to_modx(M, v));
	p.ZERO_COEFF = modx(0, M);
	return p;
}
}

namespace altruct {
namespace math {
template<>
struct polynom_mul<modx> {
	static void _mul_fft(modx* pr, int lr, const modx* p1, int l1, const modx* p2, int l2) {
		// We can do FFT because of a suitable modulus; 198 ^ (1 << 21) == 1 (mod 1012924417)
		// For a general modulus, we would need to compute several convolutions,
		// each with a suitable modulus, and then combine the results with CRT.
		// Alternatively, one can use complex numbers and break down the input
		// coefficients into 16bit or 11bit words for precission to suffice.
		auto r = convolution<modx>(p1, p1 + l1 + 1, p2, p2 + l2 + 1, modx(198, M), 1 << 21);
		r.resize(lr + 1, modx(0, M)); copy(r.begin(), r.end(), pr);
	}

	static void impl(modx* pr, int lr, const modx* p1, int l1, const modx* p2, int l2) {
		if (l2 < 16) {
			polyx::_mul_long(pr, lr, p1, l1, p2, l2);
		} else if (int64_t(l1) * l2 < 300000) {
			polyx::_mul_karatsuba(pr, lr, p1, l1, p2, l2);
		} else {
			_mul_fft(pr, lr, p1, l1, p2, l2);
		}
	}
};
}
}

TEST(polynom_modx_test, constructor) {
	const vector<modx> c = to_modx(1009, { 1, 2, 3, 4 });
	polyx p0;
	EXPECT_EQ(to_modx(1, { 0 }), p0.c);
	EXPECT_EQ(0, p0.ZERO_COEFF.v);
	EXPECT_EQ(1, p0.ZERO_COEFF.M());
	polyx p1(modx(5, 1009));
	EXPECT_EQ(to_modx(1009, { 5 }), p1.c);
	EXPECT_EQ(0, p1.ZERO_COEFF.v);
	EXPECT_EQ(1009, p1.ZERO_COEFF.M());
	polyx q1(5);
	EXPECT_EQ((vector<modx>{ 5 }), q1.c);
	EXPECT_EQ(0, q1.ZERO_COEFF.v);
	EXPECT_EQ(1, q1.ZERO_COEFF.M());
	polyx p2(c);
	EXPECT_EQ(c, p2.c);
	EXPECT_EQ(0, p2.ZERO_COEFF.v);
	EXPECT_EQ(1009, p2.ZERO_COEFF.M());
	polyx p3(p2);
	EXPECT_EQ(c, p3.c);
	EXPECT_EQ(0, p3.ZERO_COEFF.v);
	EXPECT_EQ(1009, p3.ZERO_COEFF.M());
	polyx p4(c.begin(), c.end());
	EXPECT_EQ(c, p4.c);
	EXPECT_EQ(0, p4.ZERO_COEFF.v);
	EXPECT_EQ(1009, p4.ZERO_COEFF.M());
	polyx p5(&c[0], &c[0] + c.size());
	EXPECT_EQ(c, p5.c);
	EXPECT_EQ(0, p5.ZERO_COEFF.v);
	EXPECT_EQ(1009, p5.ZERO_COEFF.M());
	polyx p6(c.end(), c.end());
	EXPECT_EQ((vector<modx> { }), p6.c);
	EXPECT_EQ(0, p6.ZERO_COEFF.v);
	EXPECT_EQ(1, p6.ZERO_COEFF.M());
	polyx p7{ modx(1, 1009), modx(2, 1009), modx(3, 1009), modx(4, 1009) };
	EXPECT_EQ(c, p7.c);
	EXPECT_EQ(0, p7.ZERO_COEFF.v);
	EXPECT_EQ(1009, p7.ZERO_COEFF.M());
	polyx q7{}; // this calls default constructor!
	EXPECT_EQ(to_modx(1, { 0 }), q7.c);
	EXPECT_EQ(0, q7.ZERO_COEFF.v);
	EXPECT_EQ(1, q7.ZERO_COEFF.M());
	polyx p8(to_modx(1009, { 1, 2, 3, 4 }));
	EXPECT_EQ(c, p8.c);
	EXPECT_EQ(0, p8.ZERO_COEFF.v);
	EXPECT_EQ(1009, p8.ZERO_COEFF.M());
	polyx q8(to_modx(1009, {}));
	EXPECT_EQ(to_modx(1, {}), q8.c);
	EXPECT_EQ(0, q8.ZERO_COEFF.v);
	EXPECT_EQ(1, q8.ZERO_COEFF.M());
	polyx p9(make_polyx(1009, { 1, 2, 3, 4 }));
	EXPECT_EQ(c, p9.c);
	EXPECT_EQ(0, p9.ZERO_COEFF.v);
	EXPECT_EQ(1009, p9.ZERO_COEFF.M());
	polyx q9(polyx(c.end(), c.end()));
	EXPECT_EQ(to_modx(1, {}), q9.c);
	EXPECT_EQ(0, q9.ZERO_COEFF.v);
	EXPECT_EQ(1, q9.ZERO_COEFF.M());
}

TEST(polynom_modx_test, swap) {
	auto p1 = make_polyx(1009, { 1, 2, 3, 4 });
	auto p2 = make_polyx(1003, { 5, 6, 7 });
	p1.swap(p2);
	EXPECT_EQ(to_modx(1003, { 5, 6, 7 }), p1.c);
	EXPECT_EQ(0, p1.ZERO_COEFF.v);
	EXPECT_EQ(1003, p1.ZERO_COEFF.M());
	EXPECT_EQ(to_modx(1009, { 1, 2, 3, 4 }), p2.c);
	EXPECT_EQ(0, p2.ZERO_COEFF.v);
	EXPECT_EQ(1009, p2.ZERO_COEFF.M());
}

TEST(polynom_modx_test, shrink_to_fit) {
	auto p = make_polyx(1009, { 1, 2, 3, 4, 0, 0 });
	EXPECT_EQ(6, p.c.size());
	p.shrink_to_fit();
	EXPECT_EQ(4, p.c.size());
	EXPECT_EQ(to_modx(1009, { 1, 2, 3, 4 }), p.c);
}

TEST(polynom_modx_test, reserve) {
	auto p = make_polyx(1009, { 1, 2, 3, 4 });
	EXPECT_EQ(4, p.c.size());
	p.reserve(6);
	EXPECT_EQ(6, p.c.size());
	EXPECT_EQ(to_modx(1009, { 1, 2, 3, 4, 0, 0 }), p.c);
}

TEST(polynom_modx_test, resize) {
	auto p = make_polyx(1009, { 1, 2, 3, 4, 5 });
	EXPECT_EQ(5, p.c.size());
	p.resize(3);
	EXPECT_EQ(3, p.c.size());
	EXPECT_EQ(to_modx(1009, { 1, 2, 3 }), p.c);
	p.resize(6);
	EXPECT_EQ(6, p.c.size());
	EXPECT_EQ(to_modx(1009, { 1, 2, 3, 0, 0, 0 }), p.c);
}

TEST(polynom_modx_test, size) {
	auto p = make_polyx(1009, { 1, 2, 3, 4 });
	EXPECT_EQ(4, p.size());
}

TEST(polynom_modx_test, at) {
	const auto p = make_polyx(1009, { 2, 3, 5, 7 });
	EXPECT_EQ(modx(2, 1009), p.at(0));
	EXPECT_EQ(modx(7, 1009), p.at(3));
	EXPECT_EQ(modx(0, 1009), p.at(4));
	EXPECT_EQ(modx(0, 1009), p.at(100));
	EXPECT_EQ(4, p.size());
}

TEST(polynom_modx_test, operator_const_brackets) {
	const auto p = make_polyx(1009, { 2, 3, 5, 7 });
	EXPECT_EQ(modx(2, 1009), p[0]);
	EXPECT_EQ(modx(7, 1009), p[3]);
	EXPECT_EQ(modx(0, 1009), p[4]);
	EXPECT_EQ(modx(0, 1009), p[100]);
	EXPECT_EQ(4, p.size());
}

TEST(polynom_modx_test, operator_brackets) {
	auto p = make_polyx(1009, {});
	p[3] = modx(3, 1009);
	EXPECT_EQ(modx(0, 1009), p[0]);
	EXPECT_EQ(modx(0, 1009), p[4]);
	EXPECT_EQ(modx(3, 1009), p[3]);
	EXPECT_EQ(modx(0, 1009), p[4]);
	EXPECT_EQ(modx(0, 1009), p[100]);
	EXPECT_EQ(101, p.size());
}

TEST(polynom_modx_test, degree) {
    const auto p1 = make_polyx(1009, {});
    EXPECT_EQ(0, p1.deg());
    const auto p2 = make_polyx(1009, { 4 });
    EXPECT_EQ(0, p2.deg());
    const auto p3 = make_polyx(1009, { 0, 3 });
    EXPECT_EQ(1, p3.deg());
    const auto p4 = make_polyx(1009, { 2, 3, 5, 7 });
    EXPECT_EQ(3, p4.deg());
    const auto p5 = make_polyx(1009, { 2, 3, 5, 7, 0, 0 });
    EXPECT_EQ(3, p5.deg());
}

TEST(polynom_modx_test, lowest) {
    const auto p1 = make_polyx(1009, {});
    EXPECT_EQ(0, p1.lowest());
    const auto p2 = make_polyx(1009, { 4 });
    EXPECT_EQ(0, p2.lowest());
    const auto p3 = make_polyx(1009, { 0, 3 });
    EXPECT_EQ(1, p3.lowest());
    const auto p4 = make_polyx(1009, { 2, 3, 5, 7 });
    EXPECT_EQ(0, p4.lowest());
    const auto p5 = make_polyx(1009, { 0, 0, 2, 3, 5, 7 });
    EXPECT_EQ(2, p5.lowest());
}

TEST(polynom_modx_test, leading_coefficient) {
	const auto p1 = make_polyx(1009, {});
	EXPECT_EQ(modx(0, 1), p1.leading_coeff());
	const auto p2 = make_polyx(1009, { 4 });
	EXPECT_EQ(modx(4, 1009), p2.leading_coeff());
	const auto p3 = make_polyx(1009, { 0, 3 });
	EXPECT_EQ(modx(3, 1009), p3.leading_coeff());
	const auto p4 = make_polyx(1009, { 2, 3, 5, 7 });
	EXPECT_EQ(modx(7, 1009), p4.leading_coeff());
	const auto p5 = make_polyx(1009, { 2, 3, 5, 7, 0, 0 });
	EXPECT_EQ(modx(7, 1009), p5.leading_coeff());
}

TEST(polynom_modx_test, is_power) {
	auto p1 = make_polyx(1009, {});
	EXPECT_FALSE(p1.is_power());
	const auto p2 = make_polyx(1009, { 4 });
	EXPECT_FALSE(p2.is_power());
	const auto p3 = make_polyx(1009, { 1 });
	EXPECT_TRUE(p3.is_power());
	const auto p4 = make_polyx(1009, { 0, 0, 0, 3 });
	EXPECT_FALSE(p4.is_power());
	const auto p5 = make_polyx(1009, { 0, 0, 0, 1 });
	EXPECT_TRUE(p5.is_power());
	const auto p6 = make_polyx(1009, { 0, 0, 0, 1, 0 });
	EXPECT_TRUE(p6.is_power());
}

TEST(polynom_modx_test, cmp) {
	const auto p0 = make_polyx(1009, {});
	const auto p1 = make_polyx(1009, { 4 });
	const auto p2 = make_polyx(1009, { 1, 3, 5, 7 });
	const auto p3 = make_polyx(1009, { 2, 3, 5, 7 });
	const auto p4 = make_polyx(1009, { 2, 3, 5, 7, 0, 0 });
	const auto p5 = make_polyx(1009, { 2, 3, 6, 7 });
	EXPECT_EQ(0, polyx::cmp(p0, p0));
	EXPECT_EQ(-1, polyx::cmp(p0, p1));
	EXPECT_EQ(+1, polyx::cmp(p1, p0));
	EXPECT_EQ(0, polyx::cmp(p1, p1));
	EXPECT_EQ(-1, polyx::cmp(p0, p2));
	EXPECT_EQ(+1, polyx::cmp(p2, p0));
	EXPECT_EQ(-1, polyx::cmp(p1, p2));
	EXPECT_EQ(+1, polyx::cmp(p2, p1));
	EXPECT_EQ(0, polyx::cmp(p2, p2));
	EXPECT_EQ(-1, polyx::cmp(p2, p3));
	EXPECT_EQ(+1, polyx::cmp(p3, p2));
	EXPECT_EQ(0, polyx::cmp(p3, p4));
	EXPECT_EQ(0, polyx::cmp(p4, p3));
	EXPECT_EQ(-1, polyx::cmp(p4, p5));
	EXPECT_EQ(+1, polyx::cmp(p5, p4));
}

TEST(polynom_modx_test, neg) {
	const auto p0 = make_polyx(1009, {});
	const auto p1 = make_polyx(1009, { 4 });
	const auto p2 = make_polyx(1009, { 1, -3, 5, 7 });
	const auto p3 = make_polyx(1009, { 2, 3, 5, -7, 0, 0 });
	polyx pr;
	polyx::neg(pr, p0);
	EXPECT_EQ(make_polyx(1009, {}), pr);
	polyx::neg(pr, p1);
	EXPECT_EQ(make_polyx(1009, { -4 }), pr);
	polyx::neg(pr, p2);
	EXPECT_EQ(make_polyx(1009, { -1, 3, -5, -7 }), pr);
	polyx::neg(pr, p3);
	EXPECT_EQ(make_polyx(1009, { -2, -3, -5, 7 }), pr);
	// inplace
	polyx::neg(pr, pr);
	EXPECT_EQ(make_polyx(1009, { 2, 3, 5, -7 }), pr);
}

TEST(polynom_modx_test, add) {
	const auto p0 = make_polyx(1009, {});
	const auto p1 = make_polyx(1009, { 4 });
	const auto p2 = make_polyx(1009, { 1, -3, 5, 7 });
	const auto p3 = make_polyx(1009, { 2, 3, 5, -7, 0, 0 });
	polyx pr;
	polyx::add(pr, p0, p0);
	EXPECT_EQ(make_polyx(1009, {}), pr);
	polyx::add(pr, p0, p1);
	EXPECT_EQ(make_polyx(1009, { 4 }), pr);
	polyx::add(pr, p1, p0);
	EXPECT_EQ(make_polyx(1009, { 4 }), pr);
	polyx::add(pr, p1, p1);
	EXPECT_EQ(make_polyx(1009, { 8 }), pr);
	polyx::add(pr, p0, p2);
	EXPECT_EQ(make_polyx(1009, { 1, -3, 5, 7 }), pr);
	polyx::add(pr, p2, p0);
	EXPECT_EQ(make_polyx(1009, { 1, -3, 5, 7 }), pr);
	polyx::add(pr, p1, p2);
	EXPECT_EQ(make_polyx(1009, { 5, -3, 5, 7 }), pr);
	polyx::add(pr, p2, p1);
	EXPECT_EQ(make_polyx(1009, { 5, -3, 5, 7 }), pr);
	polyx::add(pr, p2, p3);
	EXPECT_EQ(make_polyx(1009, { 3, 0, 10 }), pr);
	polyx::add(pr, p3, p2);
	EXPECT_EQ(make_polyx(1009, { 3, 0, 10 }), pr);
	polyx::add(pr, p3, p3);
	EXPECT_EQ(make_polyx(1009, { 4, 6, 10, -14 }), pr);
	// inplace
	pr = to_modx(1009, { 4, 6, 10, -14 });
	polyx::add(pr, pr, p1);
	EXPECT_EQ(make_polyx(1009, { 8, 6, 10, -14 }), pr);
	polyx::add(pr, p1, pr);
	EXPECT_EQ(make_polyx(1009, { 12, 6, 10, -14 }), pr);
	polyx::add(pr, pr, pr);
	EXPECT_EQ(make_polyx(1009, { 24, 12, 20, -28 }), pr);
}

TEST(polynom_modx_test, sub) {
	const auto p0 = make_polyx(1009, {});
	const auto p1 = make_polyx(1009, { 4 });
	const auto p2 = make_polyx(1009, { 1, -3, 5, 7 });
	const auto p3 = make_polyx(1009, { 2, 3, 5, -7, 0, 0 });
	polyx pr;
	polyx::sub(pr, p0, p0);
	EXPECT_EQ(make_polyx(1009, {}), pr);
	polyx::sub(pr, p0, p1);
	EXPECT_EQ(make_polyx(1009, { -4 }), pr);
	polyx::sub(pr, p1, p0);
	EXPECT_EQ(make_polyx(1009, { 4 }), pr);
	polyx::sub(pr, p1, p1);
	EXPECT_EQ(make_polyx(1009, { }), pr);
	polyx::sub(pr, p0, p2);
	EXPECT_EQ(make_polyx(1009, { -1, 3, -5, -7 }), pr);
	polyx::sub(pr, p2, p0);
	EXPECT_EQ(make_polyx(1009, { 1, -3, 5, 7 }), pr);
	polyx::sub(pr, p1, p2);
	EXPECT_EQ(make_polyx(1009, { 3, 3, -5, -7 }), pr);
	polyx::sub(pr, p2, p1);
	EXPECT_EQ(make_polyx(1009, { -3, -3, 5, 7 }), pr);
	polyx::sub(pr, p2, p3);
	EXPECT_EQ(make_polyx(1009, { -1, -6, 0, 14 }), pr);
	polyx::sub(pr, p3, p2);
	EXPECT_EQ(make_polyx(1009, { 1, 6, 0, -14 }), pr);
	polyx::sub(pr, p3, p3);
	EXPECT_EQ(make_polyx(1009, { }), pr);
	// inplace
	pr = to_modx(1009, { 1, 6, 10, -14 });
	polyx::sub(pr, pr, p1);
	EXPECT_EQ(make_polyx(1009, { -3, 6, 10, -14 }), pr);
	polyx::sub(pr, p1, pr);
	EXPECT_EQ(make_polyx(1009, { 7, -6, -10, +14 }), pr);
	polyx::sub(pr, pr, pr);
	EXPECT_EQ(make_polyx(1009, {}), pr);
}

TEST(polynom_modx_test, mul) {
	const auto p0 = make_polyx(1009, {});
	const auto p1 = make_polyx(1009, { 4 });
	const auto p2 = make_polyx(1009, { 1, -3, 5, 7 });
	const auto p3 = make_polyx(1009, { 2, 3, 5, -7, 0, 0 });
	polyx pr;
	polyx::mul(pr, p0, p0);
	EXPECT_EQ(make_polyx(1009, {}), pr);
	polyx::mul(pr, p0, p1);
	EXPECT_EQ(make_polyx(1009, {}), pr);
	polyx::mul(pr, p1, p0);
	EXPECT_EQ(make_polyx(1009, {}), pr);
	polyx::mul(pr, p1, p1);
	EXPECT_EQ(make_polyx(1009, { 16 }), pr);
	polyx::mul(pr, p0, p2);
	EXPECT_EQ(make_polyx(1009, {}), pr);
	polyx::mul(pr, p2, p0);
	EXPECT_EQ(make_polyx(1009, {}), pr);
	polyx::mul(pr, p1, p2);
	EXPECT_EQ(make_polyx(1009, { 4, -12, 20, 28 }), pr);
	polyx::mul(pr, p2, p1);
	EXPECT_EQ(make_polyx(1009, { 4, -12, 20, 28 }), pr);
	polyx::mul(pr, p2, p3);
	EXPECT_EQ(make_polyx(1009, { 2, -3, 6, 7, 67, 0, -49 }), pr);
	polyx::mul(pr, p3, p2);
	EXPECT_EQ(make_polyx(1009, { 2, -3, 6, 7, 67, 0, -49 }), pr);
	polyx::mul(pr, p3, p3);
	EXPECT_EQ(make_polyx(1009, { 4, 12, 29, 2, -17, -70, 49 }), pr);
	// inplace
	pr = make_polyx(1009, { 2, 3, 5, -7 });
	polyx::mul(pr, pr, p1);
	EXPECT_EQ(make_polyx(1009, { 8, 12, 20, -28 }), pr);
	polyx::mul(pr, p1, pr);
	EXPECT_EQ(make_polyx(1009, { 32, 48, 80, -112 }), pr);
	polyx::mul(pr, pr, pr);
	EXPECT_EQ(make_polyx(1009, { 1024, 3072, 7424, 512, -4352, -17920, 12544 }), pr);
}

template<typename P, typename F>
void do_mul(F mul, P& pr, const P& p1, const P& p2, int lr = -1) {
	int l1 = p1.deg(), l2 = p2.deg(); if (lr < 0) lr = l1 + l2;
	pr.resize(lr + 1, p1.ZERO_COEFF);
	mul(pr.c.data(), lr, p1.c.data(), l1, p2.c.data(), l2);
}

TEST(polynom_modx_test, mul_size) {
	auto p1 = make_polyx(M, {}); for (int l = 100; l >= 0; l--) p1[l] = modx(l, M) * (l + 1) / 2;
	auto p2 = make_polyx(M, {}); for (int l = 80; l >= 0; l--) p2[l] = modx(l, M) * 3 + 5;
	auto q11 = p1 * p1;
	auto q12 = p1 * p2;
	polyx q_long; do_mul(polyx::_mul_long, q_long, p1, p2);
	EXPECT_EQ(q12, q_long);
	polyx q_long_inplace = p1; do_mul(polyx::_mul_long, q_long_inplace, q_long_inplace, q_long_inplace);
	EXPECT_EQ(q11, q_long_inplace);
	polyx q_kar; do_mul(polyx::_mul_karatsuba, q_kar, p1, p2);
	EXPECT_EQ(q12, q_kar);
	polyx q_kar_inplace = p1; do_mul(polyx::_mul_karatsuba, q_kar_inplace, q_kar_inplace, q_kar_inplace);
	EXPECT_EQ(q11, q_kar_inplace);
	polyx q_fft; do_mul(polynom_mul<modx>::_mul_fft, q_fft, p1, p2);
	EXPECT_EQ(q12, q_fft);
	polyx q_fft_inplace = p1; do_mul(polynom_mul<modx>::_mul_fft, q_fft_inplace, q_fft_inplace, q_fft_inplace);
	EXPECT_EQ(q11, q_fft_inplace);
	polyx q11_150(q11.c.begin(), q11.c.begin() + 150 + 1);
	polyx q12_150(q12.c.begin(), q12.c.begin() + 150 + 1);
	polyx q_long_150; do_mul(polyx::_mul_long, q_long_150, p1, p2, 150);
	EXPECT_EQ(q12_150, q_long_150);
	polyx q_long_inplace_150 = p1; do_mul(polyx::_mul_long, q_long_inplace_150, q_long_inplace_150, q_long_inplace_150, 150);
	EXPECT_EQ(q11_150, q_long_inplace_150);
	polyx q_kar_150; do_mul(polyx::_mul_karatsuba, q_kar_150, p1, p2, 150);
	EXPECT_EQ(q12_150, q_kar_150);
	polyx q_kar_inplace_150 = p1; do_mul(polyx::_mul_karatsuba, q_kar_inplace_150, q_kar_inplace_150, q_kar_inplace_150, 150);
	EXPECT_EQ(q11_150, q_kar_inplace_150);
	polyx q_fft_150; do_mul(polynom_mul<modx>::_mul_fft, q_fft_150, p1, p2, 150);
	EXPECT_EQ(q12_150, q_fft_150);
	polyx q_fft_inplace_150 = p1; do_mul(polynom_mul<modx>::_mul_fft, q_fft_inplace_150, q_fft_inplace_150, q_fft_inplace_150, 150);
	EXPECT_EQ(q11_150, q_fft_inplace_150);
}

TEST(polynom_modx_test, quot_rem) {
	const auto p0 = make_polyx(1009, {});
	const auto p1 = make_polyx(1009, { 6 });
	const auto p2 = make_polyx(1009, { 1, -3, 0, -2, 0, 0 });
	const auto p3 = make_polyx(1009, { 12, 18, 30, -42, 36, 0, 24, 0 });
	polyx pr;
	polyx::quot_rem(pr, p0, p1);
	EXPECT_EQ(make_polyx(1009, {}), pr);
	polyx::quot_rem(pr, p0, p2);
	EXPECT_EQ(make_polyx(1009, {}), pr);
	polyx::quot_rem(pr, p1, p1);
	EXPECT_EQ(make_polyx(1009, { 1 }), pr);
	polyx::quot_rem(pr, p1, p2);
	EXPECT_EQ(make_polyx(1009, { 6 }), pr);
	polyx::quot_rem(pr, p3, p1);
	EXPECT_EQ(make_polyx(1009, { 2, 3, 5, -7, 6, 0, 4 }), pr);
	polyx::quot_rem(pr, p3, p2);
	EXPECT_EQ(make_polyx(1009, {-3, 63, 30, 15, 0, 0, -12}), pr);
	polyx::quot_rem(pr, p2, p3);
	EXPECT_EQ(make_polyx(1009, { 1, -3, 0, -2 }), pr);
	// inplace
	pr = p3;
	polyx::quot_rem(pr, pr, p2);
	EXPECT_EQ(make_polyx(1009, {-3, 63, 30, 15, 0, 0, -12}), pr);
	pr = p2;
	polyx::quot_rem(pr, pr, p3);
	EXPECT_EQ(make_polyx(1009, { 1, -3, 0, -2 }), pr);
	pr = p3;
	polyx::quot_rem(pr, pr, p3);
	EXPECT_EQ(make_polyx(1009, { 0, 0, 0, 0, 0, 0, 1 }), pr);
}

TEST(polynom_modx_test, div) {
	const auto p0 = make_polyx(1009, {});
	const auto p1 = make_polyx(1009, { 6 });
	const auto p2 = make_polyx(1009, { 1, -3, 0, -2, 0, 0 });
	const auto p3 = make_polyx(1009, { 12, 18, 30, -42, 36, 0, 24, 0 });
	polyx pr;
	polyx::div(pr, p0, p1);
	EXPECT_EQ(make_polyx(1009, {}), pr);
	polyx::div(pr, p0, p2);
	EXPECT_EQ(make_polyx(1009, {}), pr);
	polyx::div(pr, p1, p1);
	EXPECT_EQ(make_polyx(1009, { 1 }), pr);
	polyx::div(pr, p1, p2);
	EXPECT_EQ(make_polyx(1009, { }), pr);
	polyx::div(pr, p3, p1);
	EXPECT_EQ(make_polyx(1009, { 2, 3, 5, -7, 6, 0, 4 }), pr);
	polyx::div(pr, p3, p2);
	EXPECT_EQ(make_polyx(1009, { 15, 0, 0, -12}), pr);
	polyx::div(pr, p2, p3);
	EXPECT_EQ(make_polyx(1009, {}), pr);
	// inplace
	pr = p3;
	polyx::div(pr, pr, p2);
	EXPECT_EQ(make_polyx(1009, { 15, 0, 0, -12}), pr);
	pr = p2;
	polyx::div(pr, pr, p3);
	EXPECT_EQ(make_polyx(1009, {}), pr);
	pr = p3;
	polyx::div(pr, pr, p3);
	EXPECT_EQ(make_polyx(1009, { 1 }), pr);
}

TEST(polynom_modx_test, mod) {
	const auto p0 = make_polyx(1009, {});
	const auto p1 = make_polyx(1009, { 6 });
	const auto p2 = make_polyx(1009, { 1, -3, 0, -2, 0, 0 });
	const auto p3 = make_polyx(1009, { 12, 18, 30, -42, 36, 0, 24, 0 });
	polyx pr;
	polyx::mod(pr, p0, p1);
	EXPECT_EQ(make_polyx(1009, {}), pr);
	polyx::mod(pr, p0, p2);
	EXPECT_EQ(make_polyx(1009, {}), pr);
	polyx::mod(pr, p1, p1);
	EXPECT_EQ(make_polyx(1009, {}), pr);
	polyx::mod(pr, p1, p2);
	EXPECT_EQ(make_polyx(1009, { 6 }), pr);
	polyx::mod(pr, p3, p1);
	EXPECT_EQ(make_polyx(1009, {}), pr);
	polyx::mod(pr, p3, p2);
	EXPECT_EQ(make_polyx(1009, {-3, 63, 30 }), pr);
	polyx::mod(pr, p2, p3);
	EXPECT_EQ(make_polyx(1009, { 1, -3, 0, -2 }), pr);
	// inplace
	pr = p3;
	polyx::mod(pr, pr, p2);
	EXPECT_EQ(make_polyx(1009, {-3, 63, 30 }), pr);
	pr = p2;
	polyx::mod(pr, pr, p3);
	EXPECT_EQ(make_polyx(1009, { 1, -3, 0, -2 }), pr);
	pr = p3;
	polyx::mod(pr, pr, p3);
	EXPECT_EQ(make_polyx(1009, { 0, 0, 0, 0, 0, 0 }), pr);
}

TEST(polynom_modx_test, muls) {
	const auto p0 = make_polyx(1009, {});
	const auto p1 = make_polyx(1009, { 4 });
	const auto p2 = make_polyx(1009, { 1, -3, 5, 7 });
	const auto p3 = make_polyx(1009, { 2, 3, 5, -7, 0, 0 });
	polyx pr;
	polyx::mul(pr, p0, modx(11, 1009));
	EXPECT_EQ(make_polyx(1009, {}), pr);
	polyx::mul(pr, p1, modx(11, 1009));
	EXPECT_EQ(make_polyx(1009, { 44 }), pr);
	polyx::mul(pr, p2, modx(11, 1009));
	EXPECT_EQ(make_polyx(1009, { 11, -33, 55, 77 }), pr);
	polyx::mul(pr, p3, modx(11, 1009));
	EXPECT_EQ(make_polyx(1009, { 22, 33, 55, -77 }), pr);
	// inplace
	pr = make_polyx(1009, { 2, 3, 5, -7 });
	polyx::mul(pr, pr, modx(11, 1009));
	EXPECT_EQ(make_polyx(1009, { 22, 33, 55, -77 }), pr);
}

TEST(polynom_modx_test, divs) {
	const auto p0 = make_polyx(1009, {});
	const auto p1 = make_polyx(1009, { 44 });
	const auto p2 = make_polyx(1009, { 11, -33, 55, 77 });
	const auto p3 = make_polyx(1009, { 22, 33, 55, -77, 0, 0 });
	polyx pr;
	polyx::div(pr, p0, modx(11, 1009));
	EXPECT_EQ(make_polyx(1009, {}), pr);
	polyx::div(pr, p1, modx(11, 1009));
	EXPECT_EQ(make_polyx(1009, { 4 }), pr);
	polyx::div(pr, p2, modx(11, 1009));
	EXPECT_EQ(make_polyx(1009, { 1, -3, 5, 7 }), pr);
	polyx::div(pr, p3, modx(11, 1009));
	EXPECT_EQ(make_polyx(1009, { 2, 3, 5, -7 }), pr);
	// inplace
	pr = make_polyx(1009, { 22, 33, 55, -77 });
	polyx::div(pr, pr, modx(11, 1009));
	EXPECT_EQ(make_polyx(1009, { 2, 3, 5, -7 }), pr);
}

TEST(polynom_modx_test, operators_comparison) {
	const auto p1 = make_polyx(1009, { 4 });
	const auto p2 = make_polyx(1009, { 1, 3, 5, 7 });
	const auto p3 = make_polyx(1009, { 1, 3, 5, 7, 0, 0, 0 });
	ASSERT_COMPARISON_OPERATORS(0, p1, p1);
	ASSERT_COMPARISON_OPERATORS(0, p2, p2);
	ASSERT_COMPARISON_OPERATORS(0, p3, p3);
	ASSERT_COMPARISON_OPERATORS(-1, p1, p2);
	ASSERT_COMPARISON_OPERATORS(+1, p2, p1);
	ASSERT_COMPARISON_OPERATORS(0, p2, p3);
	ASSERT_COMPARISON_OPERATORS(0, p3, p2);
}

TEST(polynom_modx_test, operators_arithmetic) {
	const auto p1 = make_polyx(1009, { 4, 1 });
	const auto p2 = make_polyx(1009, { 1, -3, 5, 7 });
	const auto p3 = make_polyx(1009, { 11, -33, 55, 77 });
	EXPECT_EQ(make_polyx(1009, { 5, -2, 5, 7 }), p1 + p2);
	EXPECT_EQ(make_polyx(1009, { 3, 4, -5, -7 }), p1 - p2);
	EXPECT_EQ(make_polyx(1009, { -1, 3, -5, -7 }), -p2);
	EXPECT_EQ(make_polyx(1009, { 4, -11, 17, 33, 7 }), p1 * p2);
	EXPECT_EQ(make_polyx(1009, { }), p1 / p2);
	EXPECT_EQ(make_polyx(1009, { 4, 1 }), p1 % p2);
	EXPECT_EQ(make_polyx(1009, { 5, -2, 5, 7 }), p2 + p1);
	EXPECT_EQ(make_polyx(1009, { -3, -4, 5, 7 }), p2 - p1);
	EXPECT_EQ(make_polyx(1009, { -4, -1 }), -p1);
	EXPECT_EQ(make_polyx(1009, { 4, -11, 17, 33, 7 }), p2 * p1);
	EXPECT_EQ(make_polyx(1009, { 89, -23, 7 }), p2 / p1);
	EXPECT_EQ(make_polyx(1009, { -355 }), p2 % p1);
	EXPECT_EQ(make_polyx(1009, { 11, -33, 55, 77 }), p2 * modx(11, 1009));
	EXPECT_EQ(make_polyx(1009, { 1, -3, 5, 7 }), p3 / modx(11, 1009));
}

TEST(polynom_modx_test, operators_inplace) {
	const auto p1 = make_polyx(1009, { 4, 1 });
	const auto p2 = make_polyx(1009, { 1, -3, 5, 7 });
	const auto p3 = make_polyx(1009, { 11, -33, 55, 77 });
	polyx pr;
	pr = p1; pr += p2;
	EXPECT_EQ(make_polyx(1009, { 5, -2, 5, 7 }), pr);
	pr = p1; pr -= p2;
	EXPECT_EQ(make_polyx(1009, { 3, 4, -5, -7 }), pr);
	pr = p1; pr *= p2;
	EXPECT_EQ(make_polyx(1009, { 4, -11, 17, 33, 7 }), pr);
	pr = p1; pr /= p2;
	EXPECT_EQ(make_polyx(1009, { }), pr);
	pr = p1; pr %= p2;
	EXPECT_EQ(make_polyx(1009, { 4, 1 }), pr);
	pr = p2; pr += p1;
	EXPECT_EQ(make_polyx(1009, { 5, -2, 5, 7 }), pr);
	pr = p2; pr -= p1;
	EXPECT_EQ(make_polyx(1009, { -3, -4, 5, 7 }), pr);
	pr = p2; pr *= p1;
	EXPECT_EQ(make_polyx(1009, { 4, -11, 17, 33, 7 }), pr);
	pr = p2; pr /= p1;
	EXPECT_EQ(make_polyx(1009, { 89, -23, 7 }), pr);
	pr = p2; pr %= p1;
	EXPECT_EQ(make_polyx(1009, { -355 }), pr);
	pr = p2; pr *= modx(11, 1009);
	EXPECT_EQ(make_polyx(1009, { 11, -33, 55, 77 }), pr);
	pr = p3; pr /= modx(11, 1009);
	EXPECT_EQ(make_polyx(1009, { 1, -3, 5, 7 }), pr);
}

TEST(polynom_modx_test, operators_inplace_self) {
	const auto p1 = make_polyx(1009, { 2, -3, 5, 7 });
	polyx pr;
	pr = p1; pr += pr;
	EXPECT_EQ(make_polyx(1009, { 4, -6, 10, 14 }), pr);
	pr = p1; pr -= pr;
	EXPECT_EQ(make_polyx(1009, { 0, 0, 0, 0 }), pr);
	pr = p1; pr *= pr;
	EXPECT_EQ(make_polyx(1009, { 4, -12, 29, -2, -17, 70, 49 }), pr);
	pr = p1; pr %= pr;
	EXPECT_EQ(make_polyx(1009, { 0, 0, 0 }), pr);
	pr = p1; pr /= pr;
	EXPECT_EQ(make_polyx(1009, { 1 }), pr);
}

TEST(polynom_modx_test, eval) {
	const auto p1 = make_polyx(1009, { 7, 5, -3, 2 });
	vector<modx> ve1;
	for (int x = -3; x <= 4; x++) {
		ve1.push_back(p1.eval(modx(x, 1009)));
	}
	EXPECT_EQ(to_modx(1009, {-89, -31, -3, 7, 11, 21, 49, 107}), ve1);
}

TEST(polynom_modx_test, derivative) {
	const auto p1 = make_polyx(1009, { 7, 5, -3, 4 });
	const auto pd = p1.derivative();
	EXPECT_EQ(make_polyx(1009, { 5, -6, 12 }), pd);
}

TEST(polynom_modx_test, integral) {
	const auto p = make_polyx(1009, { 7, 8, 15, -4, 20 });
	const auto pi0 = p.integral();
	const auto pi3 = p.integral(3);
	EXPECT_EQ(make_polyx(1009, { 0, 7, 4, 5, -1, 4 }), pi0);
	EXPECT_EQ(make_polyx(1009, { 3, 7, 4, 5, -1, 4 }), pi3);
}

TEST(polynom_modx_test, casts) {
	polyx p1{ { 2, 1009 }, { 3, 1009 }, { 5, 1009 } };
	EXPECT_EQ(2, p1[0].v);
	EXPECT_EQ(1009, p1[0].M());
	EXPECT_EQ(3, p1[1].v);
	EXPECT_EQ(1009, p1[1].M());
	EXPECT_EQ(5, p1[2].v);
	EXPECT_EQ(1009, p1[2].M());
	EXPECT_EQ(0, p1[3].v);
	EXPECT_EQ(1009, p1[3].M());
    EXPECT_EQ(2, p1.deg());
	polyx e0 = zeroT<polyx>::of(p1);
	EXPECT_EQ(0, e0[0].v);
	EXPECT_EQ(1009, e0[0].M());
	EXPECT_EQ(0, e0.deg());
	polyx e1 = identityT<polyx>::of(p1);
	EXPECT_EQ(1, e1[0].v);
	EXPECT_EQ(1009, e1[0].M());
	EXPECT_EQ(0, e1.deg());
    const polyx p2 = castOf<polyx>(p1);
    EXPECT_EQ(2, p2[0].v);
    EXPECT_EQ(1009, p2[0].M());
    EXPECT_EQ(3, p2[1].v);
    EXPECT_EQ(1009, p2[1].M());
    EXPECT_EQ(5, p2[2].v);
    EXPECT_EQ(1009, p2[2].M());
    EXPECT_EQ(0, p2[3].v);
    EXPECT_EQ(1009, p2[3].M());
    EXPECT_EQ(2, p2.deg());
    const polyx p3 = castOf(e1, p1);
    EXPECT_EQ(2, p3[0].v);
    EXPECT_EQ(1009, p3[0].M());
    EXPECT_EQ(3, p3[1].v);
    EXPECT_EQ(1009, p3[1].M());
    EXPECT_EQ(5, p3[2].v);
    EXPECT_EQ(1009, p3[2].M());
    EXPECT_EQ(0, p3[3].v);
    EXPECT_EQ(1009, p3[3].M());
    EXPECT_EQ(2, p3.deg());
    const polyx p4 = castOf(e1, 4);
    EXPECT_EQ(4, p4[0].v);
    EXPECT_EQ(1009, p4[0].M());
    EXPECT_EQ(0, p4[1].v);
    EXPECT_EQ(1009, p4[1].M());
    EXPECT_EQ(0, p4.deg());
}
