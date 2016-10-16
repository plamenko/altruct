#include "structure/math/double_int.h"
#include "algorithm/random/xorshift.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;
using namespace altruct::random;

// MSVS
struct intr8 {
	uint8_t adc(uint8_t x, uint8_t y, int& carry) const {
		uint8_t r; carry = _addcarry_u8(carry, x, y, &r); return r;
	}
	uint8_t sbb(uint8_t x, uint8_t y, int& borrow) const {
		uint8_t r; borrow = _subborrow_u8(borrow, x, y, &r); return r;
	}
	uint8_t umul(uint8_t x, uint8_t y, uint8_t* hi) const {
		uint16_t r = uint16_t(x) * y;
		*hi = r >> 8; return uint8_t(r);
	}
};

// MSVS
struct intr64 {
	uint64_t adc(uint64_t x, uint64_t y, int& carry) const {
		uint64_t r; carry = _addcarry_u64(carry, x, y, &r); return r;
	}
	uint64_t sbb(uint64_t x, uint64_t y, int& borrow) const {
		uint64_t r; borrow = _subborrow_u64(borrow, x, y, &r); return r;
	}
	uint64_t umul(uint64_t x, uint64_t y, uint64_t* hi) const {
		return _umul128(x, y, hi);
	}
};

typedef prim_int<int8_t, uint8_t, intr8> l8_8;
typedef double_int<l8_8> l8_16;
typedef double_int<l8_16> l8_32;
typedef double_int<l8_32> l8_64;
typedef double_int<l8_64> l8_128;
typedef double_int<l8_128> l8_256;

typedef prim_int<int64_t, uint64_t, intr64> l64_64;
typedef double_int<l64_64> l64_128;
typedef double_int<l64_128> l64_256;
typedef double_int<l64_256> l64_512;
typedef double_int<l64_512> l64_1024;

TEST(double_int_test, correctness_8_32) {
	return; // don't test by default
	for (int x = 0; x < 65536; x++) {
		for (int y = 0; y < 65536; y++) {
			uint32_t z0 = uint32_t(x) * y;
			uint64_t z1 = l8_16::unsigned_mul_full(l8_16(x >> 8, x), l8_16(y >> 8, y)).to_uint64();
			if (z1 != uint64_t(z0)) {
				fprintf(stderr, "ERROR: %d * %d\n", x, y);
			}
		}
	}
	fprintf(stderr, "%d ms\n", clock());
	// verification of all 2^32  16 x 16-> 32 bit  multiplications  56.551s  ~ 76 Mops
}

template<typename T, typename F_OP, typename F_INIT>
void test_perf(const string& s, F_OP f_op, F_INIT f_init, int iter = 100000) {
	uint64_t cnt = 0;
	auto t = clock();
	T x, y, z;
	while (clock() - t < CLOCKS_PER_SEC * 1) {
		x = f_init(), y = f_init();
		for (int i = 0; i < iter; i++) {
			z = f_op(x, y); x = y, y = z;
		}
		cnt += iter;
	}
	double ops = cnt / (double(clock() - t) / CLOCKS_PER_SEC);
	fprintf(stderr, "%s\t %.3lf Mops  %p\n", s.c_str(), ops / 1e6, &z);
}

template<typename F_OP>
void test_perf(const string& s, F_OP f_op, int iter = 100000) {
	xorshift_64star rnd(12345);
	test_perf<uint32_t>(s + " uint32_t", f_op, [&]() {
		return uint32_t(rnd.next());
	});
	test_perf<uint64_t>(s + " uint64_t", f_op, [&]() {
		return uint64_t(rnd.next());
	});
	test_perf<l64_64>(s + " l64_64", f_op, [&]() {
		return l64_64{ rnd.next() };
	});
	test_perf<l64_128>(s + " l64_128", f_op, [&]() {
		return l64_128{ rnd.next(), rnd.next() };
	});
	test_perf<l64_256>(s + " l64_256", f_op, [&]() {
		return l64_256{ { rnd.next(), rnd.next() }, { rnd.next(), rnd.next() } };
	});
	test_perf<l64_512>(s + " l64_512", f_op, [&]() {
		return l64_512{ { { rnd.next(), rnd.next() }, { rnd.next(), rnd.next() } }, { { rnd.next(), rnd.next() }, { rnd.next(), rnd.next() } } };
	});
	test_perf<l64_1024>(s + " l64_1024", f_op, [&]() {
		return l64_1024{
			{ { { rnd.next(), rnd.next() }, { rnd.next(), rnd.next() } }, { { rnd.next(), rnd.next() }, { rnd.next(), rnd.next() } } },
		    { { { rnd.next(), rnd.next() }, { rnd.next(), rnd.next() } }, { { rnd.next(), rnd.next() }, { rnd.next(), rnd.next() } } } };
	});
}

struct f_mul { template<typename T> T operator() (const T& x, const T& y) const { return x * y; } };
struct f_add { template<typename T> T operator() (const T& x, const T& y) const { return x + y; } };
TEST(double_int_test, perf) {
	//xorshift_64star rnd(12345);
	//const int iter = 100000;
	//uint64_t cnt = 0;
	//auto t = clock();
	//l64_128 x, y, z;
	//while (clock() - t < CLOCKS_PER_SEC * 1) {
	//	x = { rnd.next(), rnd.next() }, y = { rnd.next(), rnd.next() };
	//	for (int i = 0; i < iter; i++) {
	//		z = x + y; x = y, y = z;
	//	}
	//	cnt += iter;
	//}
	//double ops = cnt / (double(clock() - t) / CLOCKS_PER_SEC);
	//fprintf(stderr, "%s\t %.3lf Mops  %p\n", "l64_128 add ", ops / 1e6, &z);
	
	return; // don't test perf by default
	test_perf("add", f_add());
	test_perf("mul", f_mul());
}

TEST(double_int_test, correctness_64_512) {
	xorshift_64star rnd(12345);
	//l64_512 x = { { { rnd.next(), rnd.next() }, { rnd.next(), rnd.next() } }, { { rnd.next(), rnd.next() }, { rnd.next(), rnd.next() } } };
	//l64_512 y = { { { rnd.next(), rnd.next() }, { rnd.next(), rnd.next() } }, { { rnd.next(), rnd.next() }, { rnd.next(), rnd.next() } } };
	l64_256 x = { { rnd.next(), rnd.next() }, { rnd.next(), rnd.next() } };
	l64_256 y = { { rnd.next(), rnd.next() }, { rnd.next(), rnd.next() } };
	auto z = x + y;
	cout << x.to_string16() << endl;
	cout << y.to_string16() << endl;
	cout << z.to_string16() << endl;
}


//TEST(double_int_test, constructor) {
//	quad x1;
//	EXPECT_EQ(0, x1.a);
//	EXPECT_EQ(0, x1.b);
//	EXPECT_EQ(5, x1.D());
//	quad x2(10);
//	EXPECT_EQ(10, x2.a);
//	EXPECT_EQ(0, x2.b);
//	EXPECT_EQ(5, x2.D());
//	quad x3(+2, -5);
//	EXPECT_EQ(+2, x3.a);
//	EXPECT_EQ(-5, x3.b);
//	EXPECT_EQ(5, x3.D());
//	quad q4(+2, -5, 7); // ignore D if static
//	EXPECT_EQ(+2, q4.a);
//	EXPECT_EQ(-5, q4.b);
//	EXPECT_EQ(5, q4.D());
//	quad q5(q4);
//	EXPECT_EQ(+2, q5.a);
//	EXPECT_EQ(-5, q5.b);
//	EXPECT_EQ(5, q5.D());
//}
//
//TEST(double_int_test, constructor_x) {
//	quadx x1;
//	EXPECT_EQ(0, x1.a);
//	EXPECT_EQ(0, x1.b);
//	EXPECT_EQ(0, x1.D());
//	quadx x2(10);
//	EXPECT_EQ(10, x2.a);
//	EXPECT_EQ(0, x2.b);
//	EXPECT_EQ(0, x2.D());
//	quadx x3(+2, -5);
//	EXPECT_EQ(+2, x3.a);
//	EXPECT_EQ(-5, x3.b);
//	EXPECT_EQ(0, x3.D());
//	quadx q4(+2, -5, 7); // ignore D if static
//	EXPECT_EQ(+2, q4.a);
//	EXPECT_EQ(-5, q4.b);
//	EXPECT_EQ(7, q4.D());
//	quadx q5(q4);
//	EXPECT_EQ(+2, q5.a);
//	EXPECT_EQ(-5, q5.b);
//	EXPECT_EQ(7, q5.D());
//}
//
//template<typename T>
//void test_comparison(bool eq, bool lt, const T& lhs, const T& rhs) {
//	ASSERT_FALSE(eq && lt);
//	EXPECT_EQ(eq, lhs == rhs);
//	EXPECT_EQ(!eq, lhs != rhs);
//	EXPECT_EQ(lt, lhs < rhs);
//	EXPECT_EQ(!(lt || eq), lhs > rhs);
//	EXPECT_EQ((lt || eq), lhs <= rhs);
//	EXPECT_EQ(!lt, lhs >= rhs);
//}
//
//TEST(double_int_test, operators_comparison) {
//	test_comparison(true, false, quad(2, 5), quad(2, 5));
//	test_comparison(false, false, quad(2, 5), quad(2, 3));
//	test_comparison(false, true, quad(2, 5), quad(2, 7));
//	test_comparison(false, true, quad(2, 5), quad(4, 5));
//	test_comparison(false, true, quad(2, 5), quad(4, 3));
//	test_comparison(false, true, quad(2, 5), quad(4, 7));
//	test_comparison(false, false, quad(2, 5), quad(1, 5));
//	test_comparison(false, false, quad(2, 5), quad(1, 3));
//	test_comparison(false, false, quad(2, 5), quad(1, 7));
//}
//
//TEST(double_int_test, operators_arithmetic) {
//	const quad x1(2, -5);
//	const quad x2(3, 4);
//	const quad x3(3, -2);
//	EXPECT_EQ(quad(5, -1), x1 + x2);
//	EXPECT_EQ(quad(-1, -9), x1 - x2);
//	EXPECT_EQ(quad(-2, 5), -x1);
//	EXPECT_EQ(quad(-94, -7), x1 * x2);
//	EXPECT_EQ(quad(4, 1), x1 / x3);
//	EXPECT_EQ(quad(5, -1), x1 % x2);
//	EXPECT_EQ(quad(5, -1), x2 + x1);
//	EXPECT_EQ(quad(1, 9), x2 - x1);
//	EXPECT_EQ(quad(-3, -4), -x2);
//	EXPECT_EQ(quad(-94, -7), x2 * x1);
//	EXPECT_EQ(quad(-6, 15), x1 * -3);
//	EXPECT_EQ(quad(1, -2), x1 / 2);
//}
//
//TEST(double_int_test, operators_inplace) {
//	const quad x1(2, -5);
//	const quad x2(3, 4);
//	const quad x3(3, -2);
//	quad qr;
//	qr = x1; qr += x2;
//	EXPECT_EQ(quad(5, -1), qr);
//	qr = x1; qr -= x2;
//	EXPECT_EQ(quad(-1, -9), qr);
//	qr = x1; qr *= x2;
//	EXPECT_EQ(quad(-94, -7), qr);
//	qr = x1; qr /= x3;
//	EXPECT_EQ(quad(4, 1), qr);
//	qr = x1; qr %= x2;
//	EXPECT_EQ(quad(5, -1), qr);
//	qr = x2; qr += x1;
//	EXPECT_EQ(quad(5, -1), qr);
//	qr = x2; qr -= x1;
//	EXPECT_EQ(quad(1, 9), qr);
//	qr = x2; qr *= x1;
//	EXPECT_EQ(quad(-94, -7), qr);
//	qr = x1; qr *= -3;
//	EXPECT_EQ(quad(-6, 15), qr);
//	qr = x1; qr /= 2;
//	EXPECT_EQ(quad(1, -2), qr);
//}
//
//TEST(double_int_test, operators_inplace_self) {
//	const quad x1(2, -5);
//	quad qr;
//	qr = x1; qr += qr;
//	EXPECT_EQ(quad(4, -10), qr);
//	qr = x1; qr -= qr;
//	EXPECT_EQ(quad(0, 0), qr);
//	qr = x1; qr *= qr;
//	EXPECT_EQ(quad(129, -20), qr);
//	qr = x1; qr /= qr;
//	EXPECT_EQ(quad(1, 0), qr);
//	qr = x1; qr %= qr;
//	EXPECT_EQ(quad(0, 0), qr);
//}
//
//TEST(double_int_test, conjugate) {
//	const quad x1(2, -5);
//	const quad x2(2, 3);
//	EXPECT_EQ(quad(2, 5), x1.conjugate());
//	EXPECT_EQ(quad(2, -3), x2.conjugate());
//}
//
//TEST(double_int_test, norm) {
//	const quad x1(2, -5);
//	const quad x2(3, 4);
//	EXPECT_EQ(-121, x1.norm());
//	EXPECT_EQ(-71, x2.norm());
//	const gaussian g1(2, -5);
//	const gaussian g2(3, 4);
//	EXPECT_EQ(29, g1.norm());
//	EXPECT_EQ(25, g2.norm());
//}
//
//TEST(double_int_test, identity) {
//	const quad q(2, -5);
//	const quad e0 = zeroT<quad>::of(q);
//	const quad e1 = identityT<quad>::of(q);
//	EXPECT_EQ(0, e0.a);
//	EXPECT_EQ(0, e0.b);
//	EXPECT_EQ(5, e0.D());
//	EXPECT_EQ(1, e1.a);
//	EXPECT_EQ(0, e1.b);
//	EXPECT_EQ(5, e1.D());
//}
//
//TEST(double_int_test, identity_x) {
//	const quadx z(2, -5, -1);
//	const quadx z0 = zeroT<quadx>::of(z);
//	const quadx z1 = identityT<quadx>::of(z);
//	EXPECT_EQ(0, z0.a);
//	EXPECT_EQ(0, z0.b);
//	EXPECT_EQ(-1, z0.D());
//	EXPECT_EQ(1, z1.a);
//	EXPECT_EQ(0, z1.b);
//	EXPECT_EQ(-1, z1.D());
//}
