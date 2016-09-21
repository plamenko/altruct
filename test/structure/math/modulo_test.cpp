#include "structure/math/modulo.h"

#include "gtest/gtest.h"

#include <functional>

using namespace std;
using namespace altruct::math;

typedef modulo<int, 1000000000, modulo_storage::CONSTANT> modc;
typedef modulo<int, 1000000007> mod;
typedef moduloX<int> modx;

TEST(modulo_test, constructor) {
	mod m1;
	EXPECT_EQ(0, m1.v);
	mod m2(10);
	EXPECT_EQ(10, m2.v);
	mod m3(m2);
	EXPECT_EQ(10, m3.v);
	mod m4(2000000008);
	EXPECT_EQ(1000000001, m4.v);
	mod m5(-6);
	EXPECT_EQ(1000000001, m5.v);
	mod m6(-1000000013);
	EXPECT_EQ(1000000001, m6.v);
}

TEST(modulo_test, operators_comparison) {
	const mod m1(10);
	const mod m2(20);
	EXPECT_EQ(false, m1 == m2);
	EXPECT_EQ(true, m1 != m2);
	EXPECT_EQ(true, m1 < m2);
	EXPECT_EQ(false, m1 > m2);
	EXPECT_EQ(true, m1 <= m2);
	EXPECT_EQ(false, m1 >= m2);
	EXPECT_EQ(false, m2 == m1);
	EXPECT_EQ(true, m2 != m1);
	EXPECT_EQ(false, m2 < m1);
	EXPECT_EQ(true, m2 > m1);
	EXPECT_EQ(false, m2 <= m1);
	EXPECT_EQ(true, m2 >= m1);
	EXPECT_EQ(true, m2 == m2);
	EXPECT_EQ(false, m2 != m2);
	EXPECT_EQ(false, m2 < m2);
	EXPECT_EQ(false, m2 > m2);
	EXPECT_EQ(true, m2 <= m2);
	EXPECT_EQ(true, m2 >= m2);
}

TEST(modulo_test, operators_arithmetic) {
	const mod m1(1000000000);
	const mod m2(2000000023);
	const mod m3(3000000000LL % mod::M());
	EXPECT_EQ(mod(-7), m1);
	EXPECT_EQ(mod(9), m2);
	EXPECT_EQ(mod(-21), m3);
	EXPECT_EQ(mod(2), m1 + m2);
	EXPECT_EQ(mod(-16), m1 - m2);
	EXPECT_EQ(mod(7), -m1);
	EXPECT_EQ(mod(-63), m1 * m2);
	EXPECT_EQ(mod(222222223), m1 / m2);
	EXPECT_EQ(mod(1), m1 % m2);
	EXPECT_EQ(mod(2), m2 + m1);
	EXPECT_EQ(mod(16), m2 - m1);
	EXPECT_EQ(mod(-9), -m2);
	EXPECT_EQ(mod(-63), m2 * m1);
	EXPECT_EQ(mod(714285718), m2 / m1);
	EXPECT_EQ(mod(9), m2 % m1);
	EXPECT_EQ(mod(3), m3 / m1);
	EXPECT_EQ(mod(333333336), m1 / m3);
}

TEST(modulo_test, operators_inplace) {
	const mod m1(1000000000);
	const mod m2(2000000023);
	const mod m3(3000000000LL % mod::M());
	mod mr;
	mr = m1; mr += m2;
	EXPECT_EQ(mod(2), mr);
	mr = m1; mr -= m2;
	EXPECT_EQ(mod(-16), mr);
	mr = m1; mr *= m2;
	EXPECT_EQ(mod(-63), mr);
	mr = m1; mr /= m2;
	EXPECT_EQ(mod(222222223), mr);
	mr = m1; mr %= m2;
	EXPECT_EQ(mod(1), mr);
	mr = m2; mr += m1;
	EXPECT_EQ(mod(2), mr);
	mr = m2; mr -= m1;
	EXPECT_EQ(mod(16), mr);
	mr = m2; mr *= m1;
	EXPECT_EQ(mod(-63), mr);
	mr = m2; mr /= m1;
	EXPECT_EQ(mod(714285718), mr);
	mr = m2; mr %= m1;
	EXPECT_EQ(mod(9), mr);
	mr = m3; mr /= m1;
	EXPECT_EQ(mod(3), m3 / m1);
	mr = m1; mr /= m3;
	EXPECT_EQ(mod(333333336), m1 / m3);
}

TEST(modulo_test, operators_inplace_self) {
	const mod m1(1000000000);
	mod mr;
	mr = m1; mr += mr;
	EXPECT_EQ(mod(-14), mr);
	mr = m1; mr -= mr;
	EXPECT_EQ(mod(0), mr);
	mr = m1; mr *= mr;
	EXPECT_EQ(mod(49), mr);
	mr = m1; mr /= mr;
	EXPECT_EQ(mod(1), mr);
	mr = m1; mr %= mr;
	EXPECT_EQ(mod(0), mr);
}

TEST(modulo_test, division) {
	// 18 directly divisible by 6
	EXPECT_EQ(modx(3, 1000), modx(18, 1000) / modx(6, 1000));
	EXPECT_EQ(modx(18, 1000), modx(3, 1000) * modx(6, 1000));

	// 7 is invertible modulo 1000
	EXPECT_EQ(modx(430, 1000), modx(10, 1000) / modx(7, 1000));
	EXPECT_EQ(modx(10, 1000), modx(430, 1000) * modx(7, 1000));

	// 48 is not invertible modulo 1000,
	// but after dividing all three (56, 48 and 1000)
	// by their GCD 8,  48/8=6 is invertible modulo 1000/8=125
	EXPECT_EQ(modx(147, 1000), modx(56, 1000) / modx(48, 1000));
	EXPECT_EQ(modx(56, 1000), modx(147, 1000) * modx(48, 1000));

	// 48 is not invertible modulo 1000,
	// and even after dividing all three (28, 48 and 1000)
	// by their GCD 4,  48/4=12 is still not invertible modulo 1000/4=250
	// hence the result is g times bigger where g = gcd(12, 250) = 2
	EXPECT_EQ(modx(147, 1000), modx(28, 1000) / modx(48, 1000));
	EXPECT_EQ(modx(28*2, 1000), modx(147, 1000) * modx(48, 1000));

	EXPECT_EQ(modx(53, 100), modx(17, 100).inv());
}

TEST(modulo_test, identity) {
	mod m1(1000000000);
	mod e0 = zeroT<mod>::of(m1);
	mod e1 = identityT<mod>::of(m1);
	EXPECT_EQ(0, e0.v);
	EXPECT_EQ(1000000007, e0.M());
	EXPECT_EQ(1, e1.v);
	EXPECT_EQ(1000000007, e1.M());
	mod mr = powT(m1, 10);
	EXPECT_EQ(282475249, mr.v);
	EXPECT_EQ(1000000007, mr.M());

}

TEST(modulo_test, int64) {
	typedef modulo<int64_t, 1> modl;
	modl::M() = 1000000000000000003LL;
	const modl m1(1000000000000000000LL);
	const modl m2(2000000000000000008LL);
	const modl m4(4000000000000000000LL);
	EXPECT_EQ(modl(-3), m1);
	EXPECT_EQ(modl(2), m2);
	EXPECT_EQ(modl(-12), m4);
	EXPECT_EQ(modl(-1), m1 + m2);
	EXPECT_EQ(modl(-5), m1 - m2);
	EXPECT_EQ(modl(3), -m1);
	EXPECT_EQ(modl(-6), m1 * m2);
	EXPECT_EQ(modl(500000000000000000LL), m1 / m2);
	EXPECT_EQ(modl(0), m1 % m2);
	EXPECT_EQ(modl(-1), m2 + m1);
	EXPECT_EQ(modl(5), m2 - m1);
	EXPECT_EQ(modl(-2), -m2);
	EXPECT_EQ(modl(-6), m2 * m1);
	EXPECT_EQ(modl(666666666666666668LL), m2 / m1);
	EXPECT_EQ(modl(2), m2 % m1);
	EXPECT_EQ(modl(4), m4 / m1);
	EXPECT_EQ(modl(250000000000000001LL), m1 / m4);
}

template<typename T, typename F>
void modulo_test_perf_impl(T a, T b, int n, const char *msg, const F& func) {
	double clocks_per_sec = 1000;
	auto T0 = clock();
	for (int i = 0; i < n; i++) func(a, b);
	double dT = (clock() - T0) / clocks_per_sec;
	double Mops = n / dT / 1000000;
	printf("%s: %0.2lf Mops  %0.2lf s   %p\n", msg, Mops, dT, &a);
}

TEST(modulo_test, perf) {
	return; // do not run perf test by default

	typedef modulo<int64_t, 1> modl;
	modl::M() = 1000000000000000003LL;
	
	int ni = 1000000000;
	int ai(12345678);
	int bi(456789);
	modulo_test_perf_impl(ai, bi, ni, "int add", [](int &a, int &b){a += b; b++; });
	modulo_test_perf_impl(ai, bi, ni, "int sub", [](int &a, int &b){a -= b; b--; });
	modulo_test_perf_impl(ai, bi, ni, "int neg", [](int &a, int &b){a = -b; b--; });
	modulo_test_perf_impl(ai, bi, ni / 3, "int mul", [](int &a, int &b){a *= b; b++; });
	modulo_test_perf_impl(ai, bi, ni / 10, "int div", [](int &a, int &b){a /= b; a += 1000000000; });

	int nmc = 100000000;
	modc amc(12345678);
	modc bmc(13456789);
	modulo_test_perf_impl(amc, bmc, nmc, "modc<int> add", [](modc &a, modc &b){a += b; b.v++; });
	modulo_test_perf_impl(amc, bmc, nmc, "modc<int> sub", [](modc &a, modc &b){a -= b; b.v++; });
	modulo_test_perf_impl(amc, bmc, nmc, "modc<int> neg", [](modc &a, modc &b){a = -b; b.v++; });
	modulo_test_perf_impl(amc, bmc, nmc / 3, "modc<int> mul", [](modc &a, modc &b){a *= b; b.v++; });
	modulo_test_perf_impl(amc, bmc, nmc / 30, "modc<int> div", [](modc &a, modc &b){a /= b; a.v++; });

	int nmi = 100000000;
	mod ami(12345678);
	mod bmi(13456789);
	modulo_test_perf_impl(ami, bmi, nmi, "mod<int> add", [](mod &a, mod &b){a += b; b.v++; });
	modulo_test_perf_impl(ami, bmi, nmi, "mod<int> sub", [](mod &a, mod &b){a -= b; b.v++; });
	modulo_test_perf_impl(ami, bmi, nmi, "mod<int> neg", [](mod &a, mod &b){a = -b; b.v++; });
	modulo_test_perf_impl(ami, bmi, nmi / 3, "mod<int> mul", [](mod &a, mod &b){a *= b; b.v++; });
	modulo_test_perf_impl(ami, bmi, nmi / 30, "mod<int> div", [](mod &a, mod &b){a /= b; a.v++; });

	int nml = 100000000;
	modl aml(12345678);
	modl bml(13456789);
	modulo_test_perf_impl(aml, bml, nml, "mod<ll> add", [](modl &a, modl &b){a += b; b.v++; });
	modulo_test_perf_impl(aml, bml, nml, "mod<ll> sub", [](modl &a, modl &b){a -= b; b.v++; });
	modulo_test_perf_impl(aml, bml, nml, "mod<ll> neg", [](modl &a, modl &b){a = -b; b.v++; });
	modulo_test_perf_impl(aml, bml, nml / 100, "mod<ll> mul", [](modl &a, modl &b){a *= b; b.v++; });
	modulo_test_perf_impl(aml, bml, nml / 300, "mod<ll> div", [](modl &a, modl &b){a /= b; a.v++; });
}
