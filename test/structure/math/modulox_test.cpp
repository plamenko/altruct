#include "structure/math/modulo.h"

#include "gtest/gtest.h"

#include <functional>

using namespace std;
using namespace altruct::math;

typedef moduloX<int> modx;

TEST(modulox_test, constructor) {
	modx m1;
	EXPECT_EQ(0, m1.v);
	EXPECT_EQ(1, m1.M());
	modx m2(10);
	EXPECT_EQ(10, m2.v);
	EXPECT_EQ(1, m2.M());
	modx m3(m2);
	EXPECT_EQ(10, m3.v);
	EXPECT_EQ(1, m3.M());
	modx m4(2000000008, 1000000007);
	EXPECT_EQ(1000000001, m4.v);
	EXPECT_EQ(1000000007, m4.M());
	modx m5(-6, 1000000011);
	EXPECT_EQ(1000000005, m5.v);
	EXPECT_EQ(1000000011, m5.M());
	modx m6(-1000000013, 1000000007);
	EXPECT_EQ(1000000001, m6.v);
	EXPECT_EQ(1000000007, m6.M());
	modx m7(m6);
	EXPECT_EQ(1000000001, m7.v);
	EXPECT_EQ(1000000007, m7.M());
}

TEST(modulox_test, operators_comparison) {
	const modx m1(10, 1000000007);
	const modx m2(20, 1000000007);
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

TEST(modulox_test, operators_arithmetic) {
	const modx m1(1000000000, 1000000007);
	const modx m2(2000000023, 1000000007);
	const modx m3(999999986, 1000000007);
	EXPECT_EQ(modx(-7, 1000000007), m1);
	EXPECT_EQ(modx(9, 1000000007), m2);
	EXPECT_EQ(modx(-21, 1000000007), m3);
	EXPECT_EQ(modx(2, 1000000007), m1 + m2);
	EXPECT_EQ(modx(-16, 1000000007), m1 - m2);
	EXPECT_EQ(modx(7, 1000000007), -m1);
	EXPECT_EQ(modx(-63, 1000000007), m1 * m2);
	EXPECT_EQ(modx(222222223, 1000000007), m1 / m2);
	EXPECT_EQ(modx(1, 1000000007), m1 % m2);
	EXPECT_EQ(modx(2, 1000000007), m2 + m1);
	EXPECT_EQ(modx(16, 1000000007), m2 - m1);
	EXPECT_EQ(modx(-9, 1000000007), -m2);
	EXPECT_EQ(modx(-63, 1000000007), m2 * m1);
	EXPECT_EQ(modx(714285718, 1000000007), m2 / m1);
	EXPECT_EQ(modx(9, 1000000007), m2 % m1);
	EXPECT_EQ(modx(3, 1000000007), m3 / m1);
	EXPECT_EQ(modx(333333336, 1000000007), m1 / m3);
}

TEST(modulox_test, operators_inplace) {
	const modx m1(1000000000, 1000000007);
	const modx m2(2000000023, 1000000007);
	const modx m3(999999986, 1000000007);
	modx mr;
	mr = m1; mr += m2;
	EXPECT_EQ(modx(2, 1000000007), mr);
	mr = m1; mr -= m2;
	EXPECT_EQ(modx(-16, 1000000007), mr);
	mr = m1; mr *= m2;
	EXPECT_EQ(modx(-63, 1000000007), mr);
	mr = m1; mr /= m2;
	EXPECT_EQ(modx(222222223, 1000000007), mr);
	mr = m1; mr %= m2;
	EXPECT_EQ(modx(1, 1000000007), mr);
	mr = m2; mr += m1;
	EXPECT_EQ(modx(2, 1000000007), mr);
	mr = m2; mr -= m1;
	EXPECT_EQ(modx(16, 1000000007), mr);
	mr = m2; mr *= m1;
	EXPECT_EQ(modx(-63, 1000000007), mr);
	mr = m2; mr /= m1;
	EXPECT_EQ(modx(714285718, 1000000007), mr);
	mr = m2; mr %= m1;
	EXPECT_EQ(modx(9, 1000000007), mr);
	mr = m3; mr /= m1;
	EXPECT_EQ(modx(3, 1000000007), m3 / m1);
	mr = m1; mr /= m3;
	EXPECT_EQ(modx(333333336, 1000000007), m1 / m3);
}

TEST(modulox_test, operators_inplace_self) {
	const modx m1(1000000000, 1000000007);
	modx mr;
	mr = m1; mr += mr;
	EXPECT_EQ(modx(-14, 1000000007), mr);
	mr = m1; mr -= mr;
	EXPECT_EQ(modx(0, 1000000007), mr);
	mr = m1; mr *= mr;
	EXPECT_EQ(modx(49, 1000000007), mr);
	mr = m1; mr /= mr;
	EXPECT_EQ(modx(1, 1000000007), mr);
	mr = m1; mr %= mr;
	EXPECT_EQ(modx(0, 1000000007), mr);
}

TEST(modulox_test, identity) {
	modx m1(1000000000, 1000000007);
	modx e0 = zeroT<modx>::of(m1);
	modx e1 = identityT<modx>::of(m1);
	EXPECT_EQ(0, e0.v);
	EXPECT_EQ(1000000007, e0.M());
	EXPECT_EQ(1, e1.v);
	EXPECT_EQ(1000000007, e1.M());
	modx mr = powT(m1, 10);
	EXPECT_EQ(282475249, mr.v);
	EXPECT_EQ(1000000007, mr.M());

}

TEST(modulox_test, int64) {
	typedef moduloX<int64_t> modxl;
	const modxl m1(1000000000000000000LL, 1000000000000000003LL);
	const modxl m2(2000000000000000008LL, 1000000000000000003LL);
	const modxl m4(4000000000000000000LL, 1000000000000000003LL);
	EXPECT_EQ(modxl(-3, 1000000000000000003LL), m1);
	EXPECT_EQ(modxl(2, 1000000000000000003LL), m2);
	EXPECT_EQ(modxl(-12, 1000000000000000003LL), m4);
	EXPECT_EQ(modxl(-1, 1000000000000000003LL), m1 + m2);
	EXPECT_EQ(modxl(-5, 1000000000000000003LL), m1 - m2);
	EXPECT_EQ(modxl(3, 1000000000000000003LL), -m1);
	EXPECT_EQ(modxl(-6, 1000000000000000003LL), m1 * m2);
	EXPECT_EQ(modxl(500000000000000000LL, 1000000000000000003LL), m1 / m2);
	EXPECT_EQ(modxl(0, 1000000000000000003LL), m1 % m2);
	EXPECT_EQ(modxl(-1, 1000000000000000003LL), m2 + m1);
	EXPECT_EQ(modxl(5, 1000000000000000003LL), m2 - m1);
	EXPECT_EQ(modxl(-2, 1000000000000000003LL), -m2);
	EXPECT_EQ(modxl(-6, 1000000000000000003LL), m2 * m1);
	EXPECT_EQ(modxl(666666666666666668LL, 1000000000000000003LL), m2 / m1);
	EXPECT_EQ(modxl(2, 1000000000000000003LL), m2 % m1);
	EXPECT_EQ(modxl(4, 1000000000000000003LL), m4 / m1);
	EXPECT_EQ(modxl(250000000000000001LL, 1000000000000000003LL), m1 / m4);
}
