#include "altruct/structure/math/polynom.h"
#include "altruct/structure/math/modulo.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

typedef modulo<int, 1000000007> mod;
typedef polynom<mod> poly;
typedef modulo<poly, 1> polymod;

TEST(polymod_test, constructor) {
	polymod::M() = poly{ 0, 0, 0, 0, 1 };
	polymod p0;
	EXPECT_EQ((poly{}), p0.v);
	polymod p1(7);
	EXPECT_EQ((poly{ 7 }), p1.v);
}

TEST(polymod_test, division) {
	// irreducible polynomial M(x) = x^2 - x^1 - x^0
	polymod::M() = poly{ -1, -1, 1 };
	const polymod x(poly{ 0, 1 });
	const polymod x20 = powT(x, 20);
	const polymod x100 = powT(x, 100);
	const polymod x120 = powT(x, 120);
	EXPECT_EQ(x120, x20 * x100);
	EXPECT_EQ(x120, x100 * x20);
	EXPECT_EQ(x100, x120 / x20);
	EXPECT_EQ(x20, x120 / x100);
}

TEST(polymod_test, fibonacci) {
	// f(n+2) - f(n+1) - f(n) = 0;
	// p(x) = x^2 - x^1 - x^0
	// f(n) = x^n % p(x)
	polymod::M() = poly{ -1, -1, 1 };
	polymod x(poly{0, 1});
	vector<mod> vf;
	for (int i = 0; i < 13; i++) {
		vf.push_back(powT(x, i).v[1]);
	}
	EXPECT_EQ((vector<mod>{0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144}), vf);
	EXPECT_EQ(mod(687995182), powT(x, 100).v[1]); // f(100) % 1000000007
}
