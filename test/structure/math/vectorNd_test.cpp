#include "altruct/structure/math/vectorNd.h"
#include "altruct/algorithm/math/modulos.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

typedef vectorNd<int, 4> vect4;
typedef vectorNd<double, 5> vect5;

TEST(vectorNd_test, constructor) {
	vect4 v1;
	EXPECT_EQ(4, v1.size());
	EXPECT_EQ((array<int, 4>{{ 0, 0, 0, 0 }}), v1.a);
	vect4 v2(1);
	EXPECT_EQ(4, v2.size());
	EXPECT_EQ((array<int, 4>{{ 1, 1, 1, 1 }}), v2.a);
	vect4 v3(v2);
	EXPECT_EQ(4, v3.size());
	EXPECT_EQ((array<int, 4>{{ 1, 1, 1, 1 }}), v3.a);
	vect5 v4(7);
	EXPECT_EQ(5, v4.size());
	EXPECT_EQ((array<double, 5>{{ 7, 7, 7, 7, 7 }}), v4.a);
	vect5 v5(7.5);
	EXPECT_EQ(5, v5.size());
	EXPECT_EQ((array<double, 5>{{ 7.5, 7.5, 7.5, 7.5, 7.5 }}), v5.a);
	vect4 v6(vect4{ 1, 2, 3, 4 });
	EXPECT_EQ(4, v6.size());
	EXPECT_EQ((array<int, 4>{{ 1, 2, 3, 4 }}), v6.a);
	vect4 v7(array<int, 4>{{ 1, 2, 3, 4 }});
	EXPECT_EQ(4, v7.size());
	EXPECT_EQ((array<int, 4>{{ 1, 2, 3, 4 }}), v7.a);
	vect4 v8{ 1, 2, 3, 4 };
	EXPECT_EQ(4, v8.size());
	EXPECT_EQ((array<int, 4>{{ 1, 2, 3, 4 }}), v8.a);
}

TEST(vectorNd_test, brackets) {
	const vect4 v1{ 2, 3, 5, 7 };
	EXPECT_EQ(2, v1[0]);
	EXPECT_EQ(3, v1[1]);
	EXPECT_EQ(5, v1[2]);
	EXPECT_EQ(7, v1[3]);
	vect4 v2{ 1, 2, 3, 4 };
	for (int i = 0; i < v2.size(); i++) v2[i] = -v2[i];
	EXPECT_EQ(-1, v2[0]);
	EXPECT_EQ(-2, v2[1]);
	EXPECT_EQ(-3, v2[2]);
	EXPECT_EQ(-4, v2[3]);
}

template<typename T>
void test_comparison(bool eq, bool lt, const T& lhs, const T& rhs) {
	ASSERT_FALSE(eq && lt);
	EXPECT_EQ(eq, lhs == rhs);
	EXPECT_EQ(!eq, lhs != rhs);
	EXPECT_EQ(lt, lhs < rhs);
	EXPECT_EQ(!(lt || eq), lhs > rhs);
	EXPECT_EQ((lt || eq), lhs <= rhs);
	EXPECT_EQ(!lt, lhs >= rhs);
}

TEST(vectNd_test, operators_comparison) {
	test_comparison(true, false, vect4{ 1, 2, 3, 4 }, vect4{ 1, 2, 3, 4 });
	test_comparison(false, false, vect4{ 1, 2, 3, 4 }, vect4{ 1, 2, 0, 4 });
	test_comparison(false, true, vect4{ 1, 2, 3, 4 }, vect4{ 1, 5, 3, 4 });
	test_comparison(false, true, vect4{ 1, 2, 3, 4 }, vect4{ 3, 2, 3, 4 });
	test_comparison(false, true, vect4{ 1, 2, 3, 4 }, vect4{ 3, 2, 0, 4 });
	test_comparison(false, true, vect4{ 1, 2, 3, 4 }, vect4{ 3, 5, 3, 4 });
	test_comparison(false, false, vect4{ 1, 2, 3, 4 }, vect4{ 0, 2, 3, 4 });
	test_comparison(false, false, vect4{ 1, 2, 3, 4 }, vect4{ 0, 2, 0, 4 });
	test_comparison(false, false, vect4{ 1, 2, 3, 4 }, vect4{ 0, 5, 3, 4 });
}

TEST(vectNd_test, operators_arithmetic) {
	const vect4 v1{ 2, -5, 3, 16 };
	const vect4 v2{ 3, 10, 12, 7 };
	const vect4 v3{ 1, 5, 4, 7 };
	EXPECT_EQ((vect4{ 5, 5, 15, 23 }), v1 + v2);
	EXPECT_EQ((vect4{ -1, -15, -9, 9 }), v1 - v2);
	EXPECT_EQ((vect4{ -2, 5, -3, -16 }), -v1);
	EXPECT_EQ((vect4{ 6, -50, 36, 112 }), v1 * v2);
	EXPECT_EQ((vect4{ 3, 2, 3, 1 }), v2 / v3);
	EXPECT_EQ((vect4{ 1, 0, 0, 7 }), v2 % v1);
	EXPECT_EQ((vect4{ 5, -2, 6, 19 }), v1 + 3);
	EXPECT_EQ((vect4{ -1, -8, 0, 13 }), v1 - 3);
	EXPECT_EQ((vect4{ -6, 15, -9, -48 }), v1 * -3);
	EXPECT_EQ((vect4{ 1, -2, 1, 8 }), v1 / 2);
	EXPECT_EQ((vect4{ 0, -1, 1, 0 }), v1 % 2);
}

TEST(vectNd_test, operators_inplace) {
	const vect4 v1{ 2, -5, 3, 16 };
	const vect4 v2{ 3, 10, 12, 7 };
	const vect4 v3{ 1, 5, 4, 7 };
	vect4 vr;
	vr = v1; vr += v2;
	EXPECT_EQ((vect4{ 5, 5, 15, 23 }), vr);
	vr = v1; vr -= v2;
	EXPECT_EQ((vect4{ -1, -15, -9, 9 }), vr);
	vr = v1; vr *= v2;
	EXPECT_EQ((vect4{ 6, -50, 36, 112 }), vr);
	vr = v2; vr /= v3;
	EXPECT_EQ((vect4{ 3, 2, 3, 1 }), vr);
	vr = v2; vr %= v1;
	EXPECT_EQ((vect4{ 1, 0, 0, 7 }), vr);
	vr = v1; vr += 3;
	EXPECT_EQ((vect4{ 5, -2, 6, 19 }), vr);
	vr = v1; vr -= 3;
	EXPECT_EQ((vect4{ -1, -8, 0, 13 }), vr);
	vr = v1; vr *= -3;
	EXPECT_EQ((vect4{ -6, 15, -9, -48 }), vr);
	vr = v1; vr /= 2;
	EXPECT_EQ((vect4{ 1, -2, 1, 8 }), vr);
	vr = v1; vr %= 2;
	EXPECT_EQ((vect4{ 0, -1, 1, 0 }), vr);
}

TEST(vectNd_test, operators_inplace_self) {
	const vect4 v1{ 1, 5, -4, 7 };
	vect4 vr;
	vr = v1; vr += vr;
	EXPECT_EQ((vect4{ 2, 10, -8, 14 }), vr);
	vr = v1; vr -= vr;
	EXPECT_EQ((vect4{ 0, 0, 0, 0 }), vr);
	vr = v1; vr *= vr;
	EXPECT_EQ((vect4{ 1, 25, 16, 49 }), vr);
	vr = v1; vr /= vr;
	EXPECT_EQ((vect4{ 1, 1, 1, 1 }), vr);
	vr = v1; vr %= vr;
	EXPECT_EQ((vect4{ 0, 0, 0, 0 }), vr);
}

TEST(vectNd_test, abs2) {
	const vect4 v0{ 0, 0, 0, 0 };
	EXPECT_EQ(0, v0.abs2());
	const vect4 v1{ 1, 5, -4, 7 };
	EXPECT_EQ(91, v1.abs2());
}

TEST(vectNd_test, identity) {
	typedef moduloX<int> modx;
	typedef vectorNd<modx, 3> vect3;
	vect3 v{ { 2, 1009 }, { 3, 1013 }, { 5, 1019 } };
	vect3 e0 = zeroT<vect3>::of(v);
	EXPECT_EQ(0, e0[0].v);
	EXPECT_EQ(1009, e0[0].M());
	EXPECT_EQ(0, e0[1].v);
	EXPECT_EQ(1013, e0[1].M());
	EXPECT_EQ(0, e0[2].v);
	EXPECT_EQ(1019, e0[2].M());
	vect3 e1 = identityT<vect3>::of(v);
	EXPECT_EQ(1, e1[0].v);
	EXPECT_EQ(1009, e1[0].M());
	EXPECT_EQ(1, e1[1].v);
	EXPECT_EQ(1013, e1[1].M());
	EXPECT_EQ(1, e1[2].v);
	EXPECT_EQ(1019, e1[2].M());
	
	vect3 v1 = e1 * modx(1000);
	vect3 v3 = powT(v1, 3);
	modx r;
	for (int i = 0; i < v3.size(); i++) {
		chinese_remainder(&r.v, &r.M(), v3[i].v, v3[i].M());
	}
	EXPECT_EQ(1000000000, r.v);
	EXPECT_EQ(1009 * 1013 * 1019, r.M());
}
