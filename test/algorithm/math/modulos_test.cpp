#include "algorithm/math/modulos.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

TEST(modulos_test, jacobi) {
	vector<vector<int>> v(1);
	for (int m = 1; m <= 50; m++) {
		v.push_back(vector<int>());
		for (int n = 0; n <= 20; n++) {
			v[m].push_back(jacobi(n, m));
		}
	}
	EXPECT_EQ((vector<int>{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}), v[1]);
	EXPECT_EQ((vector<int>{0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1}), v[3]);
	EXPECT_EQ((vector<int>{0, 1, -1, -1, 1, 0, 1, -1, -1, 1, 0, 1, -1, -1, 1, 0, 1, -1, -1, 1, 0}), v[5]);
	EXPECT_EQ((vector<int>{0, 1, 1, -1, 1, -1, -1, 0, 1, 1, -1, 1, -1, -1, 0, 1, 1, -1, 1, -1, -1}), v[7]);
	EXPECT_EQ((vector<int>{0, 1, -1, 0, 1, 0, 0, -1, -1, 0, 0, 1, 0, -1, 1, 0, 1, -1, 0, 1, 0}), v[45]);
}

TEST(modulos_test, sqrt_cipolla) {
	EXPECT_EQ(0, sqrt_cipolla(0, 17));
	EXPECT_EQ(1, sqrt_cipolla(1, 17));
	EXPECT_EQ(6, sqrt_cipolla(2, 17));
	EXPECT_EQ(15, sqrt_cipolla(4, 17));
	EXPECT_EQ(12, sqrt_cipolla(8, 17));
	EXPECT_EQ(14, sqrt_cipolla(9, 17));
	EXPECT_EQ(8, sqrt_cipolla(13, 17));
	EXPECT_EQ(7, sqrt_cipolla(15, 17));
	EXPECT_EQ(13, sqrt_cipolla(16, 17));
	typedef modulo<int, 17> mod;
	EXPECT_EQ(mod(0), sqT(sqrt_cipolla(mod(0))));
	EXPECT_EQ(mod(1), sqT(sqrt_cipolla(mod(1))));
	EXPECT_EQ(mod(2), sqT(sqrt_cipolla(mod(2))));
	EXPECT_EQ(mod(4), sqT(sqrt_cipolla(mod(4))));
	EXPECT_EQ(mod(8), sqT(sqrt_cipolla(mod(8))));
	EXPECT_EQ(mod(9), sqT(sqrt_cipolla(mod(9))));
	EXPECT_EQ(mod(13), sqT(sqrt_cipolla(mod(13))));
	EXPECT_EQ(mod(15), sqT(sqrt_cipolla(mod(15))));
	EXPECT_EQ(mod(16), sqT(sqrt_cipolla(mod(16))));
}

TEST(modulos_test, sqrt_hensel_lift) {
	EXPECT_EQ(0, sqrt_hensel_lift(0, 17, 5));
	EXPECT_EQ(1, sqrt_hensel_lift(1, 17, 5));
	EXPECT_EQ(461199, sqrt_hensel_lift(2, 17, 5));
	EXPECT_EQ(1419855, sqrt_hensel_lift(4, 17, 5));
	EXPECT_EQ(922398, sqrt_hensel_lift(8, 17, 5));
	EXPECT_EQ(1419854, sqrt_hensel_lift(9, 17, 5));
	EXPECT_EQ(499740, sqrt_hensel_lift(13, 17, 5));
	EXPECT_EQ(1318629, sqrt_hensel_lift(15, 17, 5));
	EXPECT_EQ(1419853, sqrt_hensel_lift(16, 17, 5));
	EXPECT_EQ(883131, sqrt_hensel_lift(12346, 17, 5));
}
