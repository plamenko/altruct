#include "algorithm/math/totient_sums.h"
#include "algorithm/math/ranges.h"
#include "structure/math/modulo.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

TEST(totient_sums_test, sum_phi_D_L) {
	typedef modulo<int, 1000000007> mod;

	auto id = mod(1);
	auto castT = [](int64_t n){ return mod(n % 1000000007); };
	auto vn = range<int64_t>(21);
	
	EXPECT_EQ((vector<mod>{0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}), sum_phi_D_L(0, 0, vn, id, castT));
	EXPECT_EQ((vector<mod>{0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}), sum_phi_D_L(0, 1, vn, id, castT));
	EXPECT_EQ((vector<mod>{0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}), sum_phi_D_L(0, 2, vn, id, castT));

	EXPECT_EQ((vector<mod>{0, 1, 2, 4, 6, 10, 12, 18, 22, 28, 32, 42, 46, 58, 64, 72, 80, 96, 102, 120, 128}), sum_phi_D_L(1, 0, vn, id, castT));
	EXPECT_EQ((vector<mod>{0, 1, 3, 9, 17, 37, 49, 91, 123, 177, 217, 327, 375, 531, 615, 735, 863, 1135, 1243, 1585, 1745}), sum_phi_D_L(1, 1, vn, id, castT));
	EXPECT_EQ((vector<mod>{0, 1, 5, 23, 55, 155, 227, 521, 777, 1263, 1663, 2873, 3449, 5477, 6653, 8453, 10501, 15125, 17069, 23567, 26767}), sum_phi_D_L(1, 2, vn, id, castT));

	EXPECT_EQ((vector<mod>{0, 1, 3, 8, 15, 29, 42, 69, 95, 134, 172, 237, 287, 377, 452, 552, 652, 804, 915, 1104, 1252}), sum_phi_D_L(2, 0, vn, id, castT));
	EXPECT_EQ((vector<mod>{0, 1, 5, 20, 48, 118, 196, 385, 593, 944, 1324, 2039, 2639, 3809, 4859, 6359, 7959, 10543, 12541, 16132, 19092}), sum_phi_D_L(2, 1, vn, id, castT));
	EXPECT_EQ((vector<mod>{0, 1, 9, 54, 166, 516, 984, 2307, 3971, 7130, 10930, 18795, 25995, 41205, 55905, 78405, 104005, 147933, 183897, 252126, 311326}), sum_phi_D_L(2, 2, vn, id, castT));
}

TEST(totient_sums_test, sum_phi_D_L_modx) {
	typedef moduloX<int> modx;

	auto id = modx(1, 1009);
	auto castT = [](int64_t n){ return modx(n % 1009, 1009); };
	auto vn = range<int64_t>(21);
	
	EXPECT_EQ((vector<modx>{0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}), sum_phi_D_L(0, 0, vn, id, castT));
	EXPECT_EQ((vector<modx>{0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}), sum_phi_D_L(0, 1, vn, id, castT));
	EXPECT_EQ((vector<modx>{0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}), sum_phi_D_L(0, 2, vn, id, castT));

	EXPECT_EQ((vector<modx>{0, 1, 2, 4, 6, 10, 12, 18, 22, 28, 32, 42, 46, 58, 64, 72, 80, 96, 102, 120, 128}), sum_phi_D_L(1, 0, vn, id, castT));
	EXPECT_EQ((vector<modx>{0, 1, 3, 9, 17, 37, 49, 91, 123, 177, 217, 327, 375, 531, 615, 735, 863, 126, 234, 576, 736}), sum_phi_D_L(1, 1, vn, id, castT));
	EXPECT_EQ((vector<modx>{0, 1, 5, 23, 55, 155, 227, 521, 777, 254, 654, 855, 422, 432, 599, 381, 411, 999, 925, 360, 533}), sum_phi_D_L(1, 2, vn, id, castT));

	EXPECT_EQ((vector<modx>{0, 1, 3, 8, 15, 29, 42, 69, 95, 134, 172, 237, 287, 377, 452, 552, 652, 804, 915, 95, 243}), sum_phi_D_L(2, 0, vn, id, castT));
	EXPECT_EQ((vector<modx>{0, 1, 5, 20, 48, 118, 196, 385, 593, 944, 315, 21, 621, 782, 823, 305, 896, 453, 433, 997, 930}), sum_phi_D_L(2, 1, vn, id, castT));
	EXPECT_EQ((vector<modx>{0, 1, 9, 54, 166, 516, 984, 289, 944, 67, 840, 633, 770, 845, 410, 712, 78, 619, 259, 885, 554}), sum_phi_D_L(2, 2, vn, id, castT));
}
