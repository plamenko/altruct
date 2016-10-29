#include "algorithm/math/totient_sums.h"
#include "algorithm/math/ranges.h"
#include "structure/math/modulo.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;
using namespace altruct::container;

typedef modulo<int, 1000000007> field;

TEST(totient_sums_test, dirichlet_convolution) {
	int n = 21;
	typedef moduloX<int> modx;
	vector<int> vmu(n); moebius_mu(vmu.data(), (int)vmu.size());
	auto id = [&](int n){ return modx(n, 1009); };
	auto mu = [&](int n){ return modx(vmu[n], 1009); };
	vector<modx> phi(n); dirichlet_convolution<modx>(n, id, mu, phi);
	EXPECT_EQ((vector<modx>{0, 1, 1, 2, 2, 4, 2, 6, 4, 6, 4, 10, 4, 12, 6, 8, 8, 16, 6, 18, 8}), phi);
}

TEST(totient_sums_test, dirichlet_inverse) {
	int n = 21;
	typedef moduloX<int> modx;
	auto f = [](int n){ return modx(n * (n + 2), 1009); };
	vector<modx> f_inv(n); dirichlet_inverse<modx>(n, f, f_inv);
	EXPECT_EQ((vector<modx>{0, 673, 896, 671, 635, 893, 452, 1002, 435, 670, 269, 881, 113, 651, 573, 459, 441, 861, 678, 292, 861}), f_inv);
}

TEST(totient_sums_test, sieve_m) {
	int n = 21;
	auto t = [](int n){ return n * (n + 1) / 2; };
	vector<int> actual1(n); sieve_m<int>(n, t, actual1);
	EXPECT_EQ((vector<int>{0, 1, 2, 4, 6, 10, 12, 18, 22, 28, 32, 42, 46, 58, 64, 72, 80, 96, 102, 120, 128}), actual1);

	typedef moduloX<int> modx;
	auto t2 = [](int n){ return modx(n * (n + 1) / 2, 1009); };
	auto p2 = [](int n){ return modx(n + 2, 1009); };
	vector<modx> actual2(n); sieve_m<modx>(n, t2, p2, actual2);
	EXPECT_EQ((vector<modx>{0, 673, 449, 1, 973, 77, 264, 938, 540, 840, 205, 992, 170, 509, 61, 809, 482, 934, 112, 116, 490}), actual2);
}

TEST(totient_sums_test, mertens) {
	int n = 30;
	// preprocess `U = n^(2/3)` values of `Sum[p(k) * f[k], {k, 1, U}]`
	int U = (int)isq(icbrt(n));
	sqrt_map<int, field> mm(U, n);
	vector<int> mu(U); moebius_mu(mu.data(), U);
	for (int k = 1; k < U; k++) mm[k] = mm[k - 1] + mu[k];

	vector<field> va;
	for (int k = 0; k <= n; k++) {
		mm.reset_max(k);
		va.push_back(mertens(k, mm, field(1)));
	}
	EXPECT_EQ((vector<field>{0, 1, 0, -1, -1, -2, -1, -2, -2, -2, -1, -2, -2, -3, -2, -1, -1, -2, -2, -3, -3, -2, -1, -2, -2, -2, -1, -1, -1, -2, -3}), va);
}

TEST(totient_totient_sums_test, sum_phi_D_L) {
	auto id = field(1);
	auto castT = [](int64_t n){ return field(n % 1000000007); };
	auto vn = range<int64_t>(21);
	
	EXPECT_EQ((vector<field>{0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}), sum_phi_D_L(0, 0, vn, 0, id, castT));
	EXPECT_EQ((vector<field>{0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}), sum_phi_D_L(0, 1, vn, 0, id, castT));
	EXPECT_EQ((vector<field>{0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}), sum_phi_D_L(0, 2, vn, 0, id, castT));

	EXPECT_EQ((vector<field>{0, 1, 2, 4, 6, 10, 12, 18, 22, 28, 32, 42, 46, 58, 64, 72, 80, 96, 102, 120, 128}), sum_phi_D_L(1, 0, vn, 0, id, castT));
	EXPECT_EQ((vector<field>{0, 1, 3, 9, 17, 37, 49, 91, 123, 177, 217, 327, 375, 531, 615, 735, 863, 1135, 1243, 1585, 1745}), sum_phi_D_L(1, 1, vn, 0, id, castT));
	EXPECT_EQ((vector<field>{0, 1, 5, 23, 55, 155, 227, 521, 777, 1263, 1663, 2873, 3449, 5477, 6653, 8453, 10501, 15125, 17069, 23567, 26767}), sum_phi_D_L(1, 2, vn, 0, id, castT));

	EXPECT_EQ((vector<field>{0, 1, 3, 8, 15, 29, 42, 69, 95, 134, 172, 237, 287, 377, 452, 552, 652, 804, 915, 1104, 1252}), sum_phi_D_L(2, 0, vn, 0, id, castT));
	EXPECT_EQ((vector<field>{0, 1, 5, 20, 48, 118, 196, 385, 593, 944, 1324, 2039, 2639, 3809, 4859, 6359, 7959, 10543, 12541, 16132, 19092}), sum_phi_D_L(2, 1, vn, 0, id, castT));
	EXPECT_EQ((vector<field>{0, 1, 9, 54, 166, 516, 984, 2307, 3971, 7130, 10930, 18795, 25995, 41205, 55905, 78405, 104005, 147933, 183897, 252126, 311326}), sum_phi_D_L(2, 2, vn, 0, id, castT));

	EXPECT_EQ(field(356214470), sum_phi_D_L(1, 0, 10000000, 0, id, castT));
}

TEST(totient_totient_sums_test, sum_phi_D_L_modx) {
	typedef moduloX<int> modx;

	auto id = modx(1, 1009);
	auto castT = [](int64_t n){ return modx(n % 1009, 1009); };
	auto vn = range<int64_t>(21);
	
	EXPECT_EQ((vector<modx>{0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}), sum_phi_D_L(0, 0, vn, 0, id, castT));
	EXPECT_EQ((vector<modx>{0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}), sum_phi_D_L(0, 1, vn, 0, id, castT));
	EXPECT_EQ((vector<modx>{0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}), sum_phi_D_L(0, 2, vn, 0, id, castT));

	EXPECT_EQ((vector<modx>{0, 1, 2, 4, 6, 10, 12, 18, 22, 28, 32, 42, 46, 58, 64, 72, 80, 96, 102, 120, 128}), sum_phi_D_L(1, 0, vn, 0, id, castT));
	EXPECT_EQ((vector<modx>{0, 1, 3, 9, 17, 37, 49, 91, 123, 177, 217, 327, 375, 531, 615, 735, 863, 126, 234, 576, 736}), sum_phi_D_L(1, 1, vn, 0, id, castT));
	EXPECT_EQ((vector<modx>{0, 1, 5, 23, 55, 155, 227, 521, 777, 254, 654, 855, 422, 432, 599, 381, 411, 999, 925, 360, 533}), sum_phi_D_L(1, 2, vn, 0, id, castT));

	EXPECT_EQ((vector<modx>{0, 1, 3, 8, 15, 29, 42, 69, 95, 134, 172, 237, 287, 377, 452, 552, 652, 804, 915, 95, 243}), sum_phi_D_L(2, 0, vn, 0, id, castT));
	EXPECT_EQ((vector<modx>{0, 1, 5, 20, 48, 118, 196, 385, 593, 944, 315, 21, 621, 782, 823, 305, 896, 453, 433, 997, 930}), sum_phi_D_L(2, 1, vn, 0, id, castT));
	EXPECT_EQ((vector<modx>{0, 1, 9, 54, 166, 516, 984, 289, 944, 67, 840, 633, 770, 845, 410, 712, 78, 619, 259, 885, 554}), sum_phi_D_L(2, 2, vn, 0, id, castT));

	EXPECT_EQ(modx(984, 1009), sum_phi_D_L(1, 0, 10000000, 0, id, castT));
}

TEST(totient_sums_test, sum_primes) {
	vector<int> vp(isqrt(1030) + 1);
	int m = primes(vp.data(), nullptr, (int)vp.size());
	vector<field> va1, va2;
	for (int n = 0; n < 30; n++) {
		va1.push_back(sum_primes(n, vp.data(), field(1)));
		va2.push_back(sum_primes(1000 + n, vp.data(), field(1)) - 76127);
	}
	EXPECT_EQ((vector<field>{0, 0, 2, 5, 5, 10, 10, 17, 17, 17, 17, 28, 28, 41, 41, 41, 41, 58, 58, 77, 77, 77, 77, 100, 100, 100, 100, 100, 100, 129}), va1);
	EXPECT_EQ((vector<field>{0, 0, 0, 0, 0, 0, 0, 0, 0, 1009, 1009, 1009, 1009, 2022, 2022, 2022, 2022, 2022, 2022, 3041, 3041, 4062, 4062, 4062, 4062, 4062, 4062, 4062, 4062, 4062}), va2);
}

TEST(totient_sums_test, sum_primes2) {
	vector<int> vp(1000);
	vector<char> vq(1000);
	int m = primes(vp.data(), vq.data(), (int)vq.size());
	vector<int> ve, va;
	int c = 0;
	for (int n = 0; n < vp.size(); n++) {
		ve.push_back(c += n * vq[n]);
		va.push_back(sum_primes(n, vp.data(), 1));
	}
	EXPECT_EQ(ve, va);
}
