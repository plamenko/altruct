#include "algorithm/math/sums.h"
#include "algorithm/math/primes.h"
#include "algorithm/math/polynoms.h"
#include "structure/math/polynom.h"
#include <unordered_map>

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

typedef modulo<int, 1000000007> field;
typedef polynom<field> poly;

TEST(sums_test, sum) {
	auto f = [](int k) { return k*k; };
	EXPECT_EQ(0, (sum<int>(f, 1, 0)));
	EXPECT_EQ(1, (sum<int>(f, 1, 1)));
	EXPECT_EQ(5, (sum<int>(f, 1, 2)));
	EXPECT_EQ(385, (sum<int>(f, 1, 10)));
	EXPECT_EQ(355, (sum<int>(f, 5, 10)));
}

vector<field> calc_sum_pow(int p, int n) {
	vector<field> v;
	for (int k = 0; k <= n; k++) {
		v.push_back(sum_pow<field>(p, k));
	}
	return v;
	
}
TEST(sums_test, sum_pow) {
	EXPECT_EQ((vector<field>{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }), calc_sum_pow(0, 10));
	EXPECT_EQ((vector<field>{ 0, 1, 3, 6, 10, 15, 21, 28, 36, 45, 55 }), calc_sum_pow(1, 10));
	EXPECT_EQ((vector<field>{ 0, 1, 5, 14, 30, 55, 91, 140, 204, 285, 385 }), calc_sum_pow(2, 10));
	EXPECT_EQ((vector<field>{ 0, 1, 9, 36, 100, 225, 441, 784, 1296, 2025, 3025 }), calc_sum_pow(3, 10));
	EXPECT_EQ((vector<field>{ 0, 1, 129, 2316, 18700, 96825, 376761, 1200304, 3297456, 8080425, 18080425 }), calc_sum_pow(7, 10));
}

TEST(sums_test, sum_sqrt) {
	// Sum[[n/k], {k,1,n}]
	auto f0 = [](int m) { return m; };
	vector<int> ve0, va0;
	for (int n = 0; n < 100; n++) {
		ve0.push_back(sum<int>([&](int k) { return f0(n / k); }, 1, n));
		va0.push_back(sum_sqrt<int>(f0, n));
	}
	EXPECT_EQ(ve0, va0);

	// Sum[k*[n/k], {k,1,n}]
	auto f1 = [](int k, int m) { return k * m; };
	auto sf1 = [](int k, int m) { return sum_pow<int>(1, k) * m; };
	vector<int> ve1, va1;
	for (int n = 0; n < 100; n++) {
		ve1.push_back(sum<int>([&](int k) { return f1(k, n / k); }, 1, n));
		va1.push_back(sum_sqrt2<int>(sf1, n));
	}
	EXPECT_EQ(ve1, va1);

	// Sum[k^2[n/k]^2, {k,1,n}]
	auto f2 = [](int k, int m) { return k * k * m * m; };
	auto sf2 = [](int k, int m) { return sum_pow<int>(2, k) * m * m; };
	vector<int> ve2, va2;
	for (int n = 0; n < 100; n++) {
		ve2.push_back(sum<int>([&](int k) { return f2(k, n / k); }, 1, n));
		va2.push_back(sum_sqrt2<int>(sf2, n));
	}
	EXPECT_EQ(ve2, va2);
	
	// Sum[(k+3)*[n/k+2]^2, {k,1,n}]
	auto f3 = [](int k) { return k + 3; };
	auto sf3 = [](int k) { return sum_pow<int>(1, k) + 3 * k; };
	auto g3 = [](int m) { return m + 2; };
	vector<int> ve3, va3a, va3b;
	for (int n = 0; n < 100; n++) {
		ve3.push_back(sum<int>([&](int k) { return f3(k) * g3(n / k); }, 1, n));
		va3a.push_back(sum_sqrt2m<int>(sf3, g3, n));
		va3b.push_back(sum_sqrt2m<int>(f3, sf3, g3, n, 0));
	}
	EXPECT_EQ(ve3, va3a);
	EXPECT_EQ(ve3, va3b);
}

vector<field> calc_sum_m(int e, const poly& g, int n) {
	// initialize polynomials
	poly p = powT(poly{ 0, 1 }, e);
	poly sp = polynom_sum(p);
	poly st = polynom_sum(g * p);
	// wrapping functions that evaluate polynomials
	auto _g = [&](int n){ return g(field(n)); };
	auto _p = [&](int n){ return p(field(n)); };
	auto _sp = [&](int n){ return sp(field(n)); };
	auto _st = [&](int n){ return st(field(n)); };
	
	// preprocess `U = n^(2/3)` values of `Sum[p(k) * f[k], {k, 1, U}]`
	int U = (int)isq(icbrt(n));
	unordered_map<int, field> msf(U);
	std::vector<int> mu(U); moebius_mu(&mu[0], U);
	moebius_transform(msf, U, _g, &mu[0]);
	for (int k = 1; k < U; k++) msf[k] = _p(k) * msf[k] + msf[k - 1];

	// calculate all values
	vector<field> v;
	for (int k = 0; k <= n; k++) {
		v.push_back(sum_m<field>(k, _st, _sp, msf));
	}
	return v;

}
TEST(sums_test, sum_s) {
	// Sum[k^p euler_phi(k), { k, 1, n }]
	// euler_phi(x) is Euler Totient function defined as:
	//   Sum[[GCD(x, y)==1], {y,1,x}]
	//   Sum[mu(n/d) * d, {d|x}]
	poly g_phi = poly{ 0, 1 }; // d
	EXPECT_EQ((vector<field>{0, 1, 2, 4, 6, 10, 12, 18, 22, 28, 32, 42, 46, 58, 64, 72, 80, 96, 102, 120, 128}), calc_sum_m(0, g_phi, 20));
	EXPECT_EQ((vector<field>{0, 1, 3, 9, 17, 37, 49, 91, 123, 177, 217, 327, 375, 531, 615, 735, 863, 1135, 1243, 1585, 1745}), calc_sum_m(1, g_phi, 20));
	EXPECT_EQ((vector<field>{0, 1, 5, 23, 55, 155, 227, 521, 777, 1263, 1663, 2873, 3449, 5477, 6653, 8453, 10501, 15125, 17069, 23567, 26767}), calc_sum_m(2, g_phi, 20));
	
	// Sum[k^p euler_phi2(k), { k, 1, n }]
	// euler_phi2(x) is Euler Totient function in 2D defined as:
	//   Sum[[GCD(x, y, z)==1], {y,1,x}, {z,1,y}]
	//   Sum[mu(n/d) * d * (d + 1) / 2, {d|x}]
	poly g_phi2 = poly{ 0, 1, 1 } / field(2); // d * (d + 1) / 2
	EXPECT_EQ((vector<field>{0, 1, 3, 8, 15, 29, 42, 69, 95, 134, 172, 237, 287, 377, 452, 552, 652, 804, 915, 1104, 1252}), calc_sum_m(0, g_phi2, 20));
	EXPECT_EQ((vector<field>{0, 1, 5, 20, 48, 118, 196, 385, 593, 944, 1324, 2039, 2639, 3809, 4859, 6359, 7959, 10543, 12541, 16132, 19092}), calc_sum_m(1, g_phi2, 20));
	EXPECT_EQ((vector<field>{0, 1, 9, 54, 166, 516, 984, 2307, 3971, 7130, 10930, 18795, 25995, 41205, 55905, 78405, 104005, 147933, 183897, 252126, 311326}), calc_sum_m(2, g_phi2, 20));
}

TEST(sums_test, mertens) {
	int n = 30;
	// preprocess `U = n^(2/3)` values of `Sum[p(k) * f[k], {k, 1, U}]`
	int U = (int)isq(icbrt(n));
	unordered_map<int, field> mm(U);
	std::vector<int> mu(U); moebius_mu(&mu[0], U);
	for (int k = 1; k < U; k++) mm[k] = mm[k - 1] + mu[k];
	
	vector<field> va;
	for (int k = 0; k <= n; k++) {
		va.push_back(mertens(k, mm, field(1)));
	}
	EXPECT_EQ((vector<field>{0, 1, 0, -1, -1, -2, -1, -2, -2, -2, -1, -2, -2, -3, -2, -1, -1, -2, -2, -3, -3, -2, -1, -2, -2, -2, -1, -1, -1, -2, -3}), va);
}

TEST(sums_test, sum_primes) {
	vector<int> vp(isqrt(1030) + 1);
	int m = primes(&vp[0], nullptr, (int)vp.size());
	vector<field> va1, va2;
	for (int n = 0; n < 30; n++) {
		va1.push_back(sum_primes(n, &vp[0], field(1)));
		va2.push_back(sum_primes(1000 + n, &vp[0], field(1)) - 76127);
	}
	EXPECT_EQ((vector<field>{0, 0, 2, 5, 5, 10, 10, 17, 17, 17, 17, 28, 28, 41, 41, 41, 41, 58, 58, 77, 77, 77, 77, 100, 100, 100, 100, 100, 100, 129}), va1);
	EXPECT_EQ((vector<field>{0, 0, 0, 0, 0, 0, 0, 0, 0, 1009, 1009, 1009, 1009, 2022, 2022, 2022, 2022, 2022, 2022, 3041, 3041, 4062, 4062, 4062, 4062, 4062, 4062, 4062, 4062, 4062}), va2);
}

TEST(sums_test, sum_primes2) {
	vector<int> vp(1000);
	vector<char> vq(1000);
	int m = primes(&vp[0], &vq[0], (int)vq.size());
	vector<int> ve, va;
	int c = 0;
	for (int n = 0; n < vp.size(); n++) {
		ve.push_back(c += n * vq[n]);
		va.push_back(sum_primes(n, &vp[0], 1));
	}
	EXPECT_EQ(ve, va);
}
