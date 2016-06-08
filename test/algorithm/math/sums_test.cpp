#include "algorithm/math/sums.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

typedef modulo<int, 1000000007> field;

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
