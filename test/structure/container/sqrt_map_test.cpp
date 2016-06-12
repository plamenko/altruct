#include "algorithm/math/base.h"
#include "structure/container/sqrt_map.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;
using namespace altruct::container;

TEST(sqrt_map_test, sqrt_map) {
	int n = 100;
	int q = isqrt(n);
	sqrt_map<int, double> m(q, n);
	for (int i = 1; i <= q; i++) {
		EXPECT_EQ(0, m.count(i)) << "unexpected element at " << i;
		m[i] = 1.0 / i;
	}
	for (int i = 1; i <= q; i++) {
		int k = n / i;
		EXPECT_EQ(k <= q ? 1 : 0, m.count(k)) << "unexpected element at " << k;
		m[k] = 1.0 / k;
	}

	for (int i = 1; i <= q; i++) {
		EXPECT_EQ(1, m.count(i)) << "no element at " << i;
		EXPECT_EQ(1.0 / i, m.at(i)) << "unexpected element at " << i;
	}
	for (int i = 1; i <= q; i++) {
		int k = n / i;
		EXPECT_EQ(1, m.count(k)) << "no element at " << k;
		EXPECT_EQ(1.0 / k, m.at(k)) << "unexpected element at " << k;
	}

	n = 200;
	m.reset_max(n);
	for (int i = 1; i <= q; i++) {
		EXPECT_EQ(1, m.count(i)) << "no element at " << i;
		EXPECT_EQ(1.0 / i, m.at(i)) << "unexpected element at " << i;
	}
	for (int i = 1; i <= q; i++) {
		int k = n / i;
		EXPECT_EQ(k <= q ? 1 : 0, m.count(k)) << "unexpected element at " << k;
		m[k] = 1.0 / k;
	}

	n = 50;
	m.reset_max(n);
	for (int i = 1; i <= q; i++) {
		EXPECT_EQ(1, m.count(i)) << "no element at " << i;
		EXPECT_EQ(1.0 / i, m.at(i)) << "unexpected element at " << i;
	}
	for (int i = 1; i <= q; i++) {
		int k = n / i;
		EXPECT_EQ(k <= q ? 1 : 0, m.count(k)) << "unexpected element at " << k;
		m[k] = 1.0 / k;
	}
}

void test_insert_erase(sqrt_map<int, double>& m, int key, double val, double val2) {
	EXPECT_EQ(0, m.count(key)) << "unexpected element at " << key;
	EXPECT_EQ(make_pair(key, true), m.insert({ key, val })) << "insert failed at " << key;
	EXPECT_EQ(1, m.count(key)) << "no element at " << key;
	EXPECT_EQ(val, m.at(key)) << "unexpected element at " << key;
	EXPECT_EQ(make_pair(key, false), m.insert({ key, val2 })) << "no element at " << key;
	EXPECT_EQ(1, m.count(key)) << "no element at " << key;
	EXPECT_EQ(val, m.at(key)) << "unexpected element at " << key;
	EXPECT_EQ(make_pair(key, false), m.insert({ key, val2 })) << "no element at " << key;
	EXPECT_EQ(1, m.count(key)) << "no element at " << key;
	EXPECT_EQ(val, m.at(key)) << "unexpected element at " << key;

	EXPECT_EQ(1, m.erase(key)) << "erase failed at " << key;
	EXPECT_EQ(0, m.erase(key)) << "unexpected element at " << key;
	EXPECT_EQ(0, m.erase(key)) << "unexpected element at " << key;

	EXPECT_EQ(make_pair(key, true), m.insert({ key, val2 })) << "insert failed at " << key;
	EXPECT_EQ(1, m.count(key)) << "no element at " << key;
	EXPECT_EQ(val2, m.at(key)) << "unexpected element at " << key;
}

TEST(sqrt_map_test, insert_erase) {
	int n = 100;
	int q = isqrt(n);
	sqrt_map<int, double> m(q, n);
	test_insert_erase(m, 17, 3.14, 2.71);
	test_insert_erase(m, 5, 0.61, 1.61);
}
