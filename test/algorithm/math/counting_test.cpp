#include "algorithm/collections/collections.h"
#include "algorithm/math/counting.h"
#include "structure/math/modulo.h"

#include <algorithm>
#include <vector>

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;
using namespace altruct::collections;

vector<vector<moduloX<int>>> to_modx(int M, const vector<vector<int>>& vv) {
	return transform(vv, [&](const vector<int>& v) {
		return transform(v, [&](int a) {
			return moduloX<int>(a, M);
		});
	});
}

template<typename T>
vector<T> get_column(const vector<vector<T>>& vv, int k) {
	return transform(vv, [&](const vector<T>& v){
		return (k < v.size()) ? v[k] : zeroT<T>::of(v[0]);
	});
}

TEST(counting_test, stirling_s1) {
	typedef moduloX<int> modx;
	modx id(1, 1000000007);
	vector<vector<modx>> expected_all = to_modx(id.M(), {
		{ 1 },
		{ 0, 1 },
		{ 0, -1, 1 },
		{ 0, 2, -3, 1 },
		{ 0, -6, 11, -6, 1 },
		{ 0, 24, -50, 35, -10, 1 },
		{ 0, -120, 274, -225, 85, -15, 1 },
		{ 0, 720, -1764, 1624, -735, 175, -21, 1 },
		{ 0, -5040, 13068, -13132, 6769, -1960, 322, -28, 1 },
		{ 0, 40320, -109584, 118124, -67284, 22449, -4536, 546, -36, 1 },
		{ 0, -362880, 1026576, -1172700, 723680, -269325, 63273, -9450, 870, -45, 1 },
	});

	auto actual_all = stirling_s1_all(11, 11, id);
	EXPECT_EQ(expected_all, actual_all);

	auto expected_4 = transform(expected_all, [](vector<modx> v) { v.resize(min((int)v.size(), 4)); return v; });
	auto actual_4 = stirling_s1_all(11, 4, id);
	EXPECT_EQ(expected_4, actual_4);

	for (int k = 0; k < (int)expected_all.size(); k++) {
		auto expected_n = get_column(expected_all, k);
		auto actual_n = stirling_s1_all_n_for_k((int)expected_n.size(), k, id);
		EXPECT_EQ(expected_n, actual_n) << "k = " << k;
	}

	for (int n = 0; n < (int)expected_all.size(); n++) {
		auto actual_k = stirling_s1_all_k_for_n(n, id);
		EXPECT_EQ(expected_all[n], actual_k) << "n = " << n;
	}

	for (int n = 0; n < (int)expected_all.size(); n++) {
		for (int k = 0; k <= n; k++) {
			auto actual_n_k = stirling_s1(n, k, id);
			EXPECT_EQ(expected_all[n][k], actual_n_k) << "n = " << n << ", k = " << k;
		}
	}
}

TEST(counting_test, stirling_s2) {
	typedef moduloX<int> modx;
	modx id(1, 1000000007);
	vector<vector<modx>> expected_all = to_modx(id.M(), {
		{ 1 },
		{ 0, 1 },
		{ 0, 1, 1 },
		{ 0, 1, 3, 1 },
		{ 0, 1, 7, 6, 1 },
		{ 0, 1, 15, 25, 10, 1 },
		{ 0, 1, 31, 90, 65, 15, 1 },
		{ 0, 1, 63, 301, 350, 140, 21, 1 },
		{ 0, 1, 127, 966, 1701, 1050, 266, 28, 1 },
		{ 0, 1, 255, 3025, 7770, 6951, 2646, 462, 36, 1 },
		{ 0, 1, 511, 9330, 34105, 42525, 22827, 5880, 750, 45, 1 },
	});

	auto actual_all = stirling_s2_all(11, 11, id);
	EXPECT_EQ(expected_all, actual_all);

	auto expected_4 = transform(expected_all, [](vector<modx> v) { v.resize(min((int)v.size(), 4)); return v; });
	auto actual_4 = stirling_s2_all(11, 4, id);
	EXPECT_EQ(expected_4, actual_4);

	for (int k = 0; k < (int)expected_all.size(); k++) {
		auto expected_n = get_column(expected_all, k);
		auto actual_n = stirling_s2_all_n_for_k((int)expected_n.size(), k, id);
		EXPECT_EQ(expected_n, actual_n) << "k = " << k;
	}

	for (int n = 0; n < (int)expected_all.size(); n++) {
		auto actual_k = stirling_s2_all_k_for_n(n, id);
		EXPECT_EQ(expected_all[n], actual_k) << "n = " << n;
	}

	for (int n = 0; n < (int)expected_all.size(); n++) {
		for (int k = 0; k <= n; k++) {
			auto actual_n_k = stirling_s2(n, k, id);
			EXPECT_EQ(expected_all[n][k], actual_n_k) << "n = " << n << ", k = " << k;
		}
	}
}

TEST(counting_test, partitions_p) {
	typedef moduloX<int> modx;
	modx id(1, 1000000007);
	vector<modx> expected = to_modx(1009, { { 1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 101, 135, 176, 231, 297, 385, 490, 627 } })[0];
	auto actual = partitions_p(21, id);
	EXPECT_EQ(expected, actual);
}
