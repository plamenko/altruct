#include "algorithm/math/combinatorics.h"

#include <algorithm>
#include <vector>

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

TEST(combinatorics_test, next_partition) {
	vector<int> v{ 5, 0, 0, 0, 0 };
	EXPECT_EQ(true, next_partition(v.begin(), v.end()));
	EXPECT_EQ((vector<int>{4, 1, 0, 0, 0}), v);
	EXPECT_EQ(true, next_partition(v.begin(), v.end()));
	EXPECT_EQ((vector<int>{3, 2, 0, 0, 0}), v);
	EXPECT_EQ(true, next_partition(v.begin(), v.end()));
	EXPECT_EQ((vector<int>{3, 1, 1, 0, 0}), v);
	EXPECT_EQ(true, next_partition(v.begin(), v.end()));
	EXPECT_EQ((vector<int>{2, 2, 1, 0, 0}), v);
	EXPECT_EQ(true, next_partition(v.begin(), v.end()));
	EXPECT_EQ((vector<int>{2, 1, 1, 1, 0}), v);
	EXPECT_EQ(true, next_partition(v.begin(), v.end()));
	EXPECT_EQ((vector<int>{1, 1, 1, 1, 1}), v);
	EXPECT_EQ(false, next_partition(v.begin(), v.end()));
	EXPECT_EQ((vector<int>{5, 0, 0, 0, 0}), v);

	vector<int> v2{ 5, 0, 0 };
	EXPECT_EQ(true, next_partition(v2.begin(), v2.end()));
	EXPECT_EQ((vector<int>{4, 1, 0}), v2);
	EXPECT_EQ(true, next_partition(v2.begin(), v2.end()));
	EXPECT_EQ((vector<int>{3, 2, 0}), v2);
	EXPECT_EQ(true, next_partition(v2.begin(), v2.end()));
	EXPECT_EQ((vector<int>{3, 1, 1}), v2);
	EXPECT_EQ(true, next_partition(v2.begin(), v2.end()));
	EXPECT_EQ((vector<int>{2, 2, 1}), v2);
	EXPECT_EQ(false, next_partition(v2.begin(), v2.end()));
	EXPECT_EQ((vector<int>{5, 0, 0}), v2);
}

TEST(combinatorics_test, next_combination) {
	vector<int> v{ 1, 2, 3, 4, 5 };
	EXPECT_EQ(true, next_combination(v.begin(), v.begin() + 3, v.end()));
	EXPECT_EQ((vector<int>{1, 2, 4, 3, 5}), v);
	EXPECT_EQ(true, next_combination(v.begin(), v.begin() + 3, v.end()));
	EXPECT_EQ((vector<int>{1, 2, 5, 3, 4}), v);
	EXPECT_EQ(true, next_combination(v.begin(), v.begin() + 3, v.end()));
	EXPECT_EQ((vector<int>{1, 3, 4, 2, 5}), v);
	EXPECT_EQ(true, next_combination(v.begin(), v.begin() + 3, v.end()));
	EXPECT_EQ((vector<int>{1, 3, 5, 2, 4}), v);
	EXPECT_EQ(true, next_combination(v.begin(), v.begin() + 3, v.end()));
	EXPECT_EQ((vector<int>{1, 4, 5, 2, 3}), v);
	EXPECT_EQ(true, next_combination(v.begin(), v.begin() + 3, v.end()));
	EXPECT_EQ((vector<int>{2, 3, 4, 1, 5}), v);
	EXPECT_EQ(true, next_combination(v.begin(), v.begin() + 3, v.end()));
	EXPECT_EQ((vector<int>{2, 3, 5, 1, 4}), v);
	EXPECT_EQ(true, next_combination(v.begin(), v.begin() + 3, v.end()));
	EXPECT_EQ((vector<int>{2, 4, 5, 1, 3}), v);
	EXPECT_EQ(true, next_combination(v.begin(), v.begin() + 3, v.end()));
	EXPECT_EQ((vector<int>{3, 4, 5, 1, 2}), v);
	EXPECT_EQ(false, next_combination(v.begin(), v.begin() + 3, v.end()));
	EXPECT_EQ((vector<int>{1, 2, 3, 4, 5}), v);

	EXPECT_EQ(true, next_combination(v.begin(), v.begin() + 1, v.end()));
	EXPECT_EQ((vector<int>{2, 1, 3, 4, 5}), v);
	EXPECT_EQ(true, next_combination(v.begin(), v.begin() + 1, v.end()));
	EXPECT_EQ((vector<int>{3, 1, 2, 4, 5}), v);
	EXPECT_EQ(true, next_combination(v.begin(), v.begin() + 1, v.end()));
	EXPECT_EQ((vector<int>{4, 1, 2, 3, 5}), v);
	EXPECT_EQ(true, next_combination(v.begin(), v.begin() + 1, v.end()));
	EXPECT_EQ((vector<int>{5, 1, 2, 3, 4}), v);
	EXPECT_EQ(false, next_combination(v.begin(), v.begin() + 1, v.end()));
	EXPECT_EQ((vector<int>{1, 2, 3, 4, 5}), v);

	EXPECT_EQ(false, next_combination(v.begin(), v.begin(), v.end()));
	EXPECT_EQ((vector<int>{1, 2, 3, 4, 5}), v);
	
	EXPECT_EQ(false, next_combination(v.begin(), v.end(), v.end()));
	EXPECT_EQ((vector<int>{1, 2, 3, 4, 5}), v);
}

TEST(combinatorics_test, nth_permutation) {
	char p0[] = "abcd";
	int l0 = (int)strlen(p0);
	for (int l = l0; l >= 0; l--) {
		p0[l] = 0;
		string s = p0;
		for (int k = 0; k < 100; k++) {
			string p = p0;
			nth_permutation(p.begin(), p.end(), k);
			EXPECT_EQ(s, p);
			next_permutation(s.begin(), s.end());
		}
	}
}
