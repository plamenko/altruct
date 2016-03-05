#include "algorithm/math/reduce.h"

#include "gtest/gtest.h"

#include <set>
#include <vector>

using namespace std;
using namespace altruct::math;

TEST(reduce_test, reduce_sum) {
	EXPECT_EQ(0, reduce_sum(vector<int>{}));
	EXPECT_EQ(3, reduce_sum(vector<int>{3}));
	EXPECT_EQ(9, reduce_sum(vector<int>{4, 1, 2, 2}));
	EXPECT_EQ(7, reduce_sum(set<int>{4, 1, 2}));
}

TEST(reduce_test, reduce_product) {
	EXPECT_EQ(1, reduce_product(vector<int>{}));
	EXPECT_EQ(3, reduce_product(vector<int>{3}));
	EXPECT_EQ(36, reduce_product(vector<int>{4, 1, 3, 3}));
	EXPECT_EQ(12, reduce_product(set<int>{4, 1, 3}));
	EXPECT_EQ(0, reduce_product(set<int>{4, 1, 0, 3}));
}

TEST(reduce_test, reduce_min) {
	int inf = std::numeric_limits<int>::max();
	EXPECT_EQ(+inf, reduce_min(vector<int>{}));
	EXPECT_EQ(3, reduce_min(vector<int>{3}));
	EXPECT_EQ(-3, reduce_min(vector<int>{4, 1, -3, 3}));
	EXPECT_EQ(1, reduce_min(set<int>{4, 1, 3}));
}

TEST(reduce_test, reduce_max) {
	int inf = std::numeric_limits<int>::max();
	EXPECT_EQ(-inf, reduce_max(vector<int>{}));
	EXPECT_EQ(3, reduce_max(vector<int>{3}));
	EXPECT_EQ(4, reduce_max(vector<int>{4, 1, -5, 3}));
	EXPECT_EQ(4, reduce_max(set<int>{4, 1, 3}));
}

TEST(reduce_test, reduce_mex) {
	EXPECT_EQ(0, reduce_mex(vector<int>{}));
	EXPECT_EQ(0, reduce_mex(vector<int>{1, 2, 3}));
	EXPECT_EQ(1, reduce_mex(vector<int>{0, 2, 4, 6}));
	EXPECT_EQ(2, reduce_mex(set<int>{0, 1, 4, 7, 12}));
}
