#include "algorithm/collections/collections.h"

#include "gtest/gtest.h"

#include <string>
#include <set>
#include <vector>

using namespace std;
using namespace altruct::collections;

TEST(reduce_test, transform) {
	EXPECT_EQ((vector<int>{}), transform(vector<int>(), [](int x){return x*x; }));
	EXPECT_EQ((vector<int>{4, 9, 25, 49}), transform(vector<int>{-2, 3, 5, 7}, [](int x){return x*x; }));
	EXPECT_EQ((vector<int>{4, 9, 25, 49}), transform(set<int>{-2, 3, 5, 7}, [](int x){return x*x; }));
}

TEST(reduce_test, run_length) {
	EXPECT_EQ((vector<pair<string, int>>{}), run_length(vector<string>()));
	EXPECT_EQ((vector<pair<string, int>>{{ "a", 3 }, { "b", 1 }, { "c", 2 }, { "a", 1 }}), run_length(vector<string>{"a", "a", "a", "b", "c", "c", "a"}));
}
