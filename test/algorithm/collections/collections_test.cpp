#include "algorithm/collections/collections.h"

#include "gtest/gtest.h"

#include <string>
#include <set>
#include <vector>

using namespace std;
using namespace altruct::collections;

TEST(collections_test, filter) {
	EXPECT_EQ((vector<int>{}), filter(vector<int>(), [](int x){ return x % 4 == 1; }));
	EXPECT_EQ((vector<int>{5, 13}), filter(vector<int>{2, 3, 5, 7, 11, 13}, [](int x){ return x % 4 == 1; }));
	EXPECT_EQ((vector<int>{5, 13}), filter(set<int>{2, 3, 5, 7, 11, 13}, [](int x){ return x % 4 == 1; }));
}

TEST(collections_test, transform) {
	EXPECT_EQ((vector<int>{}), transform(vector<int>(), [](int x){return x*x; }));
	EXPECT_EQ((vector<int>{4, 9, 25, 49}), transform(vector<int>{-2, 3, 5, 7}, [](int x){return x*x; }));
	EXPECT_EQ((vector<int>{4, 9, 25, 49}), transform(set<int>{-2, 3, 5, 7}, [](int x){return x*x; }));
}

TEST(collections_test, run_length) {
	EXPECT_EQ((vector<pair<string, int>>{}), run_length(vector<string>()));
	EXPECT_EQ((vector<pair<string, int>>{{ "a", 3 }, { "b", 1 }, { "c", 2 }, { "a", 1 }}), run_length(vector<string>{"a", "a", "a", "b", "c", "c", "a"}));
}

TEST(collections_test, compare) {
	string s1 = "banana";
	string s2 = "bambus";
	string s3 = "bambu bambu";
	string s4 = "bananana";
	EXPECT_EQ(+1, compare(s1.begin(), s1.end(), s2.begin(), s2.end()));
	EXPECT_EQ(-1, compare(s2.begin(), s2.end(), s1.begin(), s1.end()));
	EXPECT_EQ(+1, compare(s1.begin(), s1.end(), s3.begin(), s3.end()));
	EXPECT_EQ(-1, compare(s3.begin(), s3.end(), s1.begin(), s1.end()));
	
	EXPECT_EQ(+1, compare(s2.begin(), s2.end(), s3.begin(), s3.end()));
	EXPECT_EQ(-1, compare(s3.begin(), s3.end(), s2.begin(), s2.end()));
	EXPECT_EQ(+1, compare(s2.begin(), s2.end(), s3.begin(), s3.end(), 6));
	EXPECT_EQ(-1, compare(s3.begin(), s3.end(), s2.begin(), s2.end(), 6));
	EXPECT_EQ(0, compare(s2.begin(), s2.end(), s3.begin(), s3.end(), 5));
	EXPECT_EQ(0, compare(s3.begin(), s3.end(), s2.begin(), s2.end(), 5));
	EXPECT_EQ(0, compare(s2.begin(), s2.end(), s3.begin(), s3.end(), 1));
	EXPECT_EQ(0, compare(s3.begin(), s3.end(), s2.begin(), s2.end(), 1));

	EXPECT_EQ(0, compare(s3.begin(), s3.end(), s3.begin(), s3.end()));
	EXPECT_EQ(0, compare(s3.begin(), s3.end(), s3.begin(), s3.end(), 1000));
	EXPECT_EQ(0, compare(s3.begin(), s3.end(), s3.begin(), s3.end(), 5));
	EXPECT_EQ(0, compare(s3.begin(), s3.end(), s3.begin(), s3.end(), 0));
	
	EXPECT_EQ(-1, compare(s1.begin(), s1.end(), s4.begin(), s4.end()));
	EXPECT_EQ(-1, compare(s1.begin(), s1.end(), s4.begin(), s4.end(), 1000));
	EXPECT_EQ(-1, compare(s1.begin(), s1.end(), s4.begin(), s4.end(), 7));
	EXPECT_EQ(0, compare(s1.begin(), s1.end(), s4.begin(), s4.end(), 6));
	EXPECT_EQ(0, compare(s1.begin(), s1.end(), s4.begin(), s4.end(), 0));
}
