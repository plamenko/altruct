#include "algorithm/collections/collections.h"

#include "gtest/gtest.h"

#include <string>
#include <set>
#include <vector>

using namespace std;
using namespace altruct::collections;

TEST(collections_test, sorted) {
    EXPECT_EQ((vector<int>{}), sorted(vector<int>()));
    EXPECT_EQ((vector<int>{5, 7, 13}), sorted(vector<int>{7, 13, 5}));
    EXPECT_EQ((vector<int>{5, 7, 13}), sorted(set<int>{7, 13, 5}));
}

TEST(collections_test, reversed) {
    EXPECT_EQ((vector<int>{}), reversed(vector<int>()));
    EXPECT_EQ((vector<int>{5, 13, 7}), reversed(vector<int>{7, 13, 5}));
    EXPECT_EQ((vector<int>{13, 7, 5}), reversed(set<int>{7, 13, 5}));
}

TEST(collections_test, take) {
    EXPECT_EQ((vector<int>{}), take(vector<int>(), 3));
    EXPECT_EQ((vector<int>{2, 3, 5}), take(vector<int>{2, 3, 5, 7, 11, 13}, 3));
    EXPECT_EQ((vector<int>{2, 3, 5, 7, 11, 13}), take(set<int>{2, 3, 5, 7, 11, 13}, 100));
}

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

TEST(collections_test, reserve_more) {
	vector<int> v;
	// set initial capacity to 100
	v.reserve(100);
	EXPECT_EQ(0, v.size());
	EXPECT_EQ(100, v.capacity());
	
	// no-op, vector already has enough capacity to accomodate 10 new elements
	reserve_more(v, 10);
	EXPECT_EQ(0, v.size());
	EXPECT_EQ(100, v.capacity());
	
	// vector doesn't have enough capacity to accomodate 110 new elements, perform exponential growth
	reserve_more(v, 110);
	EXPECT_EQ(0, v.size());
	EXPECT_EQ(150, v.capacity());

	// vector doesn't have enough capacity to accomodate 1000 new elements, reserve necessary space
	reserve_more(v, 1000);
	EXPECT_EQ(0, v.size());
	EXPECT_EQ(1000, v.capacity());

	// no reallocation, vector already has enough capacity to accomodate 900 new elements
	for (int i = 0; i < 900; i++) v.push_back(i);
	EXPECT_EQ(900, v.size());
	EXPECT_EQ(1000, v.capacity());

	// vector doesn't have enough capacity to accomodate 500 new elements, perform exponential growth
	reserve_more(v, 500);
	EXPECT_EQ(900, v.size());
	EXPECT_EQ(1500, v.capacity());

	// vector doesn't have enough capacity to accomodate 5000 new elements, reserve necessary space
	reserve_more(v, 5000);
	EXPECT_EQ(900, v.size());
	EXPECT_EQ(5900, v.capacity());

	// no reallocation, vector already has enough capacity to accomodate 5000 new elements
	for (int i = 0; i < 5000; i++) v.push_back(i);
	EXPECT_EQ(5900, v.size());
	EXPECT_EQ(5900, v.capacity());
}
