#include "structure/graph/disjoint_set.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::graph;

vector<size_t> get_all(disjoint_set& ds) {
	vector<size_t> v;
	for (int i = 0; i < ds.size(); i++) {
		v.push_back(ds.find(i));
	}
	return v;
}

TEST(disjoint_set_test, constructor) {
	disjoint_set ds0;
	EXPECT_EQ(0, ds0.size());
	EXPECT_EQ((vector<size_t>{}), get_all(ds0));
	disjoint_set ds1(10);
	EXPECT_EQ(10, ds1.size());
	EXPECT_EQ((vector<size_t>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}), get_all(ds1));
}

TEST(disjoint_set_test, clear) {
	disjoint_set ds;
	ds.unite(1, 2);
	ds.unite(2, 3);
	ds.unite(4, 5);
	EXPECT_NE((vector<size_t>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}), get_all(ds));

	ds.clear(10);
	EXPECT_EQ(10, ds.size());
	EXPECT_EQ((vector<size_t>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}), get_all(ds));
	
	ds.clear();
	EXPECT_EQ(0, ds.size());
	EXPECT_EQ((vector<size_t>{}), get_all(ds));
}

TEST(disjoint_set_test, unite) {
	disjoint_set ds;
	EXPECT_TRUE(ds.unite(1, 2));
	EXPECT_TRUE(ds.unite(2, 3));
	EXPECT_FALSE(ds.unite(1, 3)); // already same
	EXPECT_TRUE(ds.unite(6, 7));
	EXPECT_EQ((vector<size_t>{0, 1, 1, 1, 4, 5, 6, 6}), get_all(ds));
}

TEST(disjoint_set_test, find) {
	disjoint_set ds;
	EXPECT_EQ(1, ds.find(1));
	EXPECT_EQ(2, ds.find(2));
	EXPECT_EQ(3, ds.find(3));
	EXPECT_EQ(4, ds.find(4));
	ds.unite(1, 2);
	EXPECT_EQ(1, ds.find(1));
	EXPECT_EQ(1, ds.find(2));
	EXPECT_EQ(3, ds.find(3));
	EXPECT_EQ(4, ds.find(4));
	ds.unite(4, 2);
	EXPECT_EQ(1, ds.find(1));
	EXPECT_EQ(1, ds.find(2));
	EXPECT_EQ(3, ds.find(3));
	EXPECT_EQ(1, ds.find(4));
	EXPECT_EQ((vector<size_t>{0, 1, 1, 3, 1}), get_all(ds));
}

TEST(disjoint_set_test, count) {
	disjoint_set ds;
	EXPECT_EQ(1, ds.count(1));
	EXPECT_EQ(1, ds.count(2));
	EXPECT_EQ(1, ds.count(3));
	EXPECT_EQ(1, ds.count(4));
	ds.unite(1, 2);
	EXPECT_EQ(2, ds.count(1));
	EXPECT_EQ(2, ds.count(2));
	EXPECT_EQ(1, ds.count(3));
	EXPECT_EQ(1, ds.count(4));
	ds.unite(4, 2);
	EXPECT_EQ(3, ds.count(1));
	EXPECT_EQ(3, ds.count(2));
	EXPECT_EQ(1, ds.count(3));
	EXPECT_EQ(3, ds.count(4));
	EXPECT_EQ((vector<size_t>{0, 1, 1, 3, 1}), get_all(ds));
}
