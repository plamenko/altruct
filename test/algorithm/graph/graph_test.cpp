#include "algorithm/graph/dinic.h"
#include "algorithm/graph/lowest_common_ancestor.h"
#include "algorithm/graph/heavy_light_decomposition.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::graph;

TEST(graph_test, dinic) {
	vector<vector<int>> cap{ { 0, 3, 5 }, { 0, 0, 2 }, { 0, 0, 0 } };
	dinic<int> d(cap);
	EXPECT_EQ(7, d.calc_max_flow(0, 2));
}

TEST(graph_test, lowest_common_ancestor) {
	vector<vector<int>> adjl{ { 1, 2 }, { 0 }, { 0, 3 }, { 2 } };
	lowest_common_ancestor lca(adjl);
	EXPECT_EQ(0, lca.ancestor(1, 3));
}

TEST(graph_test, heavy_light_decomposition) {
	vector<vector<int>> adjl{ { 1, 2 }, { 0 }, { 0, 3 }, { 2 } };
	heavy_light_decomposition_ex hld(adjl);
	EXPECT_EQ(3, hld.parent(3, 0));
	EXPECT_EQ(2, hld.parent(3, 1));
	EXPECT_EQ(0, hld.parent(3, 2));
}
