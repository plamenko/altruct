#include "algorithm/graph/iterative_dfs.h"
#include "algorithm/graph/topological_sort.h"
#include "algorithm/graph/dinic_flow.h"
#include "algorithm/graph/push_relabel_flow.h"
#include "algorithm/graph/bipartite_matching.h"
#include "algorithm/graph/lowest_common_ancestor.h"
#include "algorithm/graph/heavy_light_decomposition.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::graph;

namespace {
    vector<vector<pair<int, int>>> dag_adjl{
        { { 4, 10 } },
        { { 4, 7 }, { 0, 3 } },
        { { 0, 5 } },
        { { 2, -4 }, { 0, 6 }, { 1, 8 }, { 4, 5 } },
        {},
        { { 1, -2 }, { 6, 6 } },
        {},
        { { 6, 7 } },
        { { 10, -5 } },
        {},
        {},
    };
    auto index_f = [](const pair<int, int>& e){return e.first; };
}

TEST(graph_test, iterative_dfs) {
    // TODO
}

TEST(graph_test, in_degrees) {
    EXPECT_EQ((vector<int>{3, 2, 1, 0, 3, 0, 2, 0, 0, 0, 1}), in_degrees(dag_adjl, index_f));
}

TEST(graph_test, topological_sort) {
    EXPECT_EQ((vector<int>{ 9, 8, 10, 7, 5, 6, 3, 1, 2, 0, 4 }), topological_sort(dag_adjl, index_f));
}

template<typename MAX_FLOW_IMPL, typename T>
void test_max_flow(const vector<vector<T>>& capacities, const vector<vector<T>>& expected_flows) {
    MAX_FLOW_IMPL mfi(capacities, 1000000);
    vector<vector<T>> actual_flows = capacities;
    for (int i = 0; i < mfi.cap.size(); i++) {
        for (int j = 0; j < mfi.cap.size(); j++) {
            actual_flows[i][j] = mfi.calc_max_flow(i, j);
        }
    }
    EXPECT_EQ(expected_flows, actual_flows);
}

TEST(graph_test, dinic_flow) {
    test_max_flow<dinic_flow<int>, int>({ { 0 } }, { { 0 } });
    test_max_flow<dinic_flow<int>, int>({ { 0, 5 }, { 7, 0 } }, { { 0, 5 }, { 7, 0 } });
    test_max_flow<dinic_flow<int>, int>({ { 0, 3, 5 }, { 0, 0, 2 }, { 0, 0, 0 } }, { { 0, 3, 7 }, { 0, 0, 2 }, { 0, 0, 0 } });
    test_max_flow<dinic_flow<double>, double>({ { 0, 5, 2 }, { 7, 0, 4 }, { 1, 3, 0 } }, { { 0, 7, 6 }, { 8, 0, 6 }, { 4, 4, 0 } });
}

TEST(graph_test, push_relabel_flow) {
    test_max_flow<push_relabel_flow<int>, int>({ { 0 } }, { { 0 } });
    test_max_flow<push_relabel_flow<int>, int>({ { 0, 5 }, { 7, 0 } }, { { 0, 5 }, { 7, 0 } });
    test_max_flow<push_relabel_flow<int>, int>({ { 0, 3, 5 }, { 0, 0, 2 }, { 0, 0, 0 } }, { { 0, 3, 7 }, { 0, 0, 2 }, { 0, 0, 0 } });
    test_max_flow<push_relabel_flow<double>, double>({ { 0, 5, 2 }, { 7, 0, 4 }, { 1, 3, 0 } }, { { 0, 7, 6 }, { 8, 0, 6 }, { 4, 4, 0 } });
}

TEST(graph_test, bipartite_matching) {
    typedef std::vector<std::pair<int, int>> edges_t;
    EXPECT_EQ((edges_t()), bipartite_matching(0, edges_t()));
    EXPECT_EQ((edges_t{ { 0, 2 }, { 1, 3 } }), bipartite_matching(4, edges_t{ { 0, 2 }, { 0, 3 }, { 1, 3 } }));
    EXPECT_EQ((edges_t{ { 0, 2 }, { 1, 3 } }), bipartite_matching(4, edges_t{ { 0, 2 }, { 1, 2 }, { 1, 3 } }));
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
