#include "altruct/structure/graph/graph.h"
#include "altruct/algorithm/graph/iterative_dfs.h"
#include "altruct/algorithm/graph/topological_sort.h"
#include "altruct/algorithm/graph/tarjan_scc.h"
#include "altruct/algorithm/graph/chain_decomposition.h"
#include "altruct/algorithm/graph/chromatic_polynomial.h"
#include "altruct/algorithm/graph/transitive_closure.h"
#include "altruct/algorithm/graph/floyd_warshall.h"
#include "altruct/algorithm/graph/dijkstra.h"
#include "altruct/algorithm/graph/prim_spanning_tree.h"
#include "altruct/algorithm/graph/dinic_flow.h"
#include "altruct/algorithm/graph/push_relabel_flow.h"
#include "altruct/algorithm/graph/bipartite_matching.h"
#include "altruct/algorithm/graph/lowest_common_ancestor.h"
#include "altruct/algorithm/graph/heavy_light_decomposition.h"
#include "altruct/algorithm/graph/sat2.h"

#include "altruct/algorithm/collections/collections.h"
#include "altruct/io/iostream_overloads.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;
using namespace altruct::graph;
using namespace altruct::collections;

namespace {
    typedef weighted_edge<int> iedge;
    graph<iedge> dag_neg1({ // no cycles, has negative weights, multiple components
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
    });
    graph<iedge> cyc_neg1({ // has cycles, has negative weights, but no negative cycles
        { { 2, -2 } },
        { { 0, 4 }, { 2, 3 } },
        { { 3, 2 } },
        { { 1, -1 }, { 4, -8 } },
        { { 5, 2 } },
        { { 6, 3 } },
        { { 4, 7 } },
        { { 5, 10 }, { 5, 6 }, { 5, 11 } },
    });
    graph<iedge> cyc_pos1({ // has cycles, only positive weights
        { { 2, 2 } },
        { { 0, 4 }, { 2, 3 }, { 4, 20 } },
        { { 3, 2 } },
        { { 1, 1 }, { 4, 8 } },
        { { 5, 2 } },
        { { 6, 3 } },
        { { 4, 7 } },
        { { 5, 10 }, { 5, 6 }, { 5, 11 } },
    });
    graph<iedge> cyc_undir({ // undirected, has cycles, multiple components
/*0*/   { { 5, 21 }, { 7, 27 }, { 8, 23 } },
        { { 2, 31 }, { 4, 33 } },
        { { 1, 31 }, { 4, 35 } },
        { { 6, 28 }, { 9, 26 } },
        { { 1, 33 }, { 2, 35 }, { 9, 34 } },
/*5*/   { { 0, 21 }, { 7, 26 }, { 9, 29 } },
        { { 3, 28 }, { 9, 28 } },
        { { 0, 27 }, { 5, 26 }, { 8, 31 }, { 9, 30 } },
        { { 0, 23 }, { 7, 31 } },
        { { 4, 34 }, { 5, 29 }, { 6, 28 }, { 7, 30 }, { 3, 26 } },
/*10*/  { { 11, 45 } },
        { { 10, 45 }, { 12, 38 } },
        { { 11, 38 } },
        { { 14, 55 }, { 15, 57 } },
        { { 13, 55 }, { 15, 54 } },
/*15*/  { { 13, 57 }, { 14, 54 }, { 16, 58 }, { 17, 50 } },
        { { 15, 58 }, { 17, 55 } },
        { { 15, 50 }, { 16, 55 } },
    });
    graph<edge> cyc_undir2({ // undirected, has cycles
        { { 1 }, { 3 }, { 4 } },
        { { 0 }, { 2 }, { 4 } },
        { { 1 }, { 3 } },
        { { 2 }, { 0 }, { 5 } },
        { { 0 }, { 1 }, { 5 } },
        { { 3 }, { 4 } },
    });
    graph<edge> cyc_undir3({ // undirected, has cycles
        { { 1 }, { 3 }, { 4 } },
        { { 0 }, { 2 }, { 7 }, { 3 }, { 6 }, { 4 } },
        { { 1 }, { 3 }, { 4 } },
        { { 2 }, { 0 }, { 4 }, { 1 } },
        { { 0 }, { 5 }, { 3 }, { 2 }, { 1 } },
        { { 4 } },
        { { 1 }, { 7 } },
        { { 6 }, { 1 } },
    });
    graph<edge> cyc_undir4({ // undirected, has cycles
        { { 1 }, { 3 } },
        { { 6 }, { 7 }, { 0 }, { 2 } },
        { { 1 }, { 3 }, { 5 } },
        { { 4 }, { 2 }, { 0 } },
        { { 5 }, { 3 } },
        { { 2 }, { 4 } },
        { { 7 }, { 1 } },
        { { 6 }, { 1 } },
    });
    graph<edge> cyc_undir5({ // undirected, has cycles
        { { 1 }, { 7 } },
        { { 0 }, { 2 }, { 3 }, { 4 } },
        { { 1 }, { 3 } },
        { { 1 }, { 2 } },
        { { 1 }, { 5 }, { 6 }, { 7 }, { 11 } },
        { { 4 }, { 6 } },
        { { 4 }, { 5 } },
        { { 0 }, { 4 }, { 8 }, { 9 }, { 10 } },
        { { 7 }, { 9 } },
        { { 7 }, { 8 } },
        { { 7 }, { 11 } },
        { { 4 }, { 10 } },
    });
    graph<edge> blocks1({ // biconnected blocks with satelite trees
        { { 1 }, { 3 } },
        { { 0 }, { 2 } },
        { { 1 }, { 3 }, { 6 }, { 8 } },
        { { 0 }, { 2 }, { 4 }, { 5 }, { 9 }, { 10 } },
        { { 3 }, { 5 }, { 13 } },
        { { 3 }, { 4 } },
        { { 2 }, { 7 } },
        { { 6 }, { 8 }, { 11 }, { 12 } },
        { { 2 }, { 7 } },
        { { 3 }, { 10 } },
        { { 3 }, { 9 }, { 14 } },
        { { 7 } },
        { { 7 } },
        { { 4 } },
        { { 10 }, { 15 }, { 16 } },
        { { 14 } },
        { { 14 } },
        { { 18 }, { 19 }, { 23 } },
        { { 17 }, { 19 }, { 20 } },
        { { 17 }, { 17 }, { 21 } },
        { { 18 } },
        { { 19 }, { 22 } },
        { { 21 } },
        { { 17 }, { 24 }, { 25 } },
        { { 23 }, { 25 } },
        { { 23 }, { 24 } },
    });
    graph<edge> blocks1_bctree({
        { { 7 }, { 6 } },
        { { 6 }, { 9 } },
        { { 7 }, { 8 } },
        { { 7 }, { 10 } },
        { { 12 }, { 14 }, { 13 } },
        { { 16 } },
        { { 0 }, { 1 } },
        { { 0 }, { 2 }, { 3 } },
        { { 2 }, { 19 } },
        { { 1 }, { 17 }, { 18 } },
        { { 3 }, { 11 } },
        { { 10 }, { 20 }, { 21 } },
        { { 4 }, { 16 } },
        { { 4 }, { 22 } },
        { { 4 }, { 15 } },
        { { 23 }, { 14 } },
        { { 5 }, { 12 } },
        { { 9 } },
        { { 9 } },
        { { 8 } },
        { { 11 } },
        { { 11 } },
        { { 13 } },
        { { 15 } },
    });
    vector<int> blocks1_map{ 0, 0, 6, 7, 8, 2, 1, 9, 1, 3, 10, 17, 18, 19, 11, 20, 21, 12, 13, 14, 22, 15, 23, 16, 5, 5 };
}

TEST(graph_algorithms_test, iterative_dfs) {
    // TODO
}

TEST(graph_algorithms_test, in_degrees) {
    EXPECT_EQ((vector<int>{3, 2, 1, 0, 3, 0, 2, 0, 0, 0, 1}), in_degrees(dag_neg1));
}

TEST(graph_algorithms_test, topological_sort) {
    EXPECT_EQ((vector<int>{ 9, 8, 10, 7, 5, 6, 3, 1, 2, 0, 4 }), topological_sort(dag_neg1));
}

TEST(graph_algorithms_test, tarjan_scc) {
    EXPECT_EQ((vector<vector<int>>{{ 9 }, { 8 }, { 10 }, { 7 }, { 5 }, { 6 }, { 3 }, { 2 }, { 1 }, { 0 }, { 4 } }), tarjan_scc(dag_neg1));
    EXPECT_EQ((vector<vector<int>>{{ 7 }, { 1, 3, 2, 0 }, { 6, 5, 4 } }), tarjan_scc(cyc_neg1));
}

TEST(graph_algorithms_test, chain_decomposition) {
    const auto d1 = chain_decomposition(cyc_undir);
    const auto ve1 = cut_edges(cyc_undir, d1);
    EXPECT_EQ((chain_decomposition_t{{ { { 0, 7, 5, 0 }, { 0, 8, 7 }, { 5, 9, 7 } }, { { 9, 3, 6, 9 } }, { { 4, 2, 1, 4 } } }, {}, { { { 13, 15, 14, 13 } }, { { 15, 17, 16, 15 } } }}), d1);
    EXPECT_EQ((vector<full_edge>{{ 4, 9 }, { 10, 11 }, { 11, 12 } }), sorted(ve1));
    EXPECT_EQ((vector<int>{ 4, 9, 11, 15 }), sorted(cut_vertices(cyc_undir, d1, ve1)));
    EXPECT_EQ((vector<vector<int>>{{ 0, 7, 5, 8, 9 }, { 4, 2, 1 }, { 9, 3, 6 }, { 13, 15, 14 }, { 15, 17, 16 } }), sorted(biconnected_components(cyc_undir, d1)));
    const auto d2 = chain_decomposition(cyc_undir2);
    EXPECT_EQ((chain_decomposition_t{ { { { 0, 3, 2, 1, 0 }, { 0, 4, 5, 3 }, { 1, 4 } } } }), d2);
    const auto d3 = chain_decomposition(cyc_undir3);
    EXPECT_EQ((vector<vector<int>>{{ { 0, 3, 2, 1, 4 }, { 1, 6, 7 } } }), sorted(biconnected_components(cyc_undir3, d3)));
    const auto d4 = chain_decomposition(cyc_undir4);
    EXPECT_EQ((vector<vector<int>>{{ 0, 3, 2, 1, 5, 4 }, { 1, 7, 6 }}), sorted(biconnected_components(cyc_undir4, d4)));
    const auto d5 = chain_decomposition(cyc_undir5);
    EXPECT_EQ((vector<vector<int>>{{ 0, 7, 4, 1, 11, 10 }, { 1, 3, 2 }, { 4, 6, 5 }, { 7, 9, 8 }}), sorted(biconnected_components(cyc_undir5, d5)));
    
    const auto d6 = chain_decomposition(blocks1);
    const auto ve6 = cut_edges(blocks1, d6);
    const auto va6 = cut_vertices(blocks1, d6, ve6);
    const auto vb6 = biconnected_components(blocks1, d6);
    EXPECT_EQ((vector<full_edge>{{ 4, 13 }, { 7, 11 }, { 7, 12 }, { 10, 14 }, { 14, 15 }, { 14, 16 }, { 17, 23 }, { 18, 20 }, { 19, 21 }, { 21, 22 } }), sorted(ve6));
    EXPECT_EQ((vector<int>{ 2, 3, 4, 7, 10, 14, 17, 18, 19, 21, 23 }), sorted(va6));
    EXPECT_EQ((vector<vector<int>>{{ 0, 3, 2, 1 }, { 2, 8, 7, 6 }, { 3, 5, 4 }, { 3, 10, 9 }, { 17, 19, 18 }, { 23, 25, 24 }}), sorted(vb6));
    EXPECT_EQ(make_pair(blocks1_bctree, blocks1_map), block_cut_tree(blocks1, ve6, va6, vb6));
}

TEST(graph_algorithms_test, transitive_closure) {
    EXPECT_EQ(graph<edge>({{ 4 }, { 4, 0 }, { 0, 4 }, { 2, 0, 1, 4 }, {}, { 1, 6, 4, 0 }, {}, { 6 }, { 10 }, {}, {} }), transitive_closure(dag_neg1));
    EXPECT_EQ(graph<edge>({ { 4 }, { 0 }, { 0 }, { 1, 2 }, {}, { 1, 6 }, {}, { 6 }, { 10 }, {}, {} }), transitive_reduction(dag_neg1));
    EXPECT_EQ(graph<edge>({ { 4 }, { 0 }, { 0 }, { 1, 2 }, {}, { 1, 6 }, {}, { 6 }, { 10 }, {}, {} }), transitive_reduction(transitive_closure(dag_neg1)));
    EXPECT_EQ(graph<edge>({ { 2, 3, 1, 4, 5, 6 }, { 0, 2, 3, 4, 5, 6 }, { 3, 1, 4, 5, 6, 0 }, { 1, 4, 5, 6, 0, 2 }, { 5, 6 }, { 6, 4 }, { 4, 5 }, { 5, 6, 4 } }), transitive_closure(cyc_neg1));
}

TEST(graph_algorithms_test, floyd_warshall) {
    int INF = 1000000000;
    EXPECT_EQ((vector<vector<iedge>> {
        { { 0, 0 }, { 2, -1 }, { 2, -2 }, { 2, 0 }, { 2, -8 }, { 2, -6 }, { 2, -3 }, { -1, INF } },
        { { 0, 4 }, { 1, 0 }, { 0, 2 }, { 0, 4 }, { 0, -4 }, { 0, -2 }, { 0, 1 }, { -1, INF } },
        { { 3, 5 }, { 3, 1 }, { 2, 0 }, { 3, 2 }, { 3, -6 }, { 3, -4 }, { 3, -1 }, { -1, INF } },
        { { 1, 3 }, { 1, -1 }, { 1, 1 }, { 3, 0 }, { 4, -8 }, { 4, -6 }, { 4, -3 }, { -1, INF } },
        { { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { 4, 0 }, { 5, 2 }, { 5, 5 }, { -1, INF } },
        { { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { 6, 10 }, { 5, 0 }, { 6, 3 }, { -1, INF } },
        { { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { 4, 7 }, { 4, 9 }, { 6, 0 }, { -1, INF } },
        { { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { 5, 16 }, { 5, 6 }, { 5, 9 }, { 7, 0 } }
    }), floyd_warshall(cyc_neg1, INF));
}

TEST(graph_algorithms_test, dijkstra) {
    int INF = 1000000000;
    EXPECT_EQ((vector<iedge>{{ 0, 0 }, { 3, 5 }, { 0, 2 }, { 2, 4 }, { 3, 12 }, { 4, 14 }, { 5, 17 }, { -1, INF } }), dijkstra(cyc_pos1, 0, INF));
    EXPECT_EQ((vector<iedge>{{ 1, 4 }, { 1, 0 }, { 1, 3 }, { 2, 5 }, { 3, 13 }, { 4, 15 }, { 5, 18 }, { -1, INF } }), dijkstra(cyc_pos1, 1, INF));
    EXPECT_EQ((vector<iedge>{{ 1, 7 }, { 3, 3 }, { 2, 0 }, { 2, 2 }, { 3, 10 }, { 4, 12 }, { 5, 15 }, { -1, INF } }), dijkstra(cyc_pos1, 2, INF));
    EXPECT_EQ((vector<iedge>{{ 1, 5 }, { 3, 1 }, { 1, 4 }, { 3, 0 }, { 3, 8 }, { 4, 10 }, { 5, 13 }, { -1, INF } }), dijkstra(cyc_pos1, 3, INF));
    EXPECT_EQ((vector<iedge>{{ -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { 4, 0 }, { 4, 2 }, { 5, 5 }, { -1, INF } }), dijkstra(cyc_pos1, 4, INF));
    EXPECT_EQ((vector<iedge>{{ -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { 6, 10 }, { 5, 0 }, { 5, 3 }, { -1, INF } }), dijkstra(cyc_pos1, 5, INF));
    EXPECT_EQ((vector<iedge>{{ -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { 6, 7 }, { 4, 9 }, { 6, 0 }, { -1, INF } }), dijkstra(cyc_pos1, 6, INF));
    EXPECT_EQ((vector<iedge>{{ -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { 6, 16 }, { 7, 6 }, { 5, 9 }, { 7, 0 } }), dijkstra(cyc_pos1, 7, INF));
}

TEST(graph_algorithms_test, prim_spanning_tree) {
    int INF = 1000000000;
    EXPECT_EQ((vector<iedge>{
        { 0, 0 }, { 4, 33 }, { 1, 31 }, { 9, 26 }, { 9, 34 }, { 0, 21 }, { 9, 28 }, { 5, 26 }, { 0, 23 }, { 5, 29 },
        { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF } }), prim_spanning_tree(cyc_undir, 0, INF));
    EXPECT_EQ((vector<iedge>{
        { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF },
        { 10, 0 }, { 10, 45 }, { 11, 38 }, { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF } }), prim_spanning_tree(cyc_undir, 10, INF));
    EXPECT_EQ((vector<iedge>{
        { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF }, { -1, INF },
        { -1, INF }, { -1, INF }, { -1, INF }, { 14, 55 }, { 15, 54 }, { 17, 50 }, { 16, 0 }, { 16, 55 } }), prim_spanning_tree(cyc_undir, 16, INF));
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

TEST(graph_algorithms_test, dinic_flow) {
    test_max_flow<dinic_flow<int>, int>({ { 0 } }, { { 0 } });
    test_max_flow<dinic_flow<int>, int>({ { 0, 5 }, { 7, 0 } }, { { 0, 5 }, { 7, 0 } });
    test_max_flow<dinic_flow<int>, int>({ { 0, 3, 5 }, { 0, 0, 2 }, { 0, 0, 0 } }, { { 0, 3, 7 }, { 0, 0, 2 }, { 0, 0, 0 } });
    test_max_flow<dinic_flow<double>, double>({ { 0, 5, 2 }, { 7, 0, 4 }, { 1, 3, 0 } }, { { 0, 7, 6 }, { 8, 0, 6 }, { 4, 4, 0 } });
}

TEST(graph_algorithms_test, push_relabel_flow) {
    test_max_flow<push_relabel_flow<int>, int>({ { 0 } }, { { 0 } });
    test_max_flow<push_relabel_flow<int>, int>({ { 0, 5 }, { 7, 0 } }, { { 0, 5 }, { 7, 0 } });
    test_max_flow<push_relabel_flow<int>, int>({ { 0, 3, 5 }, { 0, 0, 2 }, { 0, 0, 0 } }, { { 0, 3, 7 }, { 0, 0, 2 }, { 0, 0, 0 } });
    test_max_flow<push_relabel_flow<double>, double>({ { 0, 5, 2 }, { 7, 0, 4 }, { 1, 3, 0 } }, { { 0, 7, 6 }, { 8, 0, 6 }, { 4, 4, 0 } });
}

TEST(graph_algorithms_test, bipartite_matching) {
    typedef std::vector<full_edge> edges_t;
    EXPECT_EQ((edges_t()), bipartite_matching(0, edges_t()));
    EXPECT_EQ((edges_t{ { 0, 2 }, { 1, 3 } }), bipartite_matching(4, edges_t{ { 0, 2 }, { 0, 3 }, { 1, 3 } }));
    EXPECT_EQ((edges_t{ { 0, 2 }, { 1, 3 } }), bipartite_matching(4, edges_t{ { 0, 2 }, { 1, 2 }, { 1, 3 } }));
}

TEST(graph_algorithms_test, lowest_common_ancestor) {
    graph<edge> t({ { 1, 2 }, { 0 }, { 0, 3 }, { 2 } });
    lowest_common_ancestor lca(t);
    EXPECT_EQ(0, lca.ancestor(1, 3));
}

TEST(graph_algorithms_test, heavy_light_decomposition) {
    graph<edge> t({ { 1, 2 }, { 0 }, { 0, 3 }, { 2 } });
    heavy_light_decomposition_ex hld(t);
    EXPECT_EQ(3, hld.parent(3, 0));
    EXPECT_EQ(2, hld.parent(3, 1));
    EXPECT_EQ(0, hld.parent(3, 2));
}

namespace {
    auto path_graph = [](int n) {
        graph<edge> g(n);
        for (int i = 1; i < n; i++) {
            g.add_edge2(i - 1, i);
        }
        return g;
    };
    auto tree_graph = [](int n) {
        graph<edge> g(n);
        for (int i = 1; i < n; i++) {
            int p = rand() % i;
            g.add_edge2(i, p);
        }
        return g;
    };
    auto cycle_graph = [](int n) {
        auto g = path_graph(n);
        if (n > 2) g.add_edge2(n - 1, 0);
        return g;
    };
    auto complete_graph = [](int n) {
        graph<edge> g(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                g.add_edge2(i, j);
            }
        }
        return g;
    };
    auto petersen_graph = []() {
        graph<edge> g(10);
        vector<full_edge> ve{
            { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 4 }, { 4, 0 },
            { 0, 5 }, { 1, 6 }, { 2, 7 }, { 3, 8 }, { 4, 9 },
            { 5, 7 }, { 6, 8 }, { 7, 9 }, { 8, 5 }, { 9, 6 },
        };
        for (auto e : ve) g.add_edge2(e.u, e.v);
        return g;
    };
    typedef polynom<int> poly;
}

TEST(graph_algorithms_test, chromatic_polynomial) {
    for (int n = 1; n < 20; n++) {
        auto p = chromatic_polynomial(path_graph(n), 1);
        auto p0 = chromatic_polynomial_T(n, 1);
        //cout << "P_" << n << ": " << "  " << clock() << " ms" << endl;
        EXPECT_EQ(p0, p) << "ERROR: P_" << n << ": " << p << endl;
    }
    for (int n = 1; n < 20; n++) {
        auto p = chromatic_polynomial(tree_graph(n), 1);
        auto p0 = chromatic_polynomial_T(n, 1);
        //cout << "T_" << n << ": " << "  " << clock() << " ms" << endl;
        EXPECT_EQ(p0, p) << "ERROR: T_" << n << ": " << p << endl;
    }
    for (int n = 3; n < 20; n++) {
        auto p = chromatic_polynomial(cycle_graph(n), 1);
        auto p0 = chromatic_polynomial_C(n, 1);
        //cout << "C_" << n << ": " << "  " << clock() << " ms" << endl;
        EXPECT_EQ(p0, p) << "ERROR: C_" << n << ": " << p << endl;
    }
    for (int n = 1; n < 20; n++) {
        auto p = chromatic_polynomial(complete_graph(n), 1);
        auto p0 = chromatic_polynomial_K(n, 1);
        //cout << "K_" << n << ": " << "  " << clock() << " ms" << endl;
        EXPECT_EQ(p0, p) << "ERROR: K_" << n << ": " << p << endl;
    }
    if (true) {
        auto p = chromatic_polynomial(cyc_undir2, 1);
        // k (k - 1) (k - 2)^2 (k (k - 3) + 4)
        auto p0 = poly{ -8, 10, -5, 1 } *chromatic_polynomial_K(3, 1);
        //cout << "Petersen related" << ": " << "  " << clock() << " ms" << endl;
        EXPECT_EQ(p0, p) << "ERROR: Petersen related" << ": " << p << endl;
    }
    if (true) {
        auto p = chromatic_polynomial(petersen_graph(), 1);
        auto p0 = poly{ -352, 775, -814, 529, -230, 67, -12, 1 } * chromatic_polynomial_K(3, 1);
        //cout << "Petersen" << ": " << "  " << clock() << " ms" << endl;
        EXPECT_EQ(p0, p) << "ERROR: Petersen" << ": " << p << endl;
    }
}

TEST(graph_algorithms_test, sat2) {
    // (x0 || -x1) (-x0 || x2)
    vector<pair<int, int>> clauses1{ { 2 * 0 + 1, 2 * 1 + 0 }, { 2 * 0 + 0, 2 * 2 + 1 } };
    EXPECT_EQ((vector<int>{0, 0, 1}), sat2(clauses1));

    vector<int> sol2{ -1, -1, -1 };
    EXPECT_EQ(true, sat2(sol2, clauses1));
    EXPECT_EQ((vector<int>{0, 0, 1}), sol2);

    vector<int> sol3{ -1, 1, -1 };
    EXPECT_EQ(true, sat2(sol3, clauses1));
    EXPECT_EQ((vector<int>{1, 1, 1}), sol3);

    // (x0 || x1) (x0 || -x1) (-x0 || x1) (-x0 || -x1)
    vector<pair<int, int>> clauses2{ { 2 * 0 + 1, 2 * 1 + 1 }, { 2 * 0 + 1, 2 * 1 + 0 }, { 2 * 0 + 0, 2 * 1 + 1 }, { 2 * 0 + 0, 2 * 1 + 0 } };
    EXPECT_EQ((vector<int>()), sat2(clauses2));

}
