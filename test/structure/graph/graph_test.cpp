#include "altruct/structure/graph/graph.h"

#include "altruct/io/iostream_overloads.h"
#include "structure_test_util.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::graph;
using namespace altruct::test_util;

namespace {
    typedef weighted_edge<int> iedge;
    graph<iedge> g({
        { { 2, -2 } },
        { { 0, 4 }, { 2, 3 } },
        { { 3, 2 } },
        { { 1, -1 }, { 4, -8 } },
        { { 5, 2 } },
        { { 6, 3 }, { 7, 10 } },
        { { 4, 7 } },
        { { 5, 10 }, { 1, 6 }, { 5, 11 } },
    });
}

TEST(graph_test, edge_types) {
    edge e1(4);
    EXPECT_EQ(4, e1.v);
    edge e2(7);
    EXPECT_EQ(7, e2.v);
    edge e3 = e2;
    ASSERT_BASIC_COMPARISON_OPERATORS(-1, e1, e2);
    ASSERT_BASIC_COMPARISON_OPERATORS(0, e2, e3);

    full_edge f1(5, 3);
    EXPECT_EQ(5, f1.u);
    EXPECT_EQ(3, f1.v);
    full_edge f2(5, 4);
    EXPECT_EQ(5, f2.u);
    EXPECT_EQ(4, f2.v);
    full_edge f3(2, 10);
    EXPECT_EQ(2, f3.u);
    EXPECT_EQ(10, f3.v);
    full_edge f4 = f3;
    ASSERT_BASIC_COMPARISON_OPERATORS(-1, f1, f2);
    ASSERT_BASIC_COMPARISON_OPERATORS(+1, f1, f3);
    ASSERT_BASIC_COMPARISON_OPERATORS(0, f4, f3);

    weighted_edge<double> w1(5, 3.5);
    EXPECT_EQ(5, w1.v);
    EXPECT_EQ(3.5, w1.w);
    weighted_edge<double> w2(5, 4.5);
    EXPECT_EQ(5, w2.v);
    EXPECT_EQ(4.5, w2.w);
    weighted_edge<double> w3(2, 10);
    EXPECT_EQ(2, w3.v);
    EXPECT_EQ(10, w3.w);
    weighted_edge<double> w4 = w3;
    ASSERT_BASIC_COMPARISON_OPERATORS(0, w1, w2); // weight not compared
    ASSERT_BASIC_COMPARISON_OPERATORS(+1, w1, w3);
    ASSERT_BASIC_COMPARISON_OPERATORS(0, w4, w3);
}

TEST(graph_test, constructor_and_size) {
    graph<iedge> g1;
    EXPECT_EQ(0, g1.size());
    EXPECT_EQ(0, g1.num_edges());

    graph<iedge> g2(10);
    EXPECT_EQ(10, g2.size());
    EXPECT_EQ(0, g2.num_edges());

    graph<iedge> g3({
        { { 2, 10 }, { 1, 50 } },
        {},
        { { 0, 100 }, { 1, 30 } },
    });
    EXPECT_EQ(3, g3.size());
    EXPECT_EQ(4, g3.num_edges());
}

TEST(graph_test, bracket_operator) {
    auto g1 = g;
    const auto g2 = g;
    EXPECT_EQ((vector<iedge>{{ 0, 4 }, { 2, 3 } }), g1[1]);
    EXPECT_EQ((vector<iedge>{{ 0, 4 }, { 2, 3 } }), g2[1]);
}

TEST(graph_test, comparison) {
    graph<iedge> g1 = g;
    graph<iedge> g2(8);
    ASSERT_BASIC_COMPARISON_OPERATORS(0, g, g1);
    ASSERT_BASIC_COMPARISON_OPERATORS(1, g, g2);
}

TEST(graph_test, mutation) {
    auto g1 = g;
    EXPECT_EQ(8, g1.size());
    EXPECT_EQ(13, g1.num_edges());
    
    EXPECT_EQ(8, g1.add_node());
    EXPECT_EQ(9, g1.size());
    EXPECT_EQ(13, g1.num_edges());

    g1.add_edge(1, { 8, 50 });
    EXPECT_EQ(9, g1.size());
    EXPECT_EQ(14, g1.num_edges());

    g1.add_edge2(8, { 3, 70 });
    EXPECT_EQ(9, g1.size());
    EXPECT_EQ(16, g1.num_edges());
    EXPECT_EQ(iedge(3, 70), g1[8].back());
    EXPECT_EQ(iedge(8, 70), g1[3].back());

    g1.add_edge(8, { 3, 60 }); // duplicate edge from 8 to 3
    EXPECT_EQ(9, g1.size());
    EXPECT_EQ(17, g1.num_edges());

    g1.delete_edge(8, 3); // deletes both edges from 8 to 3, but not from 3 to 8!
    EXPECT_EQ(9, g1.size());
    EXPECT_EQ(15, g1.num_edges());
    EXPECT_EQ(0, g1[8].size());
    EXPECT_EQ(iedge(8, 70), g1[3].back());

    auto g2 = g;
    g2.delete_node(1);
    graph<iedge> ge2({
        { { 2, -2 } },
        { { 5, 10 }, { 5, 11 } }, // the last node 7 gets moved to the deleted index 1
        { { 3, 2 } },
        { { 4, -8 } },
        { { 5, 2 } },
        { { 6, 3 }, { 1, 10 } }, // 7->1
        { { 4, 7 } },
    });
    EXPECT_EQ(ge2, g2);

    auto g3 = g;
    g3.contract(1, 3);
    graph<iedge> ge3({
        { { 2, -2 } },
        { { 0, 4 }, { 2, 3 }, { 4, -8 } },
        { { 1, 2 } },
        { { 1, 6 }, { 5, 10 } },
        { { 5, 2 } },
        { { 3, 10 }, { 6, 3 } },
        { { 4, 7 } },
    });
    EXPECT_EQ(ge3, g3);
}
