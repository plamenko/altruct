#pragma once

#include "structure/graph/graph.h"

#include <vector>
#include <algorithm>

namespace altruct {
namespace graph {

/**
 * Calculates the transitive closure of a graph.
 *
 * Complexity: O(n*m)
 *
 * @param g - a graph with `n` nodes and `m` edges
 */
template<typename E>
graph<edge> transitive_closure(const graph<E>& g) {
    graph<edge> res(g.size());
    for (int i = 0; i < g.size(); i++) {
        std::vector<int> d(g.size(), 0);
        std::vector<int> stk{ i };
        d[i] = -1;
        while (!stk.empty()) {
            int u = stk.back(); stk.pop_back();
            for (const auto& e : g[u]) {
                int v = e.v;
                if (d[v] != 0) continue;
                stk.push_back(v);
                res.add_edge(i, v);
                d[v] = 1;
            }
        }
    }
    return res;
}

/**
 * Calculates the transitive reduction of an acyclic graph (DAG).
 *
 * Note, for a graph with cycles, first compute its condensation
 * (each storngly-connected-component is condensed to a single node)
 * and then compute a transitive reduction of the resulting acyclic
 * graph with this method.
 * Each SCC can then be reduced independently by finding a cycle that
 * visits all of its nodes. However, to find a shortest such cycle
 * (i.e. smallest subgraph with the same reachability) is NP-hard.
 *
 * Complexity: O(n*m)
 *
 * @param g - an undirected graph with `n` nodes and `m` edges
 */
template<typename E>
graph<edge> transitive_reduction(const graph<E>& g) {
    graph<edge> res(g.size());
    for (int i = 0; i < g.size(); i++) {
        std::vector<int> d(g.size(), 0);
        std::vector<int> stk;
        d[i] = -1;
        for (const auto& e : g[i]) {
            int v = e.v;
            if (d[v] == 0) stk.push_back(v);
            d[v] = 1;
        }
        d[i] = -1;
        while (!stk.empty()) {
            int u = stk.back(); stk.pop_back();
            for (const auto& e : g[u]) {
                int v = e.v;
                if (d[v] == 0) stk.push_back(v);
                d[v] = -1; // reachable indirectly, hence not in the reduction
            }
        }
        for (int v = 0; v < g.size(); v++) {
            if (d[v] == 1) res.add_edge(i, v);
        }
    }
    return res;
}

} // graph
} // altruct
