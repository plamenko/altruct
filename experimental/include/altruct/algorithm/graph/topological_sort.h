#pragma once

#include "altruct/structure/graph/graph.h"

#include "altruct/algorithm/graph/iterative_dfs.h"

#include <vector>
#include <algorithm>

namespace altruct {
namespace graph {

/**
 * Calculates in-degrees for each node.
 *
 * Complexity: O(m)
 *
 * @param g - a graph with `n` nodes and `m` edges
 */
template<typename E>
std::vector<int> in_degrees(const graph<E>& g) {
    std::vector<int> deg(g.size());
    for (int i = 0; i < g.size(); i++) {
        for (const auto& e : g[i]) {
            deg[e.v]++;
        }
    }
    return deg;
}

/**
 * Calculates topological order of DAG nodes.
 *
 * Complexity: O(m)
 *
 * @param g - a graph with `n` nodes and `m` edges
 */
template<typename E>
std::vector<int> topological_sort(const graph<E>& g) {
    auto deg = in_degrees(g);
    std::vector<int> topo;
    iterative_dfs(g, [&](int root, int parent, int node, int depth) {
        if (parent == -1 && deg[node] > 0) return false; // not a proper root
        if (node < 0) topo.push_back(parent);
        return true;
    });
    reverse(topo.begin(), topo.end());
    return topo;
}

} // graph
} // altruct
