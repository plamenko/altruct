#pragma once

#include "algorithm/graph/iterative_dfs.h"

#include <vector>
#include <algorithm>

namespace altruct {
namespace graph {

/**
 * Calculates in-degrees for each node.
 *
 * Complexity: O(m)
 *
 * @param adjl - adjacency list of `n` nodes and `m` edges
 * @param index_f - a functor that extracts the outbound node index from the edge
 */
template<typename E, typename FI>
std::vector<int> in_degrees(const std::vector<std::vector<E>>& adjl, FI index_f) {
    std::vector<int> deg(adjl.size());
    for (int i = 0; i < (int)adjl.size(); i++) {
        for (const auto& e : adjl[i]) {
            deg[index_f(e)]++;
        }
    }
    return deg;
}
template<typename E = int>
std::vector<int> in_degrees(const std::vector<std::vector<int>>& adjl) {
    return in_degrees(adjl, visitor, [](int i){ return i; });
}

/**
 * Calculates topological order of DAG nodes.
 *
 * Complexity: O(m)
 *
 * @param adjl - adjacency list of `n` nodes and `m` edges
 * @param index_f - a functor that extracts the outbound node index from the edge
 */
template<typename E, typename FI>
std::vector<int> topological_sort(const std::vector<std::vector<E>>& adjl, FI index_f) {
    auto deg = in_degrees(adjl, index_f);
    std::vector<int> topo;
    iterative_dfs(adjl, [&](int root, int parent, int node, int depth) {
        if (parent == -1 && deg[node] > 0) return false; // not a proper root
        if (node < 0) topo.push_back(parent);
        return true;
    }, index_f, -1);
    reverse(topo.begin(), topo.end());
    return topo;
}
template<typename E = int>
std::vector<int> topological_sort(const std::vector<std::vector<int>>& adjl) {
    return topological_sort(adjl, visitor, [](int i){ return i; });
}

} // graph
} // altruct
