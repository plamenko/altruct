#pragma once

#include "altruct/structure/graph/graph.h"

#include <vector>
#include <algorithm>

#include "altruct/algorithm/graph/dinic_flow.h"

namespace altruct {
namespace graph {

/**
 * Bipartite matching.
 *
 * Uses Dinic's maximum flow lgorithm.
 */
template<typename E = full_edge>
std::vector<E> bipartite_matching(int nodes, const std::vector<E>& edges) {
    int src = nodes++;
    int sink = nodes++;
    std::vector<std::vector<int>> cap(nodes, std::vector<int>(nodes));
    for (const auto& e : edges) {
        cap[src][e.u] = 1;
        cap[e.u][e.v] = 1;
        cap[e.v][sink] = 1;
    }
    dinic_flow<int> d(cap);
    d.calc_max_flow(src, sink);
    std::vector<E> selected;
    for (const auto& e : edges) {
        if (d.flow[e.u][e.v]) selected.push_back(e);
    }
    return selected;
}

} // graph
} // altruct
