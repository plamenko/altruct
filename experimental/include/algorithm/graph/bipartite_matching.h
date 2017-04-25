#pragma once

#include <vector>
#include <algorithm>

#include "algorithm/graph/dinic_flow.h"

namespace altruct {
namespace graph {

/**
 * Bipartite matching.
 *
 * Uses Dinic's maximum flow lgorithm.
 */
std::vector<std::pair<int, int>> bipartite_matching(int nodes, const std::vector<std::pair<int, int>>& edges) {
    int src = nodes++;
    int sink = nodes++;
    std::vector<std::vector<int>> cap(nodes, std::vector<int>(nodes));
    for (const auto& e : edges) {
        cap[src][e.first] = 1;
        cap[e.first][e.second] = 1;
        cap[e.second][sink] = 1;
    }
    dinic_flow<int> d(cap);
    d.calc_max_flow(src, sink);
    std::vector<std::pair<int, int>> selected;
    for (const auto& e : edges) {
        if (d.flow[e.first][e.second]) selected.push_back(e);
    }
    return selected;
}

} // graph
} // altruct
