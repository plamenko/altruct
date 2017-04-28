#pragma once

#include "graph.h"

#include <set>
#include <vector>
#include <algorithm>

namespace altruct {
namespace graph {

/**
 * Calculates the distances from the source to all other vertices using Dijkstra algorithm.
 *
 * Complexity: O(m log n)
 *
 * @param g - an undirected graph with `n` nodes and `m` edges
 * @param src - the source node to calculate the shortest path tree from
 * @return `res` where `res[v] = { i, d }` and:
 *   `i` is the previous node before `v` on the shortest path from `src` to `v`,
 *   `d` is the shortest distance between `src` and `v`.
 */
template<typename E, typename W>
std::vector<E> dijkstra(const graph<E>& g, int src, W inf) {
    std::vector<E> res(g.size(), { -1, inf });
    res[src] = { src, 0 };
    std::set<std::pair<W, int>> q{ { 0, src } };
    while (!q.empty()) {
        int u = q.begin()->second;
        q.erase(q.begin());
        for (const auto& e : g[u]) {
            int v = e.v;
            auto d_new = res[u].w + e.w;
            if (d_new < res[v].w) {
                q.erase({ res[v].w, v });
                res[v] = { u, d_new };
                q.insert({ res[v].w, v });
            }
        }
    }
    return res;
}

} // graph
} // altruct
