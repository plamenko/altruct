#pragma once

#include "altruct/structure/graph/graph.h"

#include <vector>
#include <algorithm>

namespace altruct {
namespace graph {

/**
 * Calculates the distance between every pair of vertices using Floyd-Warshall algorithm.
 *
 * Note, the algorithm works for negative weights so long as there are no negative cycles.
 *
 * Complexity: O(n^3)
 *
 * @param g - an undirected graph with `n` nodes and `m` edges
 * @return `res` where `res[u][v] = { i, d }` and:
 *   `i` is the next node after `u` on the shortest path from `u` to `v`,
 *   `d` is the shortest distance between `u` and `v`.
 */
template<typename E, typename W>
std::vector<std::vector<E>> floyd_warshall(const graph<E>& g, W inf) {
    int n = g.size();
    typedef std::vector<E> row_type;
    std::vector<row_type> res(n, row_type(n, { -1, inf }));
    for (int i = 0; i < n; i++) {
        res[i][i] = { i, 0 };
        for (const auto& e : g[i]) {
            int v = e.v;
            res[i][v] = { v, std::min(res[i][v].w, e.w) };
        }
    }
    for (int i = 0; i < n; i++) {
        for (int u = 0; u < n; u++) {
            for (int v = 0; v < n; v++) {
                // try to improve the path `u -> v` by going through `i`, i.e.: `u -> i -> v`
                if (res[u][i].v == -1 || res[i][v].v == -1) {
                    continue; // not connected
                }
                auto d_new = res[u][i].w + res[i][v].w;
                if (d_new < res[u][v].w) {
                    res[u][v] = { res[u][i].v, d_new };
                }
            }
        }
    }
    return res;
}

} // graph
} // altruct
