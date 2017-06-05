#pragma once

#include "altruct/structure/graph/graph.h"

#include <set>
#include <vector>
#include <algorithm>

namespace altruct {
namespace graph {

/**
 * Calculates the minimum spanning tree using Prim's algorithm.
 *
 * Complexity: O(m log n)
 *
 * @param g - an undirected graph with `n` nodes and `m` edges
 * @param root - the root node to calculate the minimum spanning tree from
 * @return `res` where `res[v] = { p, d }` and:
 *   `p` is the parent node of `v`,
 *   `d` is the distance from `p` to `v`.
 */
template<typename E, typename W>
std::vector<E> prim_spanning_tree(const graph<E>& g, int root, W inf) {
    std::vector<E> res(g.size(), { -1, inf });
    res[root] = { root, 0 };
    std::set<std::pair<W, int>> q{ { 0, root } };
    while (!q.empty()) {
        int u = q.begin()->second;
        q.erase(q.begin());
        for (const auto& e : g[u]) {
            int v = e.v;
            auto d_new = e.w;
            if (d_new < res[v].w && res[u].v != v) {
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
