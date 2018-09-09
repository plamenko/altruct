#pragma once

#include "altruct/structure/graph/graph.h"

#include <vector>

namespace altruct {
namespace graph {

/**
 * Performs DFS and invokes `visitor` on each step.
 *
 * If `source` is specified only the component rooted at `source` will be traversed.
 * Otherwise, each component will be traversed (the graph may be disconnected).
 *
 * Complexity: O(m)
 *
 * @param g - a graph with `n` nodes and `m` edges
 * @param visitor - a function `bool visitor(int root, int parent, int node, int depth)`
 *                  @param root - the root node of the current component
 *                  @param parent - the parent of the current node,
 *                                -1 if the current node is the root
 *                                (visited in euler-tour-order)
 *                  @param node - the current node; an artificial -1 is
 *                                reported as the last node of each parent
 *                                (visited in pre-order)
 *                  @param depth - the depth of the current node
 *                  @return - true if the node should be entered
 * @param source - the source node, or -1 if not specified
 */
template<typename F, typename E>
void iterative_dfs(const graph<E>& g, F visitor, int source = -1) {
    std::vector<int> visited(g.size());
    std::vector<std::pair<int, int>> stk;
    // this outer loop allows us to handle a disconnected graph
    int of = (source != -1) ? source : 0;
    int ol = (source != -1) ? source : g.size() - 1;
    for (int o = of; o <= ol; o++) {
        if (visited[o]) continue;
        if (visitor(o, -1, o, 0)) {
            visited[o] = 1;
            stk.push_back({ o, 0 });
        }
        while (!stk.empty()) {
            int d = (int)stk.size();
            int u = stk.back().first;
            int i = stk.back().second;
            if (i < (int)g[u].size()) {
                int v = g[u][i].v;
                stk.back().second++;
                if (visited[v]) continue;
                if (visitor(o, u, v, d)) {
                    visited[v] = 1;
                    stk.push_back({ v, 0 });
                }
            } else {
                stk.pop_back();
                visitor(o, u, -1, d);
            }
        }
    }
}

template<typename E>
std::vector<int> parents(const graph<E>& g) {
    std::vector<int> vp(g.size(), -1);
    iterative_dfs(g, [&](int root, int parent, int node, int depth) {
        if (node != -1) vp[node] = parent;
        return true;
    });
    return vp;
}

} // graph
} // altruct
