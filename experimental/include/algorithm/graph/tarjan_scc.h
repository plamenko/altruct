#pragma once

#include "structure/graph/graph.h"

#include <vector>
#include <algorithm>

namespace altruct {
namespace graph {

/**
 * Calculates the set of strongly connected components in topological order.
 *
 * Complexity: O(m)
 *
 * @param g - a directed acyclic graph with `n` nodes and `m` edges
 */
template<typename E>
std::vector<std::vector<int>> tarjan_scc(const graph<E>& g) {
    std::vector<std::vector<int>> vscc;
    int cnt = 0;
    std::vector<int> low(g.size(), -1); // node's lowest reachable ancestor
    std::vector<int> idx(g.size(), -1); // node's index in the dfs traversal
    std::vector<int> act(g.size(), -1); // whether or not a node is on `sta`
    std::vector<int> sta;                  // stack of active nodes
    std::vector<std::pair<int, int>> stk;  // dfs stack
    for (int o = 0; o < g.size(); o++) {
        if (idx[o] != -1) continue;
        act[o] = low[o] = idx[o] = cnt++;
        sta.push_back(o);
        stk.push_back({ o, 0 });
        while (!stk.empty()) {
            int u = stk.back().first;
            int i = stk.back().second;
            if (i < (int)g[u].size()) {
                int v = g[u][i].v;
                stk.back().second++;
                if (act[v] != -1) {
                    low[u] = std::min(low[u], idx[v]);
                }
                if (idx[v] != -1) continue;
                act[v] = low[v] = idx[v] = cnt++;
                sta.push_back(v);
                stk.push_back({ v, 0 });
            } else {
                if (low[u] == idx[u]) {
                    std::vector<int> scc;
                    int v = -1;
                    do {
                        v = sta.back();
                        sta.pop_back();
                        act[v] = -1;
                        scc.push_back(v);
                    } while (v != u);
                    vscc.push_back(scc);
                }
                stk.pop_back();
                if (!stk.empty()) {
                    int v = u; u = stk.back().first;
                    low[u] = std::min(low[u], low[v]);
                }
            }
        }
    }
    std::reverse(vscc.begin(), vscc.end());
    return vscc;
}

} // graph
} // altruct
