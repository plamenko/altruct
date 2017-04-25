#pragma once

#include <vector>
#include <algorithm>

namespace altruct {
namespace graph {

/**
 * Calculates the transitive closure of a graph.
 *
 * Complexity: O(n*m)
 *
 * @param adjl - adjacency list of `n` nodes and `m` edges
 * @param index_f - a functor that extracts the outbound node index from the edge
 */
template<typename E, typename FI>
std::vector<std::vector<int>> transitive_closure(const std::vector<std::vector<E>>& adjl, FI index_f) {
    std::vector<std::vector<int>> res(adjl.size());
    for (int i = 0; i < (int)adjl.size(); i++) {
        std::vector<int> d(adjl.size(), 0);
        vector<int> stk{ i };
        d[i] = -1;
        while (!stk.empty()) {
            int u = stk.back(); stk.pop_back();
            for (const auto& e : adjl[u]) {
                int v = index_f(e);
                if (d[v] != 0) continue;
                stk.push_back(v);
                res[i].push_back(v);
                d[v] = 1;
            }
        }
    }
    return res;
}
template<typename E = int>
std::vector<std::vector<int>> transitive_closure(const std::vector<std::vector<int>>& adjl) {
    return transitive_closure(adjl, [](int i){ return i; });
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
 * @param adjl - adjacency list of `n` nodes and `m` edges
 * @param index_f - a functor that extracts the outbound node index from the edge
 */
template<typename E, typename FI>
std::vector<std::vector<int>> transitive_reduction(const std::vector<std::vector<E>>& adjl, FI index_f) {
    std::vector<std::vector<int>> res(adjl.size());
    for (int i = 0; i < (int)adjl.size(); i++) {
        std::vector<int> d(adjl.size(), 0);
        vector<int> stk;
        d[i] = -1;
        for (const auto& e : adjl[i]) {
            int v = index_f(e);
            if (d[v] == 0) stk.push_back(v);
            d[v] = 1;
        }
        d[i] = -1;
        while (!stk.empty()) {
            int u = stk.back(); stk.pop_back();
            for (const auto& e : adjl[u]) {
                int v = index_f(e);
                if (d[v] == 0) stk.push_back(v);
                d[v] = -1; // reachable indirectly, hence not in the reduction
            }
        }
        for (int v = 0; v < (int)adjl.size(); v++) {
            if (d[v] == 1) res[i].push_back(v);
        }
    }
    return res;
}
template<typename E = int>
std::vector<std::vector<int>> transitive_reduction(const std::vector<std::vector<int>>& adjl) {
    return transitive_reduction(adjl, [](int i){ return i; });
}

} // graph
} // altruct
