#pragma once

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
 * @param adjl - adjacency list of `n` nodes and `m` edges
 * @param index_f - a functor that extracts the edge's outbound node index
 * @param weight_f - a functor that extracts the edge's weight
 * @return `res` where `res[u][v] = { i, d }` and:
 *   `i` is the next node after `u` on the shortest path from `u` to `v`,
 *   `d` is the shortest distance between `u` and `v`.
 */
template<typename E, typename FI, typename FW, typename W = typename std::result_of<FW(E)>::type>
std::vector<std::vector<std::pair<int, W>>> floyd_warshall(const std::vector<std::vector<E>>& adjl, FI index_f, FW weight_f, W inf) {
    int n = (int)adjl.size();
    typedef std::vector<std::pair<int, W>> row_type;
    std::vector<row_type> res(n, row_type(n, { -1, inf }));
    for (int i = 0; i < n; i++) {
        res[i][i] = { i, 0 };
        for (const auto& e : adjl[i]) {
            int v = index_f(e);
            res[i][v] = { v, weight_f(e) };
        }
    }
    for (int i = 0; i < n; i++) {
        for (int u = 0; u < n; u++) {
            for (int v = 0; v < n; v++) {
                // try to improve the path `u -> v` by going through `i`, i.e.: `u -> i -> v`
                if (res[u][i].first == -1 || res[i][v].first == -1) {
                    continue; // not connected
                }
                auto d_new = res[u][i].second + res[i][v].second;
                if (d_new < res[u][v].second) {
                    res[u][v] = { res[u][i].first, d_new };
                }
            }
        }
    }
    return res;
}
template<typename W>
std::vector<std::vector<std::pair<int, W>>> floyd_warshall(const std::vector<std::vector<std::pair<int, W>>>& adjl, W inf) {
    auto index_f = [](const std::pair<int, W>& e){ return e.first; };
    auto weight_f = [](const std::pair<int, W>& e){ return e.second; };
    return floyd_warshall(adjl, index_f, weight_f, inf);
}

} // graph
} // altruct
