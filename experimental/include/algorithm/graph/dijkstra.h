#pragma once

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
 * @param adjl - adjacency list of `n` nodes and `m` edges
 * @param src - the source node to calculate the shortest path tree from
 * @param index_f - a functor that extracts the edge's outbound node index
 * @param weight_f - a functor that extracts the edge's weight
 * @return `res` where `res[v] = { i, d }` and:
 *   `i` is the previous node before `v` on the shortest path from `src` to `v`,
 *   `d` is the shortest distance between `src` and `v`.
 */
template<typename E, typename FI, typename FW, typename W = typename std::result_of<FW(E)>::type>
std::vector<std::pair<int, W>> dijkstra(const std::vector<std::vector<E>>& adjl, int src, FI index_f, FW weight_f, W inf) {
    std::vector<std::pair<int, W>>d(adjl.size(), { -1, inf });
    d[src] = { src, 0 };
    std::set<std::pair<W, int>> q{ { 0, src } };
    while (!q.empty()) {
        int u = q.begin()->second;
        q.erase(q.begin());
        for (const auto& e : adjl[u]) {
            int v = index_f(e);
            auto d_new = d[u].second + weight_f(e);
            if (d_new < d[v].second) {
                q.erase({ d[v].second, v });
                d[v] = { u, d_new };
                q.insert({ d[v].second, v });
            }
        }
    }
    return d;
}
template<typename W>
std::vector<std::pair<int, W>> dijkstra(const std::vector<std::vector<std::pair<int, W>>>& adjl, int src, W inf) {
    auto index_f = [](const std::pair<int, W>& e){ return e.first; };
    auto weight_f = [](const std::pair<int, W>& e){ return e.second; };
    return dijkstra(adjl, src, index_f, weight_f, inf);
}

} // graph
} // altruct
