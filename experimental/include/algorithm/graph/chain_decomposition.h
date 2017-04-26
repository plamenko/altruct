#pragma once

#include <vector>
#include <unordered_set>
#include <algorithm>
#include "algorithm/hash/std_tuple_hash.h"

namespace altruct {
namespace graph {

/**
 * Calculates the chain decomposition of the given undirected graph.
 *
 * Note, for each edge {u, v}, there should be a corresponding edge {v, u}.
 *
 * Complexity: O(m)
 *
 * @param adjl - adjacency list of `n` nodes and `m` edges
 * @param index_f - a functor that extracts the outbound node index from the edge
 * @return - list of chains for each component
 */
template<typename E, typename FI>
std::vector<std::vector<std::vector<int>>> chain_decomposition(const std::vector<std::vector<E>>& adjl, FI index_f) {
    std::vector<int> que; 
    std::vector<int> ord(adjl.size(), -1);
    std::vector<int> par(adjl.size(), -1);
    iterative_dfs(adjl, [&](int root, int parent, int node, int depth) {
        if (node != -1) {
            que.push_back(node);
            ord[node] = (int)que.size();
            par[node] = parent;
        }
        return true;
    }, index_f, -1);
    std::vector<std::vector<std::vector<int>>> vvc;
    std::vector<int> vis(adjl.size(), -1);
    for (int u : que) {
        if (par[u] == -1) vvc.push_back({});
        for (const auto& e : adjl[u]) {
            int v = index_f(e);
            if (par[v] == u || ord[v] < ord[u]) continue;
            std::vector<int> c{ u, v };
            vis[u] = 1;
            do {
                vis[v] = 1;
                v = par[v];
                c.push_back(v);
            } while (vis[v] == -1);
            vvc.back().push_back(c);
        }
    }
    return vvc;
}
template<typename E = int>
std::vector<std::vector<std::vector<int>>> chain_decomposition(const std::vector<std::vector<int>>& adjl) {
    return chain_decomposition(adjl, [](int i){ return i; });
}

/**
 * Calculates all the cut vertices and edges.
 * (whose removal would increase number of components in the graph)
 *
 * Complexity: O(m)
 *
 * @param adjl - adjacency list of `n` nodes and `m` edges
 * @param index_f - a functor that extracts the outbound node index from the edge
 */
template<typename E, typename FI>
std::pair<std::vector<int>, std::vector<std::pair<int, int>>> cut_vertices_and_edges(const std::vector<std::vector<E>>& adjl, FI index_f) {
    // unordered_set could be avoided by having chain_decomposition operate on edge indices
    auto vvc = chain_decomposition(adjl, index_f);
    std::unordered_set<std::pair<int, int>> cut_edges;
    for (int u = 0; u < (int)adjl.size(); u++) {
        for (const auto& e : adjl[u]) {
            int v = index_f(e);
            if (u < v) cut_edges.insert({ u, v });
        }
    }
    for (const auto& vc : vvc) {
        for (const auto& c : vc) {
            for (int i = 1; i < (int)c.size(); i++) {
                int u = c[i - 1], v = c[i];
                if (u > v) std::swap(u, v);
                cut_edges.erase({ u, v });
            }
        }
    }
    std::unordered_set<int> cut_vertices;
    for (const auto& vc : vvc) {
        for (int i = 1; i < (int)vc.size(); i++) {
            int u = vc[i].front(), v = vc[i].back();
            if (u == v) cut_vertices.insert(u);
        }
    }
    for (const auto& e : cut_edges) {
        int u = e.first, v = e.second;
        if (adjl[u].size() > 1) cut_vertices.insert(u);
        if (adjl[v].size() > 1) cut_vertices.insert(v);
    }
    std::vector<std::pair<int, int>> ve(cut_edges.begin(), cut_edges.end());
    std::vector<int> vv(cut_vertices.begin(), cut_vertices.end());
    return{ vv, ve };
}
template<typename E = int>
std::pair<std::vector<int>, std::vector<std::pair<int, int>>> cut_vertices_and_edges(const std::vector<std::vector<int>>& adjl) {
    return cut_vertices_and_edges(adjl, [](int i){ return i; });
}

} // graph
} // altruct
