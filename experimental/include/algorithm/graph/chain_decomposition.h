#pragma once

#include "graph.h"

#include <vector>
#include <unordered_set>
#include <algorithm>
#include "algorithm/hash/std_tuple_hash.h"

namespace altruct {
namespace graph {

/**
 * Calculates the chain decomposition of an undirected graph.
 *
 * Complexity: O(m)
 *
 * @param g - an undirected graph with `n` nodes and `m` edges
 * @return - a list of chains for each component of the graph 'g';
 *           each chain is represented as a list of vertex indices
 */
template<typename E>
std::vector<std::vector<std::vector<int>>> chain_decomposition(const graph<E>& g) {
    std::vector<int> que; 
    std::vector<int> ord(g.size(), -1);
    std::vector<int> par(g.size(), -1);
    iterative_dfs(g, [&](int root, int parent, int node, int depth) {
        if (node != -1) {
            que.push_back(node);
            ord[node] = (int)que.size();
            par[node] = parent;
        }
        return true;
    });
    std::vector<std::vector<std::vector<int>>> vvc;
    std::vector<int> vis(g.size(), -1);
    for (int u : que) {
        if (par[u] == -1) vvc.push_back({});
        for (const auto& e : g[u]) {
            int v = e.v;
            if (par[v] == u || ord[v] < ord[u]) continue;
            std::vector<int> c{ u, v };
            vis[u] = 1;
            while (vis[v] == -1) {
                vis[v] = 1;
                v = par[v];
                c.push_back(v);
            };
            vvc.back().push_back(c);
        }
    }
    return vvc;
}

/**
 * Calculates all the cut edges (bridges) of an undirected graph.
 *
 * Cut-edge removal would increase the number of components in the graph.
 *
 * Complexity: O(m)
 *
 * @param g - an undirected graph with `n` nodes and `m` edges
 * @param d - a `chain_decomposition` of the graph `g`
 * @return - a list of cut edges in the graph `g`
 */
template<typename E>
std::vector<full_edge> cut_edges(const graph<E>& g, const std::vector<std::vector<std::vector<int>>>& d) {
    // unordered_set could be avoided by having chain_decomposition operate on edge indices
    std::unordered_set<full_edge> ce;
    for (int u = 0; u < g.size(); u++) {
        for (const auto& e : g[u]) {
            if (u < e.v) ce.insert({ u, e.v }); // edge is a bridge ...
        }
    }
    for (const auto& vc : d) {
        for (const auto& c : vc) {
            for (int i = 1; i < (int)c.size(); i++) {
                int u = c[i - 1], v = c[i];
                if (u > v) std::swap(u, v);
                ce.erase({ u, v }); // ... if it is not in any chain
            }
        }
    }
    return std::vector<full_edge>(ce.begin(), ce.end());
}

/**
 * Calculates all the cut vertices (articulation points) of an undirected graph.
 *
 * Cut-vertex removal would increase the number of components in the graph.
 *
 * Complexity: O(m)
 *
 * @param g - an undirected graph with `n` nodes and `m` edges
 * @param d - a `chain_decomposition` of the graph `g`
 * @return - a list of cut vertices in the graph `g`
 */
template<typename E>
std::vector<int> cut_vertices(const graph<E>& g, const std::vector<std::vector<std::vector<int>>>& d) {
    std::vector<int> is_cut(g.size(), 0);
    for (int u = 0; u < g.size(); u++) {
        if (g[u].size() >= 2) is_cut[u] = 1; // vertex is an articulation point if its degree is at least 2 ...
    }
    for (const auto& vc : d) {
        for (const auto& c : vc) {
            for (int u : c) is_cut[u] = 0; // ... and is not in any biconnected component ...
        }
    }
    for (const auto& vc : d) {
        for (int i = 1; i < (int)vc.size(); i++) {
            int u = vc[i].front(), v = vc[i].back();
            if (u == v) is_cut[u] = 1; // ... or it is a first vertex of a non-first cycle (per each component)
        }
    }
    std::vector<int> cv;
    for (int u = 0; u < g.size(); u++) {
        if (is_cut[u]) cv.push_back(u);
    }
    return cv;
}

/**
 * Calculates all the biconnected components of an undirected graph.
 *
 * Note, bridge components are not considered to be biconnected and hence not returned.
 *
 * Complexity: O(m)
 *
 * @param g - an undirected graph with `n` nodes and `m` edges
 * @param d - a `chain_decomposition` of the graph `g`
 * @return - a list of biconnected components in the graph `g`
 */
template<typename E>
std::vector<std::vector<int>> biconnected_components(const graph<E>& g, const std::vector<std::vector<std::vector<int>>>& d) {
    std::vector<int> seen(g.size(), 0);
    std::vector<std::vector<int>> vbc;
    for (const auto& vc : d) {
        for (const auto& c : vc) {
            int u = c.front(), v = c.back();
            if (u == v) vbc.push_back({}), seen[u] = 0; // start of a new biconnected component
            for (int u : c) {
                if (!seen[u]) vbc.back().push_back(u), seen[u] = 1;
            }
        }
    }
    return vbc;
}

} // graph
} // altruct
