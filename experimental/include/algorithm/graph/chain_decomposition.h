#pragma once

#include "graph.h"

#include <vector>
#include <unordered_set>
#include <algorithm>
#include "algorithm/hash/std_tuple_hash.h"

namespace altruct {
namespace graph {

/**
 * Represents a chain decomposition of a graph.
 *
 * This four-level nesting consists of a list of connected components.
 * Each connected component consists of a list of bi-connected components.
 * Each bi-connected component consists of a list of chains, the first
 * chain being a cycle and the rest of them simple paths.
 * Each chain consists of a list of vertices in the order they appear in it.
 * 
 * Essentially: [component_id][biconnected_component_id][chain_id][vertex_id]
 */
typedef std::vector<std::vector<std::vector<std::vector<int>>>> chain_decomposition_t;

/**
 * Calculates the chain decomposition of an undirected graph.
 *
 * Complexity: O(m)
 *
 * @param g - an undirected graph with `n` nodes and `m` edges
 * @return - see the structure of `chain_decomposition_t` type
 */
template<typename E>
chain_decomposition_t chain_decomposition(const graph<E>& g) {
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
    chain_decomposition_t res;
    std::vector<int> bid(g.size(), -1); // biconnected component index
    for (int u : que) {
        if (par[u] == -1) res.push_back({});
        auto& component = res.back();
        for (const auto& e : g[u]) {
            int v = e.v;
            if (par[v] == u || ord[v] <= ord[u]) continue;
            std::vector<int> chain{ u, v };
            bid[u] = (int)component.size();
            while (bid[v] == -1) {
                v = par[v];
                chain.push_back(v);
            };
            for (int w : chain) {
                bid[w] = bid[v];
            }
            if (bid[v] == (int)component.size()) {
                component.push_back({});
            }
            component[bid[v]].push_back(chain);
        }
    }
    return res;
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
std::vector<full_edge> cut_edges(const graph<E>& g, const chain_decomposition_t& d) {
    // unordered_set could be avoided by having chain_decomposition operate on edge indices
    std::unordered_set<full_edge> cute;
    for (int u = 0; u < g.size(); u++) {
        for (const auto& e : g[u]) {
            if (u < e.v) cute.insert({ u, e.v }); // edge is a bridge ...
        }
    }
    for (const auto& component : d) {
        for (const auto& biconnected : component) {
            for (const auto& chain : biconnected) {
                for (int i = 1; i < (int)chain.size(); i++) {
                    int u = chain[i - 1], v = chain[i];
                    if (u > v) std::swap(u, v);
                    cute.erase({ u, v }); // ... if it is not in any chain
                }
            }
        }
    }
    return std::vector<full_edge>(cute.begin(), cute.end());
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
std::vector<int> cut_vertices(const graph<E>& g, const chain_decomposition_t& d) {
    std::vector<int> is_cut(g.size(), 0);
    for (int u = 0; u < g.size(); u++) {
        // vertex is an articulation point if its degree is at least 2 ...
        if (g[u].size() >= 2) is_cut[u] = 1;
    }
    for (const auto& component : d) {
        for (const auto& biconnected : component) {
            for (const auto& chain : biconnected) {
                // ... and is not in any biconnected component ...
                for (int u : chain) is_cut[u] = 0;
            }
        }
    }
    for (const auto& component : d) {
        for (int i = 1; i < (int)component.size(); i++) {
            // ... or it is a first vertex of a non-first biconnected-component
            is_cut[component[i][0][0]] = 1;
        }
    }
    std::vector<int> cutv;
    for (int u = 0; u < g.size(); u++) {
        if (is_cut[u]) cutv.push_back(u);
    }
    return cutv;
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
std::vector<std::vector<int>> biconnected_components(const graph<E>& g, const chain_decomposition_t& d) {
    std::vector<int> seen(g.size(), 0);
    std::vector<std::vector<int>> vbc;
    for (const auto& component : d) {
        for (const auto& biconnected : component) {
            vbc.push_back({}); // start of a new biconnected component
            seen[biconnected[0][0]] = 0;
            for (const auto& chain : biconnected) {
                for (int u : chain) {
                    if (!seen[u]) vbc.back().push_back(u), seen[u] = 1;
                }
            }
        }
    }
    return vbc;
}

} // graph
} // altruct
