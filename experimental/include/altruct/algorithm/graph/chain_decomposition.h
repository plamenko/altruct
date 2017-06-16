#pragma once

#include "altruct/structure/graph/graph.h"

#include <vector>
#include <unordered_set>
#include <algorithm>
#include "altruct/algorithm/hash/std_tuple_hash.h"

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
    // 0 - not an articulation point
    // 1 - articulation point that's part of some biconnected block
    // 2 - articulation point that's not part of any biconnected block
    std::vector<int> is_cut(g.size(), 2);
    // vertex is an articulation point if:
    for (const auto& component : d) {
        for (const auto& biconnected : component) {
            for (const auto& chain : biconnected) {
                // is not in any biconnected component ...
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
    for (int u = 0; u < g.size(); u++) {
        // ... or it is connected to some tree
        for (const auto& e : g[u]) {
            if (!is_cut[u] && is_cut[e.v] == 2) is_cut[u] = 1;
        }
    }
    for (int u = 0; u < g.size(); u++) {
        // .. and its degree is at least 2
        if (g[u].size() < 2) is_cut[u] = 0;
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

/**
 * Builds a Block-Cut Tree of an undirected graph and the corresponding node map.
 *
 * Note: The first `vb.size()` nodes of the resulting tree correspond 1-to-1 to the
 * biconnected components. The following `va.size()` nodes correspond 1-to-1 to the
 * articulation points. The remaining nodes correspond to leaf and isolated nodes.
 * I.e. nodes that were not in any biconnected component nor an articulation point.
 * Nodes in the resulting tree are connected if and only if the corresponding nodes
 * were connected in the original graph and were not in the same bi-component.
 *
 * Complexity: O(m)
 *
 * @param g - an undirected graph with `n` nodes and `m` edges
 * @param vb - list of biconnected components (each component is a list of vertices)
 * @param va - list of articulation points (cut vertices)
 * @return - block-cut tree and a map from original nodes to the block-cut tree nodes
 */
template<typename E>
std::pair<graph<E>, std::vector<int>> block_cut_tree(const graph<E>& g, const std::vector<std::vector<int>>& vb, const std::vector<int>& va) {
    graph<E> t;
    vector<int> idx(g.size(), -1);
    vector<int> non(g.size(), 1);
    // nodes for blocks (biconnected components)
    for (const auto& b : vb) {
        int i = t.add_node();
        for (auto u : b) idx[u] = i;
    }
    // nodes for articulation points
    for (int u : va) {
        idx[u] = t.add_node();
    }
    // nodes for leafs & isolated nodes
    for (int u = 0; u < g.size(); u++) {
        if (g[u].size() <= 1) idx[u] = t.add_node();
    }
    // edges for block - art.pt.
    for (int i = 0; i < (int)vb.size(); i++) {
        for (auto u : vb[i]) {
            non[u] = 0;
            if (idx[u] != i) t.add_edge2(idx[u], i);
        }
    }
    // edges for art.pt. - non-block
    for (int u = 0; u < g.size(); u++) {
        for (const auto& e : g[u]) {
            if (non[u] || non[e.v]) {
                t.add_edge(idx[u], idx[e.v]);
            }
        }
    }
    return{ t, idx };
}

/**
 * Encapsulates various biconnectivity information of an undirected graph.
 *
 * Complexity: O(m)
 *
 * @param g - an undirected graph with `n` nodes and `m` edges
 */
template<typename E = edge>
struct biconnectivity {
    typedef std::vector<int> vertices_t;
    chain_decomposition_t decomposition; // chain_decomposition
    std::vector<full_edge> cut_edges;    // cut edges (bridges)
    vertices_t cut_vertices;             // cut vertices (articulation points)
    std::vector<vertices_t> components;  // biconnected components
    graph<E> block_cut_tree;             // block-cut tree
    vertices_t bc_tree_idx;              // mapping from original nodes to the block-cut tree nodes

    biconnectivity(const graph<E>& g) {
        decomposition = ::chain_decomposition(g);
        cut_edges = ::cut_edges(g, decomposition);
        cut_vertices = ::cut_vertices(g, decomposition);
        components = ::biconnected_components(g, decomposition);
        auto ti = ::block_cut_tree<>(g, components, cut_vertices);
        block_cut_tree = ti.first;
        bc_tree_idx = ti.second;
    }
};

} // graph
} // altruct
