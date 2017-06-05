#pragma once

#include "altruct/structure/graph/graph.h"
#include "altruct/algorithm/graph/chain_decomposition.h"
#include "altruct/structure/math/polynom.h"

#include <vector>
#include <algorithm>

namespace altruct {
namespace graph {

/** Chromatic polynomial of a tree graph. */
template<typename I>
altruct::math::polynom<I> chromatic_polynomial_T(int n, I id = I(1)) {
    typedef altruct::math::polynom<I> poly;
    return poly{ altruct::math::zeroOf(id), id } *powT(poly{ -id, id }, n - 1);
}

/** Chromatic polynomial of a cycle graph. */
template<typename I>
altruct::math::polynom<I> chromatic_polynomial_C(int n, I id = I(1)) {
    typedef altruct::math::polynom<I> poly;
    return powT(poly{ -id, id }, n) + poly{ -id, id } *((n % 2 == 0) ? +1 : -1);
}

/** Chromatic polynomial of a complete graph. */
template<typename I>
altruct::math::polynom<I> chromatic_polynomial_K(int n, I id = I(1)) {
    typedef altruct::math::polynom<I> poly;
    poly p{ id };
    for (int i = 0; i < n; i++) p *= poly{ -id * i, id };
    return p;
}

/**
 * Calculates the chromatic polynomial of an undirected simple graph.
 *
 * Graph must be:
 *   undirected: for each edge (u, v), there must be a corresponding edge (v, u).
 *   simple: there should be no self loops and each edge (u, v) should occur at most once.
 *
 * Complexity: O(phi^(n+m))
 *
 * @param g - an undirected simple graph with `n` nodes and `m` edges
 * @return - chromatic polynomial of the graph 'g'
 */
template<typename I, typename E>
altruct::math::polynom<I> chromatic_polynomial(const graph<E>& g, I id = I(1)) {
    typedef altruct::math::polynom<I> poly;
    int n = g.size();
    poly k = { altruct::math::zeroOf(id), id };
    poly k1 = { -id, id };
    // Delete vertex that is connected to every other vertex (if exists)
    for (int u = 0; u < n; u++) {
        if (g[u].size() < n - 1) continue;
        auto gd = g; gd.delete_node(u);
        return k * chromatic_polynomial<I>(gd)(k1);
    }
    // Handle each biconnected component independently
    auto d = chain_decomposition(g);
    int C = (int)d.size(), AE = g.num_edges() / 2, CE = 0;
    std::vector<graph<edge>> vg;
    std::vector<int> idx(n, -1);
    for (const auto& component : d) {
        for (const auto& biconnected : component) {
            vg.push_back({}); // a graph for each  biconnected component
            idx[biconnected[0][0]] = -1;
            for (const auto& chain : biconnected) {
                CE += (int)chain.size() - 1;
                for (int u : chain) {
                    if (idx[u] == -1) idx[u] = vg.back().add_node();
                }
                for (int i = 1; i < (int)chain.size(); i++) {
                    int u = idx[chain[i - 1]], v = idx[chain[i]];
                    vg.back().add_edge(u, v);
                    vg.back().add_edge(v, u);
                }
            }
        }
    }
    if (vg.size() != 1 || C != 1 || AE != CE) {
        poly p = powT(k, C) * powT(k1, AE - CE);
        for (const auto& gg : vg) {
            p *= chromatic_polynomial<I>(gg) / k;
        }
        return p;
    } else {
        // if reached this point, the graph consists of a single biconnected
        // component and no vertex is connected to all other vertices.
        // solve recursively
        if (AE < n * n / 4) {
            // sparse: deletion-contraction
            // delete edge incident to the smallest degree vertex
            int u = 0; for (int i = 0; i < n; i++) if (g[u].size() > g[i].size()) u = i;
            int v = -1; for (const E& e : g[u]) if (v < 0 || g[v].size() > g[e.v].size()) v = e.v;
            auto gd = g; gd.delete_edge(u, v); gd.delete_edge(v, u);
            auto gc = g; gc.contract(u, v);
            return chromatic_polynomial<I>(gd) - chromatic_polynomial<I>(gc);
        } else {
            // dense: addition-contraction
            // add edge incident to the highest degree vertex
            int u = 0; for (int i = 0; i < n; i++) if (g[u].size() < g[i].size()) u = i;
            std::vector<int> used(n); used[u] = 1; for (const E& e : g[u]) used[e.v] = 1;
            int v = -1; for (int i = 0; i < n; i++) if ((v < 0 || g[v].size() < g[i].size()) && !used[i]) v = i;
            auto ga = g; ga.add_edge(u, v); ga.add_edge(v, u);
            auto gc = g; gc.contract(u, v);
            return chromatic_polynomial<I>(ga) + chromatic_polynomial<I>(gc);
        }
    }
}

} // graph
} // altruct
