#pragma once

#include "graph.h"
#include "algorithm/graph/chain_decomposition.h"
#include "structure/math/polynom.h"

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
 * Calculates the chromatic polynomial of an undirected graph.
 *
 * Complexity: O(phi^n)
 *
 * @param g - an undirected graph with `n` nodes and `m` edges
 * @return - chromatic polynomial of the graph 'g'
 */
template<typename I, typename E>
altruct::math::polynom<I> chromatic_polynomial(const graph<E>& g, I id = I(1)) {
    typedef altruct::math::polynom<I> poly;
    int n = g.size();
    poly k = { altruct::math::zeroOf(id), id };
    poly k1 = { -id, id };
    // Delete vertex that is connected to every other vertex (if it exists)
    for (int u = 0; u < n; u++) {
        if (g[u].size() < n - 1) continue;
        std::vector<int> b(n);
        int cnt = 0;
        for (const E& e : g[u]) {
            if (!b[e.v]) b[e.v] = 1, cnt++;
        }
        if (cnt < n - 1) continue;
        auto gd = g; gd.delete_node(u);
        return k * chromatic_polynomial<I>(gd)(k1);
    }
    // Handle each biconnected component independently
    auto d = chain_decomposition(g);
    int C = (int)d.size(), E2 = 0, CE = 0;
    for (int u = 0; u < n; u++) {
        E2 += (int)g[u].size(); // this counts each edge twice
    }
    std::vector<graph<edge>> vg;
    std::vector<int> idx(n, -1);
    for (const auto& vc : d) {
        for (const auto& c : vc) {
            CE += (int)c.size() - 1;
            int u = c.front(), v = c.back();
            // cycle indicates start of a new biconnected component
            if (u == v) vg.push_back({}), idx[u] = -1;
            for (int u : c) {
                if (idx[u] == -1) idx[u] = vg.back().add_node();
            }
            for (int i = 1; i < (int)c.size(); i++) {
                vg.back().add_edge(idx[c[i - 1]], idx[c[i]]);
                vg.back().add_edge(idx[c[i]], idx[c[i - 1]]);
            }
        }
    }
    if (vg.size() != 1) {
        poly p = powT(k, C) * powT(k1, E2 / 2 - CE);
        for (const auto& gg : vg) {
            p *= chromatic_polynomial<I>(gg) / k;
        }
        return p;
    } else {
        // Deletion–contraction recursion
        // select arbitrary edge (TODO: use better heuristics)
        int u = -1; for (u = 0; u < n && g[u].empty(); u++);
        int v = g[u].back().v;
        auto gd = g; gd.delete_edge(u, v); gd.delete_edge(v, u);
        auto gc = g; gc.contract(u, v);
        return chromatic_polynomial<I>(gd) -chromatic_polynomial<I>(gc);
    }
}

} // graph
} // altruct
