#pragma once

#include "altruct/structure/graph/graph.h"
#include "altruct/algorithm/graph/tarjan_scc.h"

#include <vector>
#include <algorithm>

namespace altruct {
namespace graph {

/**
 * Solves the 2-SAT problem.
 *
 * Complexity: O(n)
 *
 * @param sol - output solution consisting of an assignment to each variable
 *     such that the equation is satisfied (if such an assignment exists).
 *     `-1` means unassigned, `0` means false, `1` means true. The initial
 *     values may all be unassigned (i.e. set to `-1`), or some of them can
 *     be assigned to `0` or `1` in which case these will be respected.
 * @param clauses - list of clauses where each clause is a disjunction
 *     of two variables. Positive literal for variable `i` is encoded
 *     as `2*i+1`, whereas its negation is encoded as `2*i+0`.
 *     `-1` is treated as false. I.e. clause `{a,-1}` is equivalent to
 *     `(a or false)` which is equivalent to just `a`. This is useful
 *     when representing the single literal clauses.
 * @return - true iff the equation is satisfiable
 */
template<typename T = void>
bool sat2(std::vector<int>& sol, const std::vector<std::pair<int, int>>& clauses) {
    int n = 0; // number of variables
    for (const auto& t : clauses) {
        n = std::max(n, t.first / 2 + 1);
        n = std::max(n, t.second / 2 + 1);
    }
    if (sol.size() < n) sol.resize(n, -1);
    graph<edge> g(n * 2); // the implication graph
    for (int i = 0; i < n; i++) {
        if (sol[i] != -1) {
            int u = i * 2 + sol[i];
            g.add_edge(u ^ 1, u);
        }
    }
    for (const auto& t : clauses) {
        int u = t.first, v = t.second;
        if (u == -1 && v == -1) return false;
        if (u == -1) u = v;
        if (v == -1) v = u;
        g.add_edge(u ^ 1, v); // !u => v
        g.add_edge(v ^ 1, u); // !v => u
    }
    for (auto scc = tarjan_scc(g); !scc.empty(); scc.pop_back()) {
        for (int u : scc.back()) {
            if (sol[u / 2] == -1) sol[u / 2] = u % 2;
        }
    }
    for (const auto& t : clauses) {
        int u = t.first, v = t.second;
        if (u == -1 && v == -1) return false;
        if (u == -1) u = v;
        if (v == -1) v = u;
        if (sol[u / 2] != (u % 2) && sol[v / 2] != (v % 2)) return false;
    }
    return true;
}

template<typename T = void>
std::vector<int> sat2(const std::vector<std::pair<int, int>>& clauses) {
    std::vector<int> sol;
    if (!sat2(sol, clauses)) sol.clear();
    return sol;
}

} // graph
} // altruct
