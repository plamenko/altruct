#pragma once

#include <vector>
#include <unordered_set>
#include <algorithm>

namespace altruct {
namespace graph {

/**
 * A base class for edge types.
 * Only the destination vertex is specified.
 */
struct edge {
    int v; // the destination vertex
    edge(int v) : v(v) {}
    bool operator < (const edge& rhs) const { return v < rhs.v; }
    bool operator == (const edge& rhs) const { return v == rhs.v; }
};

/**
 * An edge representation where both of its vertices are specified.
 */
struct full_edge : public edge {
    int u; // the source vertex
    full_edge(int u, int v) : edge(v), u(u) {}
    bool operator < (const full_edge& rhs) const { return (u != rhs.u) ? (u < rhs.u) : (v < rhs.v); }
    bool operator == (const full_edge& rhs) const { return (u == rhs.u) && (v == rhs.v); }
};

/**
 * A weighted edge.
 */
template<typename W>
struct weighted_edge : public edge {
    W w; // the weight
    weighted_edge(int v, W w) : edge(v), w(w) {}
};

/**
 * A graph represented by its adjacency list.
 *
 * Note, the graph is considered to be directed. Therefore, it is the caller's
 * responsibility to make sure that for undirected graphs edges always come in
 * pairs. I.e. if the graph is to be treated as undirected, for an edge {u, v},
 * there should be a corresponding edge {v, u}, with the same accompanying data
 * such as weight and similar.
 *
 * @param E - edge type
 */
template<typename E>
class graph {
public:
    std::vector<std::vector<E>> adjl;

    graph() {}
    graph(int n) : adjl(n) {}
    graph(const std::vector<std::vector<E>>& adjl) : adjl(adjl) {}

    bool operator < (const graph<E>& rhs) const { return adjl < rhs.adjl; }
    bool operator == (const graph<E>& rhs) const { return adjl == rhs.adjl; }

    int size() const { return (int)adjl.size(); }
    int num_edges() const { int e = 0; for (const auto& l : adjl) e += (int)l.size(); return e; }
    std::vector<E>& operator[] (int u) { return adjl[u]; }
    const std::vector<E>& operator[] (int u) const { return adjl[u]; }

    int add_node() { adjl.push_back({}); return size() - 1; }
    void add_edge(int u, const E& e) { adjl[u].push_back(e); }
    void add_edge2(int u, const E& e) { adjl[u].push_back(e); adjl[e.v].push_back(e); adjl[e.v].back().v = u; }

    void delete_edge(int u, int v) {
        for (int i = 0; i < (int)adjl[u].size();) {
            if (adjl[u][i].v == v) {
                std::swap(adjl[u][i], adjl[u].back());
                adjl[u].pop_back();
            } else {
                i++;
            }
        }
    }
    void delete_node(int u) {
        int v = size() - 1;
        std::swap(adjl[u], adjl[v]);
        adjl.pop_back();
        for (auto& l : adjl) {
            for (int i = 0; i < (int)l.size();) {
                if (l[i].v == u) {
                    std::swap(l[i], l.back());
                    l.pop_back();
                } else {
                    if (l[i].v == v) l[i].v = u;
                    i++;
                }
            }
        }
    }
    void contract(int u, int v) {
        delete_edge(u, v);
        delete_edge(v, u);
        adjl[u].insert(adjl[u].end(), adjl[v].begin(), adjl[v].end());
        for (auto& l : adjl) for (auto& e : l) if (e.v == v) e.v = u;
        delete_node(v);
        for (int u = 0; u < size(); u++) deduplicate_edges(u);
    }
    void deduplicate_edges(int u) {
        sort(adjl[u].begin(), adjl[u].end());
        auto it = unique(adjl[u].begin(), adjl[u].end());
        adjl[u].erase(it, adjl[u].end());
    }
};

} // graph
} // altruct

namespace std {
template<>
struct hash<altruct::graph::full_edge> {
    size_t operator()(const altruct::graph::full_edge& e) const { return e.u * 0x9e3779b9 + e.v; }
};
}
