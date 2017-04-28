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
    
    graph() : {}
    graph(int n) : adjl(n) {}
    graph(const std::vector<std::vector<E>>& adjl) : adjl(adjl) {}

    bool operator < (const graph<E>& rhs) const { return adjl < rhs.adjl; }
    bool operator == (const graph<E>& rhs) const { return adjl == rhs.adjl; }

    int size() const { return (int)adjl.size(); }
    int add_node() { adjl.push_back({}); return size() - 1; }
    void add_edge(int u, const E& e) { adjl[u].push_back(e); }
    std::vector<E>& operator[] (int u) { return adjl[u]; }
    const std::vector<E>& operator[] (int u) const { return adjl[u]; }
};

} // graph
} // altruct

namespace std {
template<>
struct hash<altruct::graph::full_edge> {
    size_t operator()(const altruct::graph::full_edge& e) const { return e.u * 0x9e3779b9 + e.v; }
};
}
