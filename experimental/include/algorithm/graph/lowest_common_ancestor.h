#pragma once

#include "graph.h"

#include "algorithm/graph/iterative_dfs.h"
#include "structure/container/segment_tree.h"

#include <algorithm>
#include <vector>

namespace altruct {
namespace graph {

/**
 * Structure for fast LCA queries
 *
 * Space complexity: `O(n)`
 * Time complexities:
 *   build: `O(n)`
 *   query: `O(log n)`
 */
class lowest_common_ancestor {
	typedef std::pair<int, int> pii;
	static pii min_f(pii e1, pii e2) { return std::min(e1, e2); }

	std::vector<int> indices; // `indices[u]` is the index of the last occurrence of node `u` in `levels`
	altruct::container::segment_tree<pii> levels; // segment tree for range-minimum-queries on the (level, node)

public:
	// @param g - an undirected tree (or a forest, but no cycles);
	//            if there is an edge (u,v), there should also be an edge (v,u)
	template<typename E>
    lowest_common_ancestor(const graph<E>& g) :
		indices(g.size(), -1),
		levels(g.size() * 2 + 1, min_f, { std::numeric_limits<int>::max(), -1 }) {
		int s = 0;
		iterative_dfs(g, [&](int root, int parent, int node, int depth) {
			if (parent != -1) indices[parent] = s;
			levels[s++] = { depth - 1, parent };
			return true;
		});
		levels[s++] = { -1, -1 };
		levels.rebuild();
	}

	// depth of the node `u`
	int depth(int u) {
		if (u < 0 || indices[u] < 0) return -1;
		return levels.get(indices[u]).first;
	}

	// lowest common ancestor of `u` and `v`, -1 if not connected
	int ancestor(int u, int v) {
		if (u < 0 || v < 0 || indices[u] < 0 || indices[v] < 0) return -1;
		if (indices[u] > indices[v]) std::swap(u, v);
		return levels.get(indices[u], indices[v] + 1).second;
	}

	// distance between `u` and `v`, -1 if not connected
	int distance(int u, int v) {
		int a = ancestor(u, v); if (a < 0) return -1;
		return (depth(u) - depth(a)) + (depth(v) - depth(a));
	}
};

}
}
