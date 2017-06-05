#pragma once

#include "structure/graph/graph.h"

#include "algorithm/graph/iterative_dfs.h"
#include "algorithm/graph/lowest_common_ancestor.h"

#include <vector>

namespace altruct {
namespace graph {

/**
 * Heavy-light Decomposition of a tree.
 *
 * Decomposes the tree to a set of chains in such a way that
 * every possible path goes through at most `O(log n)` chains.
 *
 * Note that each edge can be uniquely represented by its lower node.
 *
 * Chains are laid out continuously on a number line.
 * Each subtree is also continuous on a number line.
 *
 * Space complexity: `O(n)`
 * Time complexities:
 *   build: `O(n)`
 *   query: `O(1)`
 */
class heavy_light_decomposition {
	std::vector<int> parents;
	std::vector<int> positions;
	std::vector<int> nodes;
	std::vector<int> sizes;

public:
    template<typename E>
	heavy_light_decomposition(graph<E>& g) :
		parents(g.size(), -1),
		positions(g.size(), -1),
		nodes(g.size(), -1),
		sizes(g.size(), 1) {
		iterative_dfs(g, [&](int root, int parent, int node, int depth) {
			if (node == -1) {
				for (const E& e : g[parent]) {
					if (e.v != parents[parent]) sizes[parent] += sizes[e.v];
				}
			} else {
				parents[node] = parent;
			}
			return true;
		});
		for (int u = 0; u < g.size(); u++) {
			if (g[u].size() < 2) continue;
			int heavy_i = 0, heavy_size = 0;
			for (int i = 0; i < g[u].size(); i++) {
				int v = g[u][i].v;
				if (v != parents[u] && heavy_size < sizes[v]) {
					heavy_i = i, heavy_size = sizes[v];
				}
			}
			std::swap(g[u][0].v, g[u][heavy_i].v);
		}
		int pos = 0, par = -2;
		iterative_dfs(g, [&](int root, int parent, int node, int depth) {
			if (node != -1) {
				if (par == -2) par = parent;
				parents[node] = par;
				positions[node] = pos;
				nodes[pos++] = node;
			} else {
				par = -2;
			}
			return true;
		});
	}

	// gives the parent node of the chain that contains the node `u`.
	// this is the parent of the last node in the chain, and is therefore
	// not considered to be part of the chain.
	int parent(int u) {
		return parents[u];
	}

	// gives the position of the node `u` in the linearized tree.
	int position(int u) {
		return positions[u];
	}

	// gives the node at the position `pos` in the linearized tree.
	int node(int pos) {
		return nodes[pos];
	}

	// gives the size of the subtree rooted at node `u`
	int subtree_size(int u) {
		return sizes[u];
	}
};

class heavy_light_decomposition_ex {
public:
	lowest_common_ancestor lca;
	heavy_light_decomposition hld;

    template<typename E>
	heavy_light_decomposition_ex(graph<E>& g) :
		lca(g),
		hld(g) {
	}

	// `k-th` parent (ancestor) of `u`;
	// -1 if `u` has less than `k` ancestors
	int parent(int u, int k = 1) {
		if (k < 0 || k > lca.depth(u)) return -1;
		int d = lca.depth(u) - k;
		while (true) {
			int p = hld.parent(u);
			if (lca.depth(p) < d) {
				return hld.node(hld.position(u) - (lca.depth(u) - d));
			}
			u = p;
		}
	}

	// Calls `visitor(chain_begin, chain_end, path_first, path_last, path_length, direction)`
	// for each linear segment `[chain_begin, chain_end)`, and the corresponding
	// path segment `[path_first, path_last]` on the path from `u` to `v`.
	// `chain_begin` and `chain_end` are indices in the HLD-linearized tree.
	// `path_first` and `path_last` are distances relative to the beginning of the
	// path (i.e. distances from `u`).
	// Note that the bounds of a path segment are both inclusive, and also that
	// `path_first` might be bigger than `path_last` if a path segment goes in the
	// direction opposite of the direction of its corresponding chain segment.
	template<typename F>
	int walk(int u, int v, F visitor) {
		int a = lca.ancestor(u, v);
		int uv_dist = lca.distance(u, v);
		// going up from `u` to `a`, note that this is in reverse direction of the chain segment
		for (int w = u; w != a;) {
			int p = hld.parent(w);
			if (lca.depth(p) < lca.depth(a)) p = a;
			int len = lca.depth(w) - lca.depth(p);
			int end_pos = hld.position(w) + 1;
			int uw_dist = lca.depth(u) - lca.depth(w);
			visitor(end_pos - len, end_pos, uw_dist + len - 1, uw_dist, uv_dist, true);
			w = p;
		}
		// going up from `v` to `a`
		for (int w = v; w != a;) {
			int p = hld.parent(w);
			if (lca.depth(p) < lca.depth(a)) p = a;
			int len = lca.depth(w) - lca.depth(p);
			int end_pos = hld.position(w) + 1;
			int uw_dist = uv_dist - (lca.depth(v) - lca.depth(w));
			visitor(end_pos - len, end_pos, uw_dist - len, uw_dist - 1, uv_dist, false);
			w = p;
		}
		return uv_dist;
	}
};

}
}
