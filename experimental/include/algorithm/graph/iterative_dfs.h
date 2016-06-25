#pragma once

#include <vector>

namespace altruct {
namespace graph {

/**
 * Performs DFS and invokes `visitor` on each step.
 *
 * If `source` is specified only the component rooted at `source` will be traversed.
 * Otherwise, each component will be traversed (the graph may be disconnected).
 *
 * @param adjl - adjacency list
 * @param visitor - a function `bool visitor(root, parent, node, depth)`
 *                  @param root - the root node of the current component
 *                  @param parent - the parent of the current node,
 *                                -1 if the current node is the root
 *                                (visited in euler-tour-order)
 *                  @param node - the current node; an artificial -1 is
 *                                reported as the last node of each parent
 *                                (visited in pre-order)
 *                  @param depth - the depth of the current node
 *                  @return - true if the node should be entered
 * @param source - the source node, or -1 if not specified
 */
template<typename F>
void iterative_dfs(const std::vector<std::vector<int>>& adjl, F visitor, int source = -1) {
	std::vector<int> visited(adjl.size());
	std::vector<std::pair<int, int>> stk;
	// this outer loop allows us to handle a disconnected graph
	for (int o = 0; o < (int)adjl.size(); o++) {
		if (source != -1) o = source;
		// handle the component rooted at `o`
		if (visited[o]) continue;
		visited[o] = 1;
		stk.push_back({ o, 0 });
		visitor(o, -1, o, 0);
		while (!stk.empty()) {
			int d = (int)stk.size();
			int u = stk.back().first;
			int i = stk.back().second;
			if (i < (int)adjl[u].size()) {
				int v = adjl[u][i];
				stk.back().second++;
				if (visited[v]) continue;
				if (visitor(o, u, v, d)) {
					visited[v] = 1;
					stk.push_back({ v, 0 });
				}
			} else {
				stk.pop_back();
				visitor(o, u, -1, d);
			}
		}
		if (source != -1) break;
	}
}

}
}
