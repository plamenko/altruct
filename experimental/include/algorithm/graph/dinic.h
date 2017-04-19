#pragma once

#include <queue>
#include <vector>
#include <limits>

namespace altruct {
namespace graph {

/**
 * Uses Dinic's algorithm to find the maximum flow.
 */
template<typename T>
class dinic {
public:
    T infinity;
	std::vector<std::vector<int>> adjl; // adjl[a] is the adjacency list of a.
	std::vector<std::vector<T>> cap;    // cap[a][b] is the capacity from a to b.
	std::vector<std::vector<T>> flow;   // flow[a][b] is the occupied flow from a to b.
	std::vector<int> level;             // level[a] is the level in level graph of a.

    dinic(const std::vector<std::vector<T>>& cap, T infinity = std::numeric_limits<T>::max()) : cap(cap), infinity(infinity) {
		adjl.resize(cap.size());
		for (int u = 0; u < (int)cap.size(); u++) {
			for (int v = 0; v < u; v++) {
                if (cap[u][v] > 0 || cap[v][u] > 0) {
					adjl[u].push_back(v);
					adjl[v].push_back(u);
				}
			}
		}
	}

	T calc_max_flow(int source, int sink) {
		T f = 0;
        if (source == sink) return f;
        flow.assign(adjl.size(), std::vector<T>(adjl.size()));
		while (build_level_graph(source, sink)) {
			f += construct_blocking_flow(source, sink);
		}
		return f;
	}

private:
	bool build_level_graph(int source, int sink) {
		std::queue<int> que;
		level.assign(adjl.size(), 0);
		level[source] = 1;
		que.push(source);
		while (!que.empty()) {
			int u = que.front(); que.pop();
			for (int i = 0; i < adjl[u].size(); ++i) {
				int v = adjl[u][i];
				if (!level[v] && (cap[u][v] > flow[u][v] || flow[v][u] > 0)) {
					level[v] = level[u] + 1;
					que.push(v);
				}
			}
		}
		return level[sink] != 0;
	}

	T construct_blocking_flow(int source, int sink) {
		T ret = 0;
		std::vector<int> prev(adjl.size());
		std::vector<int> visited(adjl.size());
		std::vector<std::pair<int, int>> stk;
		visited[source] = 1;
		stk.push_back({ source, 0 });
		while (!stk.empty()) {
			int u = stk.back().first;
			int i = stk.back().second;
			if (u == sink) {
				T f = get_path_flow(source, sink, prev);
				int bottleneck = update_path_flow(source, sink, prev, f);
				while (!stk.empty() && stk.back().first != bottleneck) stk.pop_back();
				ret += f;
				continue;
			}
			if (i < (int)adjl[u].size()) {
				int v = adjl[u][i];
				stk.back().second++;
				if (visited[v] || level[v] != level[u] + 1) continue;
				if (cap[u][v] > flow[u][v] || flow[v][u] > 0) {
					visited[v] = 1;
					stk.push_back({ v, 0 });
					prev[v] = u;
				}
			} else {
				stk.pop_back();
			}
		}
		return ret;
	}

	T get_path_flow(int source, int sink, const std::vector<int>& prev) {
        T f = infinity;
		for (int u, v = sink; u = prev[v], v != source; v = u) {
			T df = cap[u][v] - flow[u][v];
			f = std::min(f, df > 0 ? df : flow[v][u]);
		}
		return f;
	}

	int update_path_flow(int source, int sink, const std::vector<int>& prev, T f) {
		int bottleneck = 0;
		for (int u, v = sink; u = prev[v], v != source; v = u) {
			T df = cap[u][v] - flow[u][v];
			if (df > 0) {
				flow[u][v] += f;
				if (cap[u][v] == flow[u][v]) bottleneck = u;
			} else {
				flow[v][u] -= f;
				if (flow[v][u] == 0) bottleneck = u;
			}
		}
		return bottleneck;
	}
};

} // graph
} //altruct
