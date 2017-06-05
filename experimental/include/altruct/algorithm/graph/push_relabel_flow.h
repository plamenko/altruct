#pragma once

#include <deque>
#include <vector>
#include <limits>

namespace altruct {
namespace graph {

/**
 * Uses Push-relabel algorithm to find the maximum flow.
 */
template<typename T>
class push_relabel_flow {
public:
    enum selection_rule { RELABEL_TO_FRONT, LARGEST_LABEL };

    T infinity;
    std::vector<std::vector<int>> adjl; // adjl[a] is the adjacency list of a.
    std::vector<std::vector<T>> cap;    // cap[a][b] is the capacity from a to b.
    std::vector<std::vector<T>> flow;   // flow[a][b] is the occupied flow from a to b.
    
    std::vector<int> height;            // height of node.
    std::vector<T> excess;              // flow into node minus flow from node.
    std::vector<int> seen;              // neighbours seen since last relabel.

    // infinity = std::numeric_limits<T>::max() is no good because of overflows!!!
    push_relabel_flow(const std::vector<std::vector<T>>& cap, T infinity) : cap(cap), infinity(infinity) {
        int n = (int)cap.size();
        adjl.resize(n);
        for (int u = 0; u < n; u++) {
            for (int v = 0; v < u; v++) {
                if (cap[u][v] > 0 || cap[v][u] > 0) {
                    adjl[u].push_back(v);
                    adjl[v].push_back(u);
                }
            }
        }
    }

    T calc_max_flow(int source, int sink, int rule = RELABEL_TO_FRONT) {
        int n = (int)cap.size();
        if (source == sink) return 0;
        
        flow.assign(n, std::vector<T>(n, 0));
        height.assign(n, 0);
        excess.assign(n, 0);
        seen.assign(n, 0);

        height[source] = n;        // longest path from source to sink is less than n long
        excess[source] = infinity; // send as much flow as possible to neighbours of source
        for (int v = 0; v < n; v++) {
            push(source, v);
        }

        if (rule == RELABEL_TO_FRONT) {
            std::deque<int> que;
            for (int v = 0; v < n; v++) {
                if (v != source && v != sink) que.push_back(v);
            }
            for (auto it = que.begin(); it != que.end();) {
                int u = *it;
                int old_height = height[u];
                discharge(u);
                if (height[u] > old_height) {
                    // relabel to front strategy
                    que.erase(it);
                    que.push_front(u);
                    it = que.begin();
                } else {
                    ++it;
                }
            }
        } else if (rule == LARGEST_LABEL) {
            std::vector<int> que;
            for (int v = 0; v < n; v++) {
                if (v != source && v != sink) que.push_back(v);
            }
            for (auto it = que.begin(); it != que.end();) {
                int u = *it;
                int old_height = height[u];
                discharge(u);
                if (height[u] > old_height) {
                    // largest-label strategy
                    // TODO: it is not optimal to perform a full sort here
                    sort(que.begin(), que.end(), [&](int u1, int u2){ return height[u1] > height[u2]; });
                    it = que.begin();
                } else {
                    ++it;
                }
            }
        }

        T f = 0;
        for (auto df : flow[source]) {
            f += df;
        }
        return f;
    }

private:

    void push(int u, int v) {
        auto send = std::min(excess[u], cap[u][v] - flow[u][v]);
        flow[u][v] += send;
        flow[v][u] -= send;
        excess[u] -= send;
        excess[v] += send;
    }

    void relabel(int u) {
        // find smallest new height making a push possible, if such a push is possible at all
        int min_height = std::numeric_limits<int>::max();
        for (auto v : adjl[u]) {
            if (cap[u][v] - flow[u][v] > 0){
                min_height = std::min(min_height, height[v]);
                height[u] = min_height + 1;
            }
        }
    }

    void discharge(int u) {
        int it = seen[u];
        while (excess[u] > 0) {
            if (it < adjl[u].size()) { // check next neighbour
                int v = adjl[u][it];
                if (cap[u][v] - flow[u][v] > 0 && height[u] > height[v]) {
                    push(u, v);
                } else {
                    it++;
                }
            } else { // we have checked all neighbours. must relabel
                relabel(u);
                it = 0;
            }
        }
        seen[u] = it;
    }
};

} // graph
} // altruct
