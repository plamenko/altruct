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
    T infinity;
    std::vector<std::vector<int>> adjl; // adjl[a] is the adjacency list of a.
    std::vector<std::vector<T>> cap;    // cap[a][b] is the capacity from a to b.
    std::vector<std::vector<T>> flow;   // flow[a][b] is the occupied flow from a to b.
    
    std::vector<int> height;            // height of node.
    std::vector<T> excess;              // flow into node minus flow from node.
    std::vector<int> seen;              // neighbours seen since last relabel.

    push_relabel_flow(const std::vector<std::vector<T>>& cap, T infinity = std::numeric_limits<T>::max()) : cap(cap), infinity(infinity) {
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
        int n = (int)adjl.size();
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

        deque<int> que;
        for (int v = 0; v < n; v++) {
            if (v != source && v != sink) que.push_back(v);
        }

        auto it = que.begin();
        while (it != que.end()) {
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
        int n = (int)adjl.size();
        int min_height = n + 1;
        for (int v = 0; v < n; v++) {
            if (cap[u][v] - flow[u][v] > 0){
                min_height = std::min(min_height, height[v]);
                height[u] = min_height + 1;
            }
        }
    }

    void discharge(int u) {
        int n = (int)adjl.size();
        while (excess[u] > 0) {
            if (seen[u] < n) { // check next neighbour
                int v = seen[u];
                if (cap[u][v] - flow[u][v] > 0 && height[u] > height[v]) {
                    push(u, v);
                } else {
                    seen[u]++;
                }
            } else { // we have checked all neighbours. must relabel
                relabel(u);
                seen[u] = 0;
            }
        }
    }
};

} // graph
} //altruct
