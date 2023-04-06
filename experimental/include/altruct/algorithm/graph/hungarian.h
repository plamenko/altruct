#pragma once

#include "altruct/structure/graph/graph.h"

#include <vector>
#include <algorithm>

namespace altruct {
namespace graph {

/**
 * Minimum/Maximum-cost bipartite matching.
 * Complexity: O(n^3)
 * Usage:
 *   hungarian<> h;
 *   h.resize(n);
 *   for (x)
 *     for (y)
 *       h.cost(x, y) = ...;
 *   auto r = h.calc_max_cost();
 */
template<typename I = int>
class hungarian {
public:
    // n workers and n jobs
    void resize(int n) {
        this->n = n;
        while ((1 << logSZ) < n) logSZ++;
        int sz = 1 << logSZ;
        c.resize(sz * sz);
        S.resize(sz);
        T.resize(sz);
        lx.resize(sz);
        ly.resize(sz);
        slack.resize(sz);
        slackx.resize(sz);
        prev.resize(sz);
        q.resize(sz);
        xy.resize(sz);
        yx.resize(sz);
    }

    I& cost(int x, int y) {
        return c[(x << logSZ) + y];
    }

    I calc_min_cost() {
        negate_cost();
        I ret = -calc_max_cost();
        negate_cost();
        return ret;
    }

    I calc_max_cost() {
        calc_max_cost_impl();
        I ret = 0;
        for (int x = 0; x < n; x++) {
            ret += cost(x, xy[x]);
        }
        return ret;
    }

    const std::vector<int>& get_matches_for_x() const {
        return xy;
    }

    const std::vector<int>& get_matches_for_y() const {
        return yx;
    }

private:
    // input
    int n = 0, logSZ = 0;          // size
    std::vector<I> c;              // n*n cost matrix
    // working state
    std::vector<int> S, T;         // sets S and T in algorithm
    std::vector<I> lx, ly;         // labels of X and Y parts
    std::vector<I> slack;          // as in the algorithm description
    std::vector<int> slackx;       // slackx[y] such a vertex, that l(slackx[y]) + l(y) - w(slackx[y],y) = slack[y]
    std::vector<int> prev;         // array for memorizing alternating paths
    std::vector<int> q;            // bfs queue
    // resulting matching
    std::vector<int> xy;           // xy[x] - vertex that is matched with x,
    std::vector<int> yx;           // yx[y] - vertex that is matched with y

    void negate_cost() {
        for (int x = 0; x < n; x++) {
            for (int y = 0; y < n; y++) {
                cost(x, y) = -cost(x, y);
            }
        }
    }

    void add_to_tree(int x, int prevx) {
        // x - current vertex, prevx - vertex from X before x in the alternating path,
        // so we add edges (prevx, xy[x]), (xy[x], x)
        S[x] = 1;                       // add x to S
        prev[x] = prevx;                // we need this when augmenting
        for (int y = 0; y < n; y++) {   // update slacks, because we add new vertex to S
            I s = lx[x] + ly[y] - cost(x, y);
            if (s < slack[y]) {
                slack[y] = s;
                slackx[y] = x;
            }
        }
    }

    void calc_max_cost_impl() {
        int x, y;

        for (x = 0; x < n; x++) xy[x] = -1;
        for (y = 0; y < n; y++) yx[y] = -1;
        for (y = 0; y < n; y++) ly[y] = 0;
        for (x = 0; x < n; x++) {
            lx[x] = cost(x, 0);
            for (y = 0; y < n; y++) {
                lx[x] = max(lx[x], cost(x, y));
            }
        }

        for (int max_match = 0; max_match < n; max_match++) {
            for (x = 0; x < n; x++) S[x] = 0;                        // init set S, T
            for (y = 0; y < n; y++) T[y] = 0;
            for (x = 0; x < n; x++) prev[x] = -1;                    // init prev - for the alternating tree

            int wr = 0, rd = 0;                                      // bfs queue read and write indices
            for (x = 0; x < n; x++) {                                // finding root of the tree (x)
                if (xy[x] == -1) {
                    q[wr++] = x;
                    prev[x] = -2;
                    S[x] = 1;
                    break;
                }
            }

            for (y = 0; y < n; y++) {                                // initializing slack array
                slack[y] = lx[x] + ly[y] - cost(x, y);
                slackx[y] = x;
            }

            for (;;) {                                               // main cycle
                while (rd < wr) {                                    // building tree with bfs cycle
                    x = q[rd++];                                     // current vertex from X part
                    for (y = 0; y < n; y++) {                        // iterate through all edges in equality graph
                        if (!T[y] && cost(x, y) == lx[x] + ly[y]) {
                            if (yx[y] == -1) break;                  // an exposed vertex in Y found, so augmenting path exists!
                            T[y] = 1;                                // else just add y to T,
                            q[wr++] = yx[y];                         // add vertex yx[y], which is matched with y, to the queue
                            add_to_tree(yx[y], x);                   // add edges (x,y) and (y,yx[y]) to the tree
                        }
                    }
                    if (y < n) break;                                // augmenting path found!
                }
                if (y < n) break;                                    // augmenting path found!

                // augmenting path not found, so improve labeling
                I delta = 0;
                for (y = 0; y < n; y++)                              // initial delta (instead of INF)
                    if (!T[y]) { delta = slack[y]; break; }
                for (y = 0; y < n; y++)                              // calculate delta using slack
                    if (!T[y]) delta = min(delta, slack[y]);
                for (y = 0; y < n; y++)                              // update slack array
                    if (!T[y]) slack[y] -= delta;
                for (x = 0; x < n; x++)                              // update X labels
                    if (S[x]) lx[x] -= delta;
                for (y = 0; y < n; y++)                              // update Y labels
                    if (T[y]) ly[y] += delta;

                wr = rd = 0;
                for (y = 0; y < n; y++) {
                    // in this cycle we add edges that were added to the equality graph as a
                    // result of improving the labeling, we add edge (slackx[y], y) to the tree if
                    // and only if !T[y] && slack[y] == 0, also with this edge we add another one
                    // (y, yx[y]) or augment the matching, if y was exposed
                    if (!T[y] && slack[y] == 0) {
                        if (yx[y] == -1) {                            // exposed vertex in Y found - augmenting path exists!
                            x = slackx[y];
                            break;
                        }
                        else {
                            T[y] = 1;                                // else just add y to T,
                            if (!S[yx[y]]) {
                                q[wr++] = yx[y];                     // add vertex yx[y], which is matched with y, to the queue
                                add_to_tree(yx[y], slackx[y]);       // and add edges (x,y) and (y, yx[y]) to the tree
                            }
                        }
                    }
                }
                if (y < n) break;                                    // augmenting path found!
            }

            for (int ty; x != -2; x = prev[x], y = ty) {             // inverse edges along augmenting path
                ty = xy[x];
                yx[y] = x;
                xy[x] = y;
            }
        }
    }
};

} // graph
} // altruct
