#pragma once

#include <vector>
#include <cstddef>

namespace altruct {
namespace graph {

class disjoint_set {
    struct atom {
        size_t parent, rank, count;
        atom(size_t parent = 0, size_t rank = 0, size_t count = 1) : parent(parent), rank(rank), count(count) {}
    };

    size_t dist;
    std::vector<atom> va;

    void ensure(size_t sz) {
        while (va.size() < sz) {
            va.push_back(atom(va.size(), 0, 1));
            dist++;
        }
    }

public:
    disjoint_set(size_t sz = 0) {
        clear(sz);
    }

    void clear(size_t sz = 0) {
        dist = 0;
        va.clear();
        ensure(sz);
    }

    bool unite(size_t x, size_t y) {
        x = find(x);
        y = find(y);
        if (x == y) {
            return false;
        }
        if (va[x].rank < va[y].rank) {
            va[x].parent = y, va[y].count += va[x].count;
        }
        else {
            va[y].parent = x, va[x].count += va[y].count;
        }
        if (va[x].rank == va[y].rank) {
            va[x].rank++;
        }
        dist--;
        return true;
    }

    size_t find(size_t x) {
        ensure(x + 1);
        size_t y = x, t;
        while (y != va[y].parent) {
            y = va[y].parent;
        }
        while (x != va[x].parent) {
            t = x, x = va[t].parent, va[t].parent = y;
        }
        return y;
    }

    size_t count(size_t x) {
        return va[find(x)].count;
    }

    size_t distinct() const {
        return dist;
    }

    size_t size() const {
        return va.size();
    }
};

} // graph
} // altruct
