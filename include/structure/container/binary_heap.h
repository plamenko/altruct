#pragma once

#include <vector>
#include <iterator>
#include <functional>

namespace altruct {
namespace container {

/**
 * Binary heap.
 *
 * Space complexity: `O(n)`.
 * Time complexities:
 *   build:     `O(n)`
 *   insert:    `O(log n)`
 *   pop_front: `O(log n)`
 *   front:     `O(1)`
 *
 * param T   - element type
 * param cmp - comparison functor
 */
template<typename T>
class binary_heap {
public:
    std::vector<T> v;
    std::function<bool(T, T)> cmp;

    binary_heap(std::function<bool(T, T)> cmp = std::less<T>()) : cmp(cmp)  {}
    binary_heap(size_t sz, std::function<bool(T, T)> cmp = std::less<T>()) : v(sz), cmp(cmp) {}
    binary_heap(std::vector<T>& v, std::function<bool(T, T)> cmp = std::less<T>()) : v(v), cmp(cmp) {
        rebuild();
    }
    template<typename It>
    binary_heap(It begin, It end, std::function<bool(T, T)> cmp = std::less<T>()) : v(begin, end), cmp(cmp) {
        rebuild();
    }

    size_t size() const {
        return v.size();
    }

    void insert(const T& val) {
        v.push_back(val);
        heapify_up(v.size() - 1);
    }

    const T& front() const {
        return v.front();
    }

    void pop_front() {
        std::swap(v.front(), v.back());
        v.pop_back();
        heapify_down(0, v.size());
    }

    /** rebuilds the heap */
    void rebuild() {
        for (size_t index = v.size() / 2; index > 0; index--) {
            heapify_down(index - 1, v.size());
        }
    }

    /** converts the heap to a sorted array (heap order is destroyed) */
    void sort() {
        for (size_t index = size(); index > 0; index--) {
            std::swap(v[0], v[index - 1]);
            heapify_down(0, index - 1);
        }
        reverse(v.begin(), v.end());
    }

private:
    void heapify_up(size_t index) {
        while (index > 0) {
            size_t parent = (index - 1) / 2;
            if (!cmp(v[index], v[parent])) break;
            std::swap(v[index], v[parent]);
            index = parent;
        }
    }

    void heapify_down(size_t index, size_t sz) {
        size_t best = index;
        while (true) {
            size_t left = index * 2 + 1, right = index * 2 + 2;
            if (left < sz && cmp(v[left], v[best])) best = left;
            if (right < sz && cmp(v[right], v[best])) best = right;
            if (best == index) break;
            std::swap(v[index], v[best]);
            index = best;
        }
    }
};

} // container
} // altruct
