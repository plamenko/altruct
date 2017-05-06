#pragma once

#include "structure/container/binary_search_tree.h"

#include <cstdlib>
#include <functional>

namespace altruct {
namespace container {

/**
 * Treap.
 *
 * Note, balancing is performed based on random numbers.
 *
 * Space complexity: `O(n)`.
 * Time complexities:
 *   find:   `O(h)`
 *   insert: `O(h)`
 *   erase:  `O(h)`
 * Where `h` is height that is proportional to `log(n)` with
 * very high probability.
 *
 * param K   - key type
 * param T   - value type (for maps this is going to be bst_entry<K, V>)
 * param DUP - duplicate handling mode
 * param CMP - comparison functor type
 * param ALLOC - allocator type
 */
template<typename K, typename T = K, int DUP = bst_duplicate_handling::IGNORE, typename CMP = std::less<K>, typename RAND = std::function<int()>, typename ALLOC = allocator<bst_node<T>>>
class treap : public binary_search_tree<K, T, DUP, CMP, ALLOC> {
public:
    RAND rnd;

    treap(const CMP& cmp = CMP(), const RAND& rnd = rand, const ALLOC& alloc = ALLOC()) :
        binary_search_tree(cmp, alloc), rnd(rnd) {
    }

    template<typename It>
    treap(It begin, It end, const CMP& cmp = CMP(), const RAND& rnd = rand, const ALLOC& alloc = ALLOC()) :
        binary_search_tree(cmp, alloc), rnd(rnd) {
        for (It it = begin; it != end; ++it) {
            insert(*it);
        }
    }

    treap(std::initializer_list<T> list) :
        treap(list.begin(), list.end()) {}

    treap(treap&& rhs) :
        treap() {
        swap(rhs);
    }

    treap(const treap& rhs) :
        treap(rhs.cbegin(), rhs.cend(), rhs.cmp, rhs.rnd, rhs.alloc) {
    }

    treap& operator=(treap&& rhs) {
        swap(rhs);
        return *this;
    }

    treap& operator=(const treap& rhs) {
        swap(treap(rhs));
        return *this;
    }

    iterator insert(const T& val, int cnt = 1) {
        auto it = binary_search_tree::insert(val, cnt);
        auto ptr = remove_const(it);
        ptr->balance = rnd();
        while (ptr->balance < ptr->parent->balance) {
            if (ptr->parent->left == ptr) {
                rotate_right(ptr->parent);
            } else {
                rotate_left(ptr->parent);
            }
        }
        return it;
    }

    iterator erase(const K& key, int cnt = std::numeric_limits<int>::max()) {
        if (DUP == bst_duplicate_handling::STORE) {
            return erase(lower_bound(key), upper_bound(key));
        } else {
            return erase(find(key), cnt);
        }
    }

    iterator erase(const_iterator b, const_iterator e, int cnt = std::numeric_limits<int>::max()) {
        while (b != e) erase(b++);
        return remove_const(b);
    }

    iterator erase(const_iterator it, int cnt = std::numeric_limits<int>::max()) {
        auto ptr = remove_const(it);
        if (ptr->is_nil()) return ptr;
        while (!ptr->left->is_nil() && !ptr->right->is_nil()) {
            if (ptr->left->balance < ptr->right->balance) {
                rotate_right(ptr);
            } else {
                rotate_left(ptr);
            }
        }
        return binary_search_tree::erase(it, cnt);
    }
};

} // container
} // altruct
