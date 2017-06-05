#pragma once

#include "altruct/structure/container/binary_search_tree.h"

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
template<typename K, typename T = K, int DUP = bst_duplicate_handling::IGNORE, typename CMP = std::less<K>, typename RAND = std::function<int()>, typename ALLOC = std::allocator<bst_node<T>>>
class treap : public binary_search_tree<K, T, DUP, CMP, ALLOC> {
protected:
    typedef binary_search_tree<K, T, DUP, CMP, ALLOC> bst_t;
    typedef typename bst_t::node_ptr node_ptr;
    typedef typename bst_t::const_node_ptr const_node_ptr;
public:
    typedef K key_type;
    typedef T value_type;
    typedef typename bst_t::iterator iterator;
    typedef typename bst_t::const_iterator const_iterator;
    typedef typename bst_t::reverse_iterator reverse_iterator;
    typedef typename bst_t::const_reverse_iterator const_reverse_iterator;

public:
    RAND rnd;

    treap(const CMP& cmp = CMP(), const RAND& rnd = rand, const ALLOC& alloc = ALLOC()) :
        bst_t(cmp, alloc), rnd(rnd) {
    }

    template<typename It>
    treap(It begin, It end, const CMP& cmp = CMP(), const RAND& rnd = rand, const ALLOC& alloc = ALLOC()) :
        bst_t(cmp, alloc), rnd(rnd) {
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
        treap rhs_copy(rhs);
        swap(rhs_copy);
        return *this;
    }

    void swap(treap& rhs) {
        bst_t::swap(rhs);
        std::swap(rnd, rhs.rnd);
    }

    iterator insert(const T& val, int cnt = 1) {
        return retrace_up(bst_t::insert(val, cnt));
    }

    iterator erase(const K& key, int cnt = std::numeric_limits<int>::max()) {
        if (DUP == bst_duplicate_handling::STORE) {
            return erase(bst_t::lower_bound(key), bst_t::upper_bound(key));
        } else {
            return erase(bst_t::find(key), cnt);
        }
    }

    iterator erase(const_iterator b, const_iterator e, int cnt = std::numeric_limits<int>::max()) {
        while (b != e) erase(b++);
        return bst_t::remove_const(b);
    }

    iterator erase(const_iterator it, int cnt = std::numeric_limits<int>::max()) {
        return bst_t::erase(retrace_down(it), cnt);
    }

protected:
    iterator retrace_up(const_iterator it) {
        auto ptr = bst_t::remove_const(it);
        if (ptr->is_nil()) return ptr;
        ptr->balance = rnd();
        while (ptr->balance < ptr->parent->balance) {
            if (ptr->parent->left == ptr) {
                bst_t::rotate_right(ptr->parent);
            } else {
                bst_t::rotate_left(ptr->parent);
            }
        }
        return ptr;
    }

    iterator retrace_down(const_iterator it) {
        auto ptr = bst_t::remove_const(it);
        if (ptr->is_nil()) return ptr;
        while (!ptr->left->is_nil() && !ptr->right->is_nil()) {
            if (ptr->left->balance < ptr->right->balance) {
                bst_t::rotate_right(ptr);
            } else {
                bst_t::rotate_left(ptr);
            }
        }
        return ptr;
    }
};

} // container
} // altruct
