#pragma once

#include "altruct/structure/container/binary_search_tree.h"

#include <cstdlib>
#include <functional>

namespace altruct {
namespace container {

/**
 * Bidirectional iterator.
 */
template<typename T, typename IT>
struct treap_iterator : public std::iterator<std::bidirectional_iterator_tag, T> {
    IT it;
    treap_iterator() { }
    treap_iterator(IT it) : it(it) { }
    int count() const { return it.count(); }
    int size() const { return it.size(); }
    T& operator * () { return *it; }
    T* operator -> () { return &*it; }
    bool operator == (const treap_iterator& rhs) const { return it == rhs.it; }
    bool operator != (const treap_iterator& rhs) const { return it != rhs.it; }
    treap_iterator& operator--() { --it; return *this; }
    treap_iterator operator--(int) { return it--; }
    treap_iterator& operator++() { ++it; return *this; }
    treap_iterator operator++(int) { return it++; }

    // `pos` and `add` are not suitable when duplicate mode is set to COUNT
    int pos() const { return it.pos(); }
    treap_iterator add(difference_type off) { return it.add(off); }

    treap_iterator parent() const { return it.parent(); }
    treap_iterator right() const { return it.right(); }
    treap_iterator left() const { return it.left(); }
};

/**
 * Bidirectional const iterator.
 */
template<typename T, typename CIT, typename IT>
struct treap_const_iterator : public std::iterator<std::bidirectional_iterator_tag, const T> {
    CIT it;
    treap_const_iterator() { }
    treap_const_iterator(CIT it) : it(it) { }
    treap_const_iterator(treap_iterator<T, IT> it) : it(it.it) { }
    int count() const { return it.count(); }
    int size() const { return it.size(); }
    const T& operator * () const { return *it; }
    const T* operator -> () const { return &*it; }
    bool operator == (const treap_const_iterator& rhs) const { return it == rhs.it; }
    bool operator != (const treap_const_iterator& rhs) const { return it != rhs.it; }
    treap_const_iterator& operator--() { --it; return *this; }
    treap_const_iterator operator--(int) { return it--; }
    treap_const_iterator& operator++() { ++it; return *this; }
    treap_const_iterator operator++(int) { return it++; }

    // `pos` and `add` are not suitable when duplicate mode is set to COUNT
    int pos() const { return it.pos(); }
    treap_const_iterator add(difference_type off) { return it.add(off); }

    treap_const_iterator parent() const { return it.parent(); }
    treap_const_iterator right() const { return it.right(); }
    treap_const_iterator left() const { return it.left(); }
};

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
 * param RAND - random generator type
 * param CMP - comparison functor type
 * param ALLOC - allocator type
 */
template<typename K, typename T = K, int DUP = bst_duplicate_handling::IGNORE, typename CMP = std::less<K>, typename RAND = std::function<int()>, typename ALLOC = std::allocator<bst_node<T>>>
class treap {
protected:
    typedef binary_search_tree<K, T, DUP, CMP, ALLOC> bst_t;
    bst_t tree;
    RAND rnd;

public:
    typedef K key_type;
    typedef T value_type;
    typedef treap_iterator<T, typename bst_t::iterator> iterator;
    typedef treap_const_iterator<T, typename bst_t::const_iterator, typename bst_t::iterator> const_iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    treap(const CMP& cmp = CMP(), const RAND& rnd = rand, const ALLOC& alloc = ALLOC()) :
        tree(cmp, alloc), rnd(rnd) {
    }

    template<typename It>
    treap(It begin, It end, const CMP& cmp = CMP(), const RAND& rnd = rand, const ALLOC& alloc = ALLOC()) :
        treap(cmp, rnd, alloc) {
        for (It it = begin; it != end; ++it) {
            insert(*it);
        }
    }

    treap(std::initializer_list<T> list) :
        treap(list.begin(), list.end()) {
    }

    treap(treap&& rhs) :
        treap() {
        swap(rhs);
    }

    treap(const treap& rhs) :
        tree(rhs.tree),
        rnd(rhs.rnd) {
    }

    treap& operator=(treap&& rhs) {
        swap(rhs);
        return *this;
    }

    treap& operator=(const treap& rhs) {
        tree = rhs.tree;
        rnd = rhs.rnd;
        return *this;
    }

    void swap(treap& rhs) {
        tree.swap(rhs.tree);
        std::swap(rnd, rhs.rnd);
    }

    void clear() {
        tree.clear();
    }

    bool empty() const {
        return tree.empty();
    }

    int size() const {
        return tree.size();
    }

public: // iterators
    iterator root() { return tree.root(); }
    const_iterator root() const { return tree.root(); }
    const_iterator croot() const { return tree.croot(); }
    iterator begin() { return tree.begin(); }
    const_iterator begin() const { return tree.begin(); }
    const_iterator cbegin() const { return tree.cbegin(); }
    iterator end() { return tree.end(); }
    const_iterator end() const { return tree.end(); }
    const_iterator cend() const { return tree.cend(); }
    reverse_iterator rbegin() { return reverse_iterator(end()); }
    const_reverse_iterator rbegin() const { return const_reverse_iterator(end()); }
    const_reverse_iterator crbegin() const { return const_reverse_iterator(cend()); }
    reverse_iterator rend() { return reverse_iterator(begin()); }
    const_reverse_iterator rend() const { return const_reverse_iterator(begin()); }
    const_reverse_iterator crend() const { return const_reverse_iterator(cbegin()); }

public: // relational operators
    bool operator == (const treap& rhs) const { return tree == rhs.tree; }
    bool operator < (const treap& rhs) const { return tree < rhs.tree; }
    bool operator != (const treap& rhs) const { return !(*this == rhs); }
    bool operator >(const treap& rhs) const { return (rhs < *this); }
    bool operator <= (const treap& rhs) const { return !(rhs < *this); }
    bool operator >= (const treap& rhs) const { return !(*this < rhs); }

public: // query & update
    int count_less_or_equal(const K& key) const {
        return tree.count_less_or_equal(key);
    }

    int count_less(const K& key) const {
        return tree.count_less(key);
    }

    int count(const K& key) const {
        return tree.count(key);
    }

    const_iterator find_kth(int k) const {
        return tree.find_kth(k);
    }

    iterator find_kth(int k) {
        return tree.find_kth(k);
    }

    const_iterator find(const K& key) const {
        return tree.find(key);
    }

    iterator find(const K& key) {
        return tree.find(key);
    }

    const_iterator lower_bound(const K& key) const {
        return tree.lower_bound(key);
    }

    iterator lower_bound(const K& key) {
        return tree.lower_bound(key);
    }

    const_iterator upper_bound(const K& key) const {
        return tree.upper_bound(key);
    }

    iterator upper_bound(const K& key) {
        return tree.upper_bound(key);
    }

    std::pair<const_iterator, const_iterator> equal_range(const K& key) const {
        return tree.equal_range(key);
    }

    std::pair<iterator, iterator> equal_range(const K& key) {
        return tree.equal_range(key);
    }

    iterator insert(const T& val, int cnt = 1) {
        return retrace_up(tree.insert(val, cnt));
    }

    iterator insert_before(const_iterator it, const T& val, int cnt = 1) {
        return retrace_up(tree.insert_before(it.it, val, cnt));
    }

    iterator erase(const K& key, int cnt = std::numeric_limits<int>::max()) {
        if (DUP == bst_duplicate_handling::STORE) {
            return erase(lower_bound(key), upper_bound(key));
        } else {
            return erase(find(key), cnt);
        }
    }

    iterator erase(const_iterator b, const_iterator e, int cnt = std::numeric_limits<int>::max()) {
        iterator res;
        while (b != e) res = erase(b++, cnt);
        return res;
    }

    iterator erase(const_iterator it, int cnt = std::numeric_limits<int>::max()) {
        return tree.erase(retrace_down(it.it), cnt);
    }

protected:
    typename bst_t::iterator retrace_up(typename bst_t::iterator it) {
        if (it == tree.end()) return it;
        it.balance() = rnd();
        while (it.balance() < it.parent().balance()) {
            if (it.parent().left() == it) {
                tree.rotate_right(it.parent());
            } else {
                tree.rotate_left(it.parent());
            }
        }
        return it;
    }

    typename bst_t::const_iterator retrace_down(typename bst_t::const_iterator it) {
        if (it == tree.end()) return it;
        while (it.left() != tree.end() && it.right() != tree.end()) {
            if (it.left().balance() < it.right().balance()) {
                tree.rotate_right(it);
            } else {
                tree.rotate_left(it);
            }
        }
        return it;
    }
};

} // container
} // altruct
