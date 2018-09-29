#pragma once

#include "altruct/algorithm/collections/collections.h"
#include "altruct/structure/container/treap.h"

#include <cstdlib>
#include <functional>

namespace altruct {
namespace container {

/**
 * Random access iterator.
 */
template<typename T, typename IT>
struct rope_iterator : public std::iterator<std::random_access_iterator_tag, T> {
    typedef std::iterator<std::random_access_iterator_tag, T> it_t;
    typedef typename it_t::difference_type difference_type;
    IT it;
    rope_iterator() { }
    rope_iterator(IT it) : it(it) { }
    int count() const { return it.count(); }
    int pos() const { return it.pos(); }
    T& operator * () { return *it; }
    T* operator -> () { return &*it; }
    bool operator == (const rope_iterator& rhs) const { return it == rhs.it; }
    bool operator != (const rope_iterator& rhs) const { return it != rhs.it; }
    bool operator < (const rope_iterator& rhs) const { return (*this - rhs) < 0; }
    bool operator > (const rope_iterator& rhs) const { return (*this - rhs) > 0; }
    bool operator <= (const rope_iterator& rhs) const { return (*this - rhs) <= 0; }
    bool operator >= (const rope_iterator& rhs) const { return (*this - rhs) >= 0; }
    rope_iterator& operator--() { --it; return *this; }
    rope_iterator operator--(int) { return it--; }
    rope_iterator& operator++() { ++it; return *this; }
    rope_iterator operator++(int) { return it++; }
    rope_iterator& operator+=(difference_type off) { return *this = it.add(off); }
    rope_iterator operator+(difference_type off) const { rope_iterator it = *this; return it += off; }
    rope_iterator& operator-=(difference_type off) { return *this += -off; }
    rope_iterator operator-(difference_type off) const { rope_iterator it = *this; return it -= off; }
    difference_type operator-(const rope_iterator& rhs) const { return pos() - rhs.pos(); }
    T& operator[](difference_type off) const { return *(*this + off); }
};

/**
 * Random access const iterator.
 */
template<typename T, typename CIT, typename IT>
struct rope_const_iterator : public std::iterator<std::random_access_iterator_tag, const T> {
    typedef std::iterator<std::random_access_iterator_tag, const T> it_t;
    typedef typename it_t::difference_type difference_type;
    CIT it;
    rope_const_iterator() { }
    rope_const_iterator(CIT it) : it(it) { }
    rope_const_iterator(rope_iterator<T, IT> it) : it(it.it) { }
    int count() const { return it.count(); }
    int pos() const { return it.pos(); }
    const T& operator * () const { return *it; }
    const T* operator -> () const { return &*it; }
    bool operator == (const rope_const_iterator& rhs) const { return it == rhs.it; }
    bool operator != (const rope_const_iterator& rhs) const { return it != rhs.it; }
    bool operator < (const rope_const_iterator& rhs) const { return (*this - rhs) < 0; }
    bool operator > (const rope_const_iterator& rhs) const { return (*this - rhs) > 0; }
    bool operator <= (const rope_const_iterator& rhs) const { return (*this - rhs) <= 0; }
    bool operator >= (const rope_const_iterator& rhs) const { return (*this - rhs) >= 0; }
    rope_const_iterator& operator--() { --it; return *this; }
    rope_const_iterator operator--(int) { return it--; }
    rope_const_iterator& operator++() { ++it; return *this; }
    rope_const_iterator operator++(int) { return it++; }
    rope_const_iterator& operator+=(difference_type off) { return *this = it.add(off); }
    rope_const_iterator operator+(difference_type off) const { rope_const_iterator it = *this; return it += off; }
    rope_const_iterator& operator-=(difference_type off) { return *this += -off; }
    rope_const_iterator operator-(difference_type off) const { rope_const_iterator it = *this; return it -= off; }
    difference_type operator-(const rope_const_iterator& rhs) const { return pos() - rhs.pos(); }
    const T& operator[](difference_type off) const { return *(*this + off); }
};

/**
 * Rope (array tree)
 *
 * Provides random-index operations in logarithmic time.
 *
 * Space complexity: `O(n)`.
 * Time complexities:
 *   at        is O(log n)
 *   push_back is O(log n)
 *   insert    is O(log n)
 *   erase     is O(log n)
 *
 * param T   - value type
 * param RAND - random generator type
 * param ALLOC - allocator type
 */
template<typename T, typename RAND = std::function<int()>, typename ALLOC = std::allocator<bst_node<T>>>
class rope {
protected:
    // comparator_t treats all the elements as equal so that the underlying tree doesn't reorder them
    struct comparator_t { bool operator()(const T&, const T&) const { return false; } };
    typedef treap<T, T, bst_duplicate_handling::STORE, comparator_t, RAND, ALLOC> treap_t;
    treap_t tree;

public:
    typedef T value_type;
    typedef rope_iterator<T, typename treap_t::iterator> iterator;
    typedef rope_const_iterator<T, typename treap_t::const_iterator, typename treap_t::iterator> const_iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    rope(const RAND& rnd = rand, const ALLOC& alloc = ALLOC()) :
        tree(comparator_t(), rnd, alloc) {
    }

    template<typename It>
    rope(It begin, It end, const RAND& rnd = rand, const ALLOC& alloc = ALLOC()) :
        rope(rnd, alloc) {
        for (It it = begin; it != end; ++it) {
            push_back(*it);
        }
    }

    rope(std::initializer_list<T> list) :
        rope(list.begin(), list.end()) {
    }

    rope(rope&& rhs) : rope() {
        swap(rhs);
    }

    rope(const rope& rhs) :
        tree(rhs.tree) {
    }

    rope& operator=(rope&& rhs) {
        swap(rhs);
        return *this;
    }

    rope& operator=(const rope& rhs) {
        tree = rhs.tree;
        return *this;
    }

    void swap(rope& rhs) {
        tree.swap(rhs.tree);
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
    bool operator == (const rope& rhs) const { return collections::compare(cbegin(), cend(), rhs.cbegin(), rhs.cend()) == 0; }
    bool operator < (const rope& rhs) const { return collections::compare(cbegin(), cend(), rhs.cbegin(), rhs.cend()) < 0; }
    bool operator != (const rope& rhs) const { return !(*this == rhs); }
    bool operator >  (const rope& rhs) const { return (rhs < *this); }
    bool operator <= (const rope& rhs) const { return !(rhs < *this); }
    bool operator >= (const rope& rhs) const { return !(*this < rhs); }

public: // query & update
    const_iterator find_kth(int k) const {
        return tree.find_kth(k);
    }

    iterator find_kth(int k) {
        return tree.find_kth(k);
    }

    void push_back(const T& val) {
        insert(size(), val);
    }

    void insert(int pos, const T& val) {
        insert(find_kth(pos), val);
    }

    void insert(const_iterator it, const T& val) {
        tree.insert_before(it.it, val, 1);
    }

    void erase(int b, int e) {
        erase(find_kth(b), find_kth(e));
    }

    void erase(const_iterator b, const_iterator e) {
        while (b != e) erase(b++);
    }

    void erase(int pos) {
        erase(find_kth(pos));
    }

    void erase(const_iterator it) {
        tree.erase(it.it, 1);
    }

    const T& at(int pos) const {
        return *find_kth(pos);
    }

    T& at(int pos) {
        return *find_kth(pos);
    }

    const T& operator[](int pos) const {
        return at(pos);
    }

    T& operator[](int pos) {
        return at(pos);
    }
};

} // container
} // altruct
