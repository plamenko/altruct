#pragma once

#include "structure/container/treap.h"

#include <cstdlib>
#include <functional>

namespace altruct {
namespace container {

/**
 * Random access iterator.
 */
template<typename T>
struct rope_iterator : public std::iterator<std::random_access_iterator_tag, T> {
    typedef bst_node<T>* node_ptr;
    typedef std::iterator<std::random_access_iterator_tag, T> it_t;
    typedef typename it_t::difference_type difference_type;
    node_ptr ptr;
    rope_iterator(node_ptr ptr = nullptr) : ptr(ptr) { }
    int count() const { return ptr->count(); }
    int pos() const { return bst_iterator_util<T>::inorder_pos(ptr); }
    T& operator * () { return ptr->val; }
    T* operator -> () { return &ptr->val; }
    bool operator == (const rope_iterator& rhs) const { return ptr == rhs.ptr; }
    bool operator != (const rope_iterator& rhs) const { return ptr != rhs.ptr; }
    bool operator < (const rope_iterator& rhs) const { return (*this - rhs) < 0; }
    bool operator > (const rope_iterator& rhs) const { return (*this - rhs) > 0; }
    bool operator <= (const rope_iterator& rhs) const { return (*this - rhs) <= 0; }
    bool operator >= (const rope_iterator& rhs) const { return (*this - rhs) >= 0; }
    rope_iterator& operator--() { ptr = bst_iterator_util<T>::inorder_prev(ptr); return *this; }
    rope_iterator operator--(int) { auto old = *this; --*this; return old; }
    rope_iterator& operator++() { ptr = bst_iterator_util<T>::inorder_next(ptr); return *this; }
    rope_iterator operator++(int) { auto old = *this; ++*this; return old; }
    rope_iterator& operator+=(difference_type off) { ptr = bst_iterator_util<T>::inorder_add(ptr, int(off)); return *this; }
    rope_iterator operator+(difference_type off) const { rope_iterator it = *this; return it += off; }
    rope_iterator& operator-=(difference_type off) { return *this += -off; }
    rope_iterator operator-(difference_type off) const { rope_iterator it = *this; return it -= off; }
    difference_type operator-(const rope_iterator& rhs) const { return pos() - rhs.pos(); }
    T& operator[](difference_type off) const { return *(*this + off); }
};

/**
 * Random access const iterator.
 */
template<typename T>
struct rope_const_iterator : public std::iterator<std::random_access_iterator_tag, const T> {
    typedef const bst_node<T>* const_node_ptr;
    typedef std::iterator<std::random_access_iterator_tag, T> it_t;
    typedef typename it_t::difference_type difference_type;
    const_node_ptr ptr;
    rope_const_iterator(const_node_ptr ptr = nullptr) : ptr(ptr) { }
    rope_const_iterator(rope_iterator<T> it) : ptr(it.ptr) { }
    int count() const { return ptr->count(); }
    int pos() const { return bst_iterator_util<T>::inorder_pos(ptr); }
    const T& operator * () const { return ptr->val; }
    const T* operator -> () const { return &ptr->val; }
    bool operator == (const rope_const_iterator& rhs) const { return ptr == rhs.ptr; }
    bool operator != (const rope_const_iterator& rhs) const { return ptr != rhs.ptr; }
    bool operator < (const rope_const_iterator& rhs) const { return (*this - rhs) < 0; }
    bool operator > (const rope_const_iterator& rhs) const { return (*this - rhs) > 0; }
    bool operator <= (const rope_const_iterator& rhs) const { return (*this - rhs) <= 0; }
    bool operator >= (const rope_const_iterator& rhs) const { return (*this - rhs) >= 0; }
    rope_const_iterator& operator--() { ptr = bst_iterator_util<T>::inorder_prev(ptr); return *this; }
    rope_const_iterator operator--(int) { auto old = *this; --*this; return old; }
    rope_const_iterator& operator++() { ptr = bst_iterator_util<T>::inorder_next(ptr); return *this; }
    rope_const_iterator operator++(int) { auto old = *this; ++*this; return old; }
    rope_const_iterator& operator+=(difference_type off) { ptr = bst_iterator_util<T>::inorder_add(ptr, int(off)); return *this; }
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
 * param K   - key type
 * param T   - value type (for maps this is going to be bst_entry<K, V>)
 * param DUP - duplicate handling mode
 * param CMP - comparison functor type
 * param ALLOC - allocator type
 */
template<typename T, typename RAND = std::function<int()>, typename ALLOC = std::allocator<bst_node<T>>>
class rope : protected treap<T, T, bst_duplicate_handling::COUNT, std::less<T>, RAND, ALLOC> {
protected:
    typedef treap<T, T, bst_duplicate_handling::COUNT, std::less<T>, RAND, ALLOC> treap_t;
    typedef typename treap_t::node_ptr node_ptr;
    typedef typename treap_t::const_node_ptr const_node_ptr;
public:
    typedef T value_type;
    typedef rope_iterator<T> iterator;
    typedef rope_const_iterator<T> const_iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    rope(const RAND& rnd = rand, const ALLOC& alloc = ALLOC()) :
        treap_t(std::less<T>(), rnd, alloc) {
    }

    template<typename It>
    rope(It begin, It end, const RAND& rnd = rand, const ALLOC& alloc = ALLOC()) :
        treap_t(std::less<T>(), rnd, alloc) {
        for (It it = begin; it != end; ++it) {
            push_back(*it);
        }
    }

    rope(std::initializer_list<T> list) :
        rope(list.begin(), list.end()) {}

    rope(rope&& rhs) :
        rope() {
        swap(rhs);
    }

    rope(const rope& rhs) :
        rope(rhs.cbegin(), rhs.cend(), rhs.rnd, rhs.alloc) {
    }

    rope& operator=(rope&& rhs) {
        swap(rhs);
        return *this;
    }

    rope& operator=(const rope& rhs) {
        rope rhs_copy(rhs);
        swap(rhs_copy);
        return *this;
    }

    void swap(rope& rhs) {
        treap_t::swap(rhs);
    }

    void clear() {
        treap_t::clear();
    }

    bool empty() const {
        return treap_t::empty();
    }

    int size() const {
        return treap_t::size();
    }

public: // iterators
    iterator begin() { return ++end(); }
    iterator end() { return treap_t::nil; }
    const_iterator cbegin() const { return ++cend(); }
    const_iterator cend() const { return treap_t::nil; }
    reverse_iterator rbegin() { return reverse_iterator(end()); }
    reverse_iterator rend() { return reverse_iterator(begin()); }
    const_reverse_iterator crbegin() const { return const_reverse_iterator(cend()); }
    const_reverse_iterator crend() const { return const_reverse_iterator(cbegin()); }

public: // relational operators
    bool operator == (const rope& rhs) const {
        return treap_t::equal(cbegin(), cend(), rhs.cbegin(), rhs.cend());
    }
    bool operator < (const rope& rhs) const {
        return std::lexicographical_compare(cbegin(), cend(), rhs.cbegin(), rhs.cend());
    }
    bool operator != (const rope& rhs) const { return !(*this == rhs); }
    bool operator >  (const rope& rhs) const { return (rhs < *this); }
    bool operator <= (const rope& rhs) const { return !(rhs < *this); }
    bool operator >= (const rope& rhs) const { return !(*this < rhs); }

public: // query & update
    const_iterator find_kth(int k) const {
        return treap_t::find_kth(k).ptr;
    }

    iterator find_kth(int k) {
        return treap_t::remove_const(treap_t::find_kth(k));
    }

    void push_back(const T& val) {
        insert(size(), val);
    }

    void insert(int pos, const T& val) {
        insert(find_kth(pos), val);
    }

    void insert(const_iterator it, const T& val) {
        bool go_left = it.ptr->left->is_nil(); if (!go_left) --it;
        treap_t::retrace_up(treap_t::insert_node(treap_t::remove_const(it.ptr), go_left, val, 1));
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
        treap_t::erase_node(treap_t::remove_const(treap_t::retrace_down(it.ptr)), 1);
    }

    const T& at(int pos) const {
        return find_kth(pos).ptr->val;
    }

    T& at(int pos) {
        return treap_t::remove_const(find_kth(pos).ptr)->val;
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
