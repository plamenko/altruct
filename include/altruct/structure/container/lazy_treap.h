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
struct lazy_treap_iterator : public std::iterator<std::bidirectional_iterator_tag, T> {
    typedef std::iterator<std::bidirectional_iterator_tag, T> it_t;
    typedef typename it_t::difference_type difference_type;
    IT it;
    lazy_treap_iterator() { }
    lazy_treap_iterator(IT it) : it(it) { }
    int count() const { return it.count(); }
    int size() const { return it.size(); }
    T& operator * () { return *it; }
    T* operator -> () { return &*it; }
    bool operator == (const lazy_treap_iterator& rhs) const { return it == rhs.it; }
    bool operator != (const lazy_treap_iterator& rhs) const { return it != rhs.it; }
    lazy_treap_iterator& operator--() { --it; return *this; }
    lazy_treap_iterator operator--(int) { return it--; }
    lazy_treap_iterator& operator++() { ++it; return *this; }
    lazy_treap_iterator operator++(int) { return it++; }

    // `pos` and `add` are not suitable when duplicate mode is set to COUNT
    int pos() const { return it.pos(); }
    lazy_treap_iterator add(difference_type off) { return it.add(off); }

    lazy_treap_iterator parent() const { return it.parent(); }
    lazy_treap_iterator right() const { return it.right(); }
    lazy_treap_iterator left() const { return it.left(); }
};

/**
 * Bidirectional const iterator.
 */
template<typename T, typename CIT, typename IT>
struct lazy_treap_const_iterator : public std::iterator<std::bidirectional_iterator_tag, const T> {
    typedef std::iterator<std::bidirectional_iterator_tag, const T> it_t;
    typedef typename it_t::difference_type difference_type;
    CIT it;
    lazy_treap_const_iterator() { }
    lazy_treap_const_iterator(CIT it) : it(it) { }
    lazy_treap_const_iterator(IT it) : it(it) { }
    lazy_treap_const_iterator(lazy_treap_iterator<T, IT> it) : it(it.it) { }
    int count() const { return it.count(); }
    int size() const { return it.size(); }
    const T& operator * () const { return *it; }
    const T* operator -> () const { return &*it; }
    bool operator == (const lazy_treap_const_iterator& rhs) const { return it == rhs.it; }
    bool operator != (const lazy_treap_const_iterator& rhs) const { return it != rhs.it; }
    lazy_treap_const_iterator& operator--() { --it; return *this; }
    lazy_treap_const_iterator operator--(int) { return it--; }
    lazy_treap_const_iterator& operator++() { ++it; return *this; }
    lazy_treap_const_iterator operator++(int) { return it++; }

    // `pos` and `add` are not suitable when duplicate mode is set to COUNT
    int pos() const { return it.pos(); }
    lazy_treap_const_iterator add(difference_type off) { return it.add(off); }

    lazy_treap_const_iterator parent() const { return it.parent(); }
    lazy_treap_const_iterator right() const { return it.right(); }
    lazy_treap_const_iterator left() const { return it.left(); }
};

/**
 * Treap with lazy propagation.
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
 * param T - element type
 * param f_up - associative functor for upward propagation; commutativity is not required.
 * param f_down - associative functor for lazy downward propagation; commutativity is not required.
 * param CMP - comparison functor type
 * param RAND - random generator type
 * param ALLOC - allocator type
 */
template<
    typename T,
    typename CMP = std::less<T>,
    typename F_UP = std::function<void(T& parent, const T& left, const T& right)>,
    typename F_DOWN = std::function<void(T& parent, T& left, T& right)>,
    typename RAND = std::function<int()>,
    typename ALLOC = std::allocator<bst_node<T>>>
class lazy_treap {
protected:
    typedef binary_search_tree<T, T, bst_duplicate_handling::STORE, CMP, ALLOC> bst_t;
    bst_t tree;
    RAND rnd;
    F_UP f_up;     // for up-propagation on update
    F_DOWN f_down; // for down-propagation for lazy updating
    T id, dummy;

public:
    typedef T key_type;
    typedef T value_type;
    typedef lazy_treap_iterator<T, typename bst_t::iterator> iterator;
    typedef lazy_treap_const_iterator<T, typename bst_t::const_iterator, typename bst_t::iterator> const_iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    lazy_treap(const F_UP& f_up, const F_DOWN& f_down, const CMP& cmp = CMP(), const RAND& rnd = rand, const ALLOC& alloc = ALLOC()) :
        tree(cmp, alloc), rnd(rnd), f_up(f_up), f_down(f_down) {
    }

    template<typename It>
    lazy_treap(It begin, It end, const F_UP& f_up, const F_DOWN& f_down, const CMP& cmp = CMP(), const RAND& rnd = rand, const ALLOC& alloc = ALLOC()) :
        lazy_treap(f_up, f_down, cmp, rnd, alloc) {
        for (It it = begin; it != end; ++it) {
            insert(*it);
        }
    }

    lazy_treap(lazy_treap&& rhs) :
        lazy_treap(rhs.f_up, rhs.f_down) {
        swap(rhs);
    }

    lazy_treap(const lazy_treap& rhs) :
        tree(rhs.tree),
        rnd(rhs.rnd),
        f_up(rhs.f_up),
        f_down(rhs.f_down) {
    }

    lazy_treap& operator=(lazy_treap&& rhs) {
        swap(rhs);
        return *this;
    }

    lazy_treap& operator=(const lazy_treap& rhs) {
        tree = rhs.tree;
        rnd = rhs.rnd;
        f_up = rhs.f_up;
        f_down = rhs.f_down;
        return *this;
    }

    void swap(lazy_treap& rhs) {
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
    bool operator == (const lazy_treap& rhs) const { return tree == rhs.tree; }
    bool operator < (const lazy_treap& rhs) const { return tree < rhs.tree; }
    bool operator != (const lazy_treap& rhs) const { return !(*this == rhs); }
    bool operator >(const lazy_treap& rhs) const { return (rhs < *this); }
    bool operator <= (const lazy_treap& rhs) const { return !(rhs < *this); }
    bool operator >= (const lazy_treap& rhs) const { return !(*this < rhs); }

public: // query & update
    //int count_less_or_equal(const key_type& key) const {
    //    return tree.count_less_or_equal(key);
    //}

    //int count_less(const key_type& key) const {
    //    return tree.count_less(key);
    //}

    //int count(const key_type& key) const {
    //    return tree.count(key);
    //}

    //const_iterator find_kth(int k) const {
    //    return tree.find_kth(k);
    //}

    //iterator find_kth(int k) {
    //    return tree.find_kth(k);
    //}

    // TODO: propagate down on iterator operations?
    //   this can perhaps be achieved by having enter_child hooks in bst?
    //   alternatively, we could do that on dereferencing iterator?
    //   what is the perf hit?

    iterator find(const key_type& key) {
        iterator res = end();
        for (iterator it = root(); it != end();) {
            propagate_down(it);
            if (tree.compare(*it, key)) {
                it = it.right();
            } else if (tree.compare(key, *it)) {
                it = it.left();
            } else {
                res = it;
                it = it.left();
            }
        }
        return res;
    }

    iterator lower_bound(const key_type& key) {
        iterator res = end();
        for (iterator it = root(); it != end();) {
            propagate_down(it);
            if (tree.compare(*it, key)) {
                it = it.right();
            } else {
                res = it;
                it = it.left();
            }
        }
        return res;
    }

    iterator upper_bound(const key_type& key) {
        iterator res = end();
        for (iterator it = root(); it != end();) {
            propagate_down(it);
            if (tree.compare(key, *it)) {
                res = it;
                it = it.left();
            } else {
                it = it.right();
            }
        }
        return res;
    }

    std::pair<iterator, iterator> equal_range(const key_type& key) {
        return{ lower_bound(key), upper_bound(key) };
    }

    iterator insert(const value_type& val, int cnt = 1) {
        upper_bound(val); // to propagate down
        return propagate_up(retrace_up(tree.insert(val, cnt)));
    }

    iterator insert_before(const_iterator it, const value_type& val, int cnt = 1) {
        return propagate_up(retrace_up(tree.insert_before(it.it, val, cnt)));
    }

    iterator erase(const key_type& key, int cnt = std::numeric_limits<int>::max()) {
        return erase(lower_bound(key), upper_bound(key));
    }

    iterator erase(const_iterator b, const_iterator e, int cnt = std::numeric_limits<int>::max()) {
        iterator res;
        while (b != e) res = erase(b++, cnt);
        return res;
    }

    iterator erase(const_iterator it, int cnt = std::numeric_limits<int>::max()) {
        return propagate_up(tree.erase(retrace_down(it).it, cnt));
    }

    T get(const_iterator _b, const_iterator _e) {
        auto b = remove_const(_b), e = remove_const(_e);
        if (b == e--) return id;
        auto a = remove_const(lowest_common_ancestor(b, e));
        T rb = id;
        for (bool was_left = true; b != a; b = b.parent()) {
            if (was_left) dummy = rb, rb = *b, f_up(rb, dummy, _val_or_id(b.right())); // aggregate b and its right child
            was_left = (b == b.parent().left());
        }
        T re = id;
        for (bool was_right = true; e != a; e = e.parent()) {
            if (was_right) dummy = re, re = *e, f_up(re, _val_or_id(e.left(), dummy)); // aggregate e and its left child
            was_right = (e == e.parent().right());
        }
        T r = *a; f_up(r, rb, re); // aggregate a
        return r;
    }

    template<typename F_UPDATE>
    void update(const_iterator _b, const_iterator _e, const F_UPDATE& f) {
        auto b = remove_const(_b), e = remove_const(_e);
        if (b == e--) return;
        auto a = remove_const(lowest_common_ancestor(b, e));
        for (bool was_left = true; b != a; b = b.parent()) {
            if (was_left) f(*b), f_down(*b, dummy, _val_or_dummy(b.right())); // update b and its right child
            was_left = (b == b.parent().left());
        }
        for (bool was_right = true; e != a; e = e.parent()) {
            if (was_right) f(*e), f_down(*e, _val_or_dummy(e.left()), dummy); // update e and its left child
            was_right = (e == e.parent().right());
        }
        f(*a), f_down(*a, dummy, dummy); // update a
    }

    void propagate_down_to(const_iterator it) {
        if (it == end()) return;
        propagate_down_to(it.parent());
        propagate_down(it);
    }

public: // const casting logic
    iterator remove_const(const_iterator it) {
        // the iterator itself is const, but this method is not so this is safe
        return tree.remove_const(it.it);
    }

    const lazy_treap* const_this() const {
        // it is always safe to cast from non-const to const
        return static_cast<const lazy_treap*>(this);
    }

protected: // lowest_common_ancestor
    const_iterator lowest_common_ancestor(const_iterator b, const_iterator e) {
        int db = depth(b), de = depth(e);
        for (; db > de; db--) b = b.parent();
        for (; de > db; de--) e = e.parent();
        while (b != e) b = b.parent(), e = e.parent();
        return b;
    }

    int depth(const_iterator it) {
        int d = 0;
        for (; it != end(); it = it.parent()) d++;
        return d;
    }

protected: // balancing
    iterator retrace_up(const_iterator _it) {
        auto it = remove_const(_it).it;
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

    iterator retrace_down(const_iterator _it) {
        auto it = remove_const(_it).it;
        if (it == tree.end()) return it;
        propagate_down(it);
        while (it.left() != tree.end() && it.right() != tree.end()) {
            if (it.left().balance() < it.right().balance()) {
                propagate_down(it.left());
                tree.rotate_right(it);
            } else {
                propagate_down(it.right());
                tree.rotate_left(it);
            }
        }
        return it;
    }

public: // aggregation
    iterator propagate_up(const_iterator _it) {
        auto it = remove_const(_it);
        for (auto it2 = it; it2 != tree.end(); it2 = it2.parent()) {
            f_up(*it2, _val_or_id(it2.left()), _val_or_id(it2.right()));
        }
        return it;
    }

    iterator propagate_down(const_iterator _it) {
        auto it = remove_const(_it);
        if (it == tree.end()) return it;
        f_down(*it, _val_or_dummy(it.left()), _val_or_dummy(it.right()));
        return it;
    }

private:
    const T& _val_or_id(iterator it) {
        return (it != tree.end()) ? *it : id;
    }
    
    T& _val_or_dummy(iterator it) {
        return (it != tree.end()) ? *it : dummy;
    }
};

} // container
} // altruct
