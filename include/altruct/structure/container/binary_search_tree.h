#pragma once

#include <iterator>
#include <type_traits>
#include <algorithm>

namespace altruct {
namespace container {

/**
 * A helper template for extracting key from the value.
 */
template<typename K, typename T>
struct bst_key {
    static const K& of(const T& val);
};

/**
 * A key helper template for sets.
 */
template<typename K>
struct bst_key<K, K> {
    static const K& of(const K& val) {
        return val;
    }
};

/**
 * A key helper template for maps.
 */
template<typename K, typename V>
struct bst_key<K, std::pair<const K, V>> {
    static const K& of(const std::pair<const K, V>& entry) {
        return entry.first;
    }
};

/**
 * Binary-search-tree node structure.
 *
 * Contains a reserved `balance` field to be used by subclasses for balancing.
 * This field may be used as height for AVL tree, color for red-black tree,
 * random number for treap, etc.
 *
 * A node whose parent points to itself is a special `nil` node as described in
 * `binary_search_tree` class.
 */
template<typename T>
struct bst_node {
    bst_node* parent;
    bst_node* left;
    bst_node* right;
    int balance; // used for balancing
    int size;    // size of a subtree rooted at this node
    T val;

    bst_node(const T& val) : val(val) {}

    bool is_nil() const { return parent == this; }
    int count() const { return is_nil() ? 0 : size - left->size - right->size; }
};

/**
 * Iterator utilities for inorder traversal.
 */
template<typename T>
struct bst_iterator_util {
    typedef bst_node<T>* node_ptr;
    typedef const bst_node<T>* const_node_ptr;

    // inorder previous; for nil this returns the last node
    static const_node_ptr inorder_prev(const_node_ptr ptr) {
        const_node_ptr ch = ptr->left;
        if (ch->is_nil()) {
            while (ch = ptr, ptr = ch->parent, ch != ptr->right);
        } else {
            while (ptr = ch, ch = ptr->right, !ch->is_nil());
        }
        return ptr;
    }
    static node_ptr inorder_prev(node_ptr ptr) {
        return remove_const(inorder_prev(const_node_ptr(ptr)));
    }

    // inorder next; for nil this returns the first node
    static const_node_ptr inorder_next(const_node_ptr ptr) {
        const_node_ptr ch = ptr->right;
        if (ch->is_nil()) {
            while (ch = ptr, ptr = ch->parent, ch != ptr->left);
        } else {
            while (ptr = ch, ch = ptr->left, !ch->is_nil());
        }
        return ptr;
    }
    static node_ptr inorder_next(node_ptr ptr) {
        return remove_const(inorder_next(const_node_ptr(ptr)));
    }

    // inorder add; for nil this returns nil
    static const_node_ptr inorder_add(const_node_ptr ptr, int off) {
        const_node_ptr nil;
        int pos = bst_iterator_util<T>::inorder_pos(ptr, &nil);
        return bst_iterator_util<T>::inorder_kth(nil->left, pos + off);
    }
    static node_ptr inorder_add(node_ptr ptr, int off) {
        return remove_const(inorder_add(const_node_ptr(ptr), off));
    }

    // inorder position; for nil this returns the size
    static int inorder_pos(const_node_ptr ptr, const_node_ptr* out_nil = nullptr) {
        int k = ptr->left->size; bool was_right = false;
        for (; !ptr->is_nil(); ptr = ptr->parent) {
            if (was_right) k += ptr->size - ptr->right->size;
            was_right = ptr->parent->right == ptr;
        }
        if (out_nil) *out_nil = ptr;
        return k;
    }
    static int inorder_pos(node_ptr ptr, node_ptr* out_nil = nullptr) {
        return inorder_pos(const_node_ptr(ptr), (const_node_ptr*)out_nil);
    }

    // node pointer at k-th position within the subtree rooted at ptr
    static const_node_ptr inorder_kth(const_node_ptr ptr, int k) {
        while (!ptr->is_nil()) {
            if (k < ptr->left->size) {
                ptr = ptr->left;
            } else if ((k -= ptr->size - ptr->right->size) >= 0) {
                ptr = ptr->right;
            } else {
                return ptr;
            }
        }
        return ptr;
    }
    static node_ptr inorder_kth(node_ptr ptr, int k) {
        return remove_const(inorder_kth(const_node_ptr(ptr), k));
    }

private:
    static node_ptr remove_const(const_node_ptr ptr) { return const_cast<node_ptr>(ptr); }
};

/**
 * Bidirectional iterator.
 */
template<typename T>
struct bst_iterator : public std::iterator<std::bidirectional_iterator_tag, T> {
    typedef std::iterator<std::bidirectional_iterator_tag, T> it_t;
    typedef typename it_t::difference_type difference_type;
    typedef bst_node<T>* node_ptr;
    node_ptr ptr;
    bst_iterator(node_ptr ptr = nullptr) : ptr(ptr) { }

    int count() const { return ptr->count(); }
    int size() const { return ptr->size; }
    int& balance() { return ptr->balance; }
    T& operator * () { return ptr->val; }
    T* operator -> () { return &ptr->val; }
    bool operator == (const bst_iterator& rhs) const { return ptr == rhs.ptr; }
    bool operator != (const bst_iterator& rhs) const { return ptr != rhs.ptr; }
    bst_iterator& operator--() { ptr = bst_iterator_util<T>::inorder_prev(ptr); return *this; }
    bst_iterator operator--(int) { auto old = *this; --*this; return old; }
    bst_iterator& operator++() { ptr = bst_iterator_util<T>::inorder_next(ptr); return *this; }
    bst_iterator operator++(int) { auto old = *this; ++*this; return old; }

    // `pos` and `add` are not suitable when duplicate mode is set to COUNT
    int pos() const { return bst_iterator_util<T>::inorder_pos(ptr); }
    bst_iterator add(difference_type off) { return bst_iterator_util<T>::inorder_add(ptr, int(off)); }

    bst_iterator parent() const { return ptr->parent; }
    bst_iterator right() const { return ptr->right; }
    bst_iterator left() const { return ptr->left; }
};

/**
 * Bidirectional const iterator.
 */
template<typename T>
struct bst_const_iterator : public std::iterator<std::bidirectional_iterator_tag, const T> {
    typedef std::iterator<std::bidirectional_iterator_tag, const T> it_t;
    typedef typename it_t::difference_type difference_type;
    typedef const bst_node<T>* const_node_ptr;
    const_node_ptr ptr;
    bst_const_iterator(const_node_ptr ptr = nullptr) : ptr(ptr) { }
    bst_const_iterator(bst_iterator<T> it) : ptr(it.ptr) { }

    int count() const { return ptr->count(); }
    int size() const { return ptr->size; }
    int& balance() { return *const_cast<int*>(&ptr->balance); }
    const T& operator * () const { return ptr->val; }
    const T* operator -> () const { return &ptr->val; }
    bool operator == (const bst_const_iterator& rhs) const { return ptr == rhs.ptr; }
    bool operator != (const bst_const_iterator& rhs) const { return ptr != rhs.ptr; }
    bst_const_iterator& operator--() { ptr = bst_iterator_util<T>::inorder_prev(ptr); return *this; }
    bst_const_iterator operator--(int) { auto old = *this; --*this; return old; }
    bst_const_iterator& operator++() { ptr = bst_iterator_util<T>::inorder_next(ptr); return *this; }
    bst_const_iterator operator++(int) { auto old = *this; ++*this; return old; }
    bst_const_iterator& operator+=(difference_type off) { ptr = bst_iterator_util<T>::inorder_add(ptr, int(off)); return *this; }

    // `pos` and `add` are not suitable when duplicate mode is set to COUNT
    int pos() const { return bst_iterator_util<T>::inorder_pos(ptr); }
    bst_const_iterator add(difference_type off) { return bst_iterator_util<T>::inorder_add(ptr, int(off)); }

    bst_const_iterator parent() const { return ptr->parent; }
    bst_const_iterator right() const { return ptr->right; }
    bst_const_iterator left() const { return ptr->left; }
};

/**
 * Duplicate handling mode.
 */
namespace bst_duplicate_handling {
    enum {
        IGNORE, // duplicate keys are stored only once, with count being 1
        COUNT,  // duplicate keys are stored once, with count being tracked
        STORE,  // duplicate keys are stored separately in the insertion order
    };
}

/**
 * Binary search tree.
 *
 * Note, no balancing is being performed. This is left for subclasses.
 *
 * Space complexity: `O(n)`.
 * Time complexities:
 *   find:   `O(h)`
 *   insert: `O(h)`
 *   erase:  `O(h)`
 * Where `h` is height that may be proportional to the tree size `n`.
 * Subclasses may choose to perform balancing in order to improve that.
 *
 * param K   - key type
 * param T   - value type (for maps this is going to be std::pair<const K, V>)
 * param DUP - duplicate handling mode
 * param CMP - comparison functor type
 * param ALLOC - allocator type
 */
template<typename K, typename T = K, int DUP = bst_duplicate_handling::IGNORE, typename CMP = std::less<K>, typename ALLOC = std::allocator<bst_node<T>>>
class binary_search_tree {
public: // public types
    typedef K key_type;
    typedef T value_type;
    //typedef typename std::conditional<
    //    std::is_same<K, T>::value,
    //    bst_const_iterator<T>,
    //    bst_iterator<T> >::type iterator;
    typedef bst_iterator<T> iterator;
    typedef bst_const_iterator<T> const_iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

protected: // member variables
    typedef bst_node<T>* node_ptr;
    typedef const bst_node<T>* const_node_ptr;

    CMP cmp;
    ALLOC alloc;

    // `nil` node is conceptually equivalent to a null pointer,
    // but it is instead an actual node with some handy properties:
    // * `nil->parent` is always the `nil` node itself, and
    //   is the only such node which allows for simple check.
    // * `nil->left` and `nil->right` always point to the root
    //   which further implies the following handy properties:
    // * There is no need to explicitly keep a pointer to the root
    //   as one can simply do `nil->left` or `nil->right`.
    // * There is no need for book keeping when the root changes
    //   as one can simply link parent and child as in `make_link`.
    // * `nil` conceptually closes a circle by being the node after
    //   the last node and before the first node. This then implies:
    // *   end == nil
    // *   last == inorder_prev(nil)
    // *   begin == inorder_next(nil)
    // *   nil == inorder_prev(begin)
    // *   nil == inorder_next(last)
    node_ptr nil;

    node_ptr root_ptr() {
        return nil->left;
    }

    const_node_ptr root_ptr() const {
        return nil->left;
    }

public: // constructor & size
    ~binary_search_tree() {
        _destroy();
    }

    binary_search_tree(const CMP& cmp = CMP(), const ALLOC& alloc = ALLOC()) :
        cmp(cmp), alloc(alloc) {
        _init();
    }

    template<typename It>
    binary_search_tree(It begin, It end, const CMP& cmp = CMP(), const ALLOC& alloc = ALLOC()) :
        binary_search_tree(cmp, alloc) {
        for (It it = begin; it != end; ++it) {
            insert(*it);
        }
    }

    binary_search_tree(std::initializer_list<T> list) :
        binary_search_tree(list.begin(), list.end()) {
    }

    binary_search_tree(binary_search_tree&& rhs) :
        binary_search_tree() {
        swap(rhs);
    }

    binary_search_tree(const binary_search_tree& rhs) :
        binary_search_tree(rhs.cmp, rhs.alloc) {
        make_link(nil, clone_subtree(rhs.root_ptr()), true);
    }

    binary_search_tree& operator=(binary_search_tree&& rhs) {
        swap(rhs);
        return *this;
    }

    binary_search_tree& operator=(const binary_search_tree& rhs) {
        binary_search_tree rhs_copy(rhs);
        swap(rhs_copy);
        return *this;
    }

    void swap(binary_search_tree& rhs) {
        std::swap(cmp, rhs.cmp);
        std::swap(alloc, rhs.alloc);
        std::swap(nil, rhs.nil);
    }

    void clear() {
        _free(root_ptr());
        nil->left = nil;
        nil->right = nil;
    }

    bool empty() const {
        return root_ptr()->size == 0;
    }

    int size() const {
        return root_ptr()->size;
    }

public: // iterators
    iterator root() { return root_ptr(); }
    const_iterator root() const { return root_ptr(); }
    const_iterator croot() const { return root_ptr(); }
    iterator begin() { return ++end(); }
    const_iterator begin() const { return ++end(); }
    const_iterator cbegin() const { return ++cend(); }
    iterator end() { return nil; }
    const_iterator end() const { return nil; }
    const_iterator cend() const { return nil; }
    reverse_iterator rbegin() { return reverse_iterator(end()); }
    const_reverse_iterator rbegin() const { return const_reverse_iterator(end()); }
    const_reverse_iterator crbegin() const { return const_reverse_iterator(cend()); }
    reverse_iterator rend() { return reverse_iterator(begin()); }
    const_reverse_iterator rend() const { return const_reverse_iterator(begin()); }
    const_reverse_iterator crend() const { return const_reverse_iterator(cbegin()); }

public: // relational operators
    bool compare(const T& v1, const T& v2) const { return cmp(_key(v1), _key(v2)); }
    bool operator == (const binary_search_tree& rhs) const {
        auto b1 = cbegin(), e1 = cend(), b2 = rhs.cbegin(), e2 = rhs.cend();
        for (; b1 != e1 && b2 != e2; ++b1, ++b2) {
            if (compare(*b1, *b2)) return false;
            if (compare(*b2, *b1)) return false;
        }
        return b1 == e1 && b2 == e2;
    }
    bool operator < (const binary_search_tree& rhs) const {
        auto b1 = cbegin(), e1 = cend(), b2 = rhs.cbegin(), e2 = rhs.cend();
        for (; b1 != e1 && b2 != e2; ++b1, ++b2) {
            if (compare(*b1, *b2)) return true;
            if (compare(*b2, *b1)) return false;
        }
        return b1 == e1 && b2 != e2;
    }
    bool operator != (const binary_search_tree& rhs) const { return !(*this == rhs); }
    bool operator >  (const binary_search_tree& rhs) const { return (rhs < *this); }
    bool operator <= (const binary_search_tree& rhs) const { return !(rhs < *this); }
    bool operator >= (const binary_search_tree& rhs) const { return !(*this < rhs); }

public: // query & update
    int count_less_or_equal(const K& key) const {
        int k = 0;
        for (const_node_ptr ptr = root_ptr(); !ptr->is_nil();) {
            if (cmp(key, _key(ptr->val))) {
                ptr = ptr->left;
            } else {
                k += ptr->size - ptr->right->size;
                ptr = ptr->right;
            }
        }
        return k;
    }

    int count_less(const K& key) const {
        int k = 0;
        for (const_node_ptr ptr = root_ptr(); !ptr->is_nil();) {
            if (cmp(_key(ptr->val), key)) {
                k += ptr->size - ptr->right->size;
                ptr = ptr->right;
            } else {
                ptr = ptr->left;
            }
        }
        return k;
    }

    int count(const K& key) const {
        if (DUP == bst_duplicate_handling::STORE) {
            return count_less_or_equal(key) - count_less(key);
        } else {
            return find(key).count();
        }
    }

    const_iterator find_kth(int k) const {
        return bst_iterator_util<T>::inorder_kth(root_ptr(), k);
    }

    iterator find_kth(int k) {
        return remove_const(const_this()->find_kth(k));
    }

    const_iterator find(const K& key) const {
        const_node_ptr res = nil;
        for (const_node_ptr ptr = root_ptr(); !ptr->is_nil();) {
            if (cmp(_key(ptr->val), key)) {
                ptr = ptr->right;
            } else if (cmp(key, _key(ptr->val))) {
                ptr = ptr->left;
            } else {
                if (DUP == bst_duplicate_handling::STORE) {
                    res = ptr;
                    ptr = ptr->left;
                } else {
                    return ptr;
                }
            }
        }
        return res;
    }

    iterator find(const K& key) {
        return remove_const(const_this()->find(key));
    }

    const_iterator lower_bound(const K& key) const {
        const_node_ptr res = nil;
        for (const_node_ptr ptr = root_ptr(); !ptr->is_nil();) {
            if (cmp(_key(ptr->val), key)) {
                ptr = ptr->right;
            } else {
                res = ptr;
                ptr = ptr->left;
            }
        }
        return res;
    }

    iterator lower_bound(const K& key) {
        return remove_const(const_this()->lower_bound(key));
    }

    const_iterator upper_bound(const K& key) const {
        const_node_ptr res = nil;
        for (const_node_ptr ptr = root_ptr(); !ptr->is_nil();) {
            if (cmp(key, _key(ptr->val))) {
                res = ptr;
                ptr = ptr->left;
            } else {
                ptr = ptr->right;
            }
        }
        return res;
    }

    iterator upper_bound(const K& key) {
        return remove_const(const_this()->upper_bound(key));
    }

    std::pair<const_iterator, const_iterator> equal_range(const K& key) const {
        return{ lower_bound(key), upper_bound(key) };
    }

    std::pair<iterator, iterator> equal_range(const K& key) {
        return{ lower_bound(key), upper_bound(key) };
    }

    iterator insert(const T& val, int cnt = 1) {
        node_ptr ptr, par = nil; bool go_left = true;
        for (ptr = root_ptr(); !ptr->is_nil();) {
            par = ptr;
            if (cmp(_key(val), _key(ptr->val))) {
                go_left = true;
                ptr = ptr->left;
            } else if (cmp(_key(ptr->val), _key(val))) {
                go_left = false;
                ptr = ptr->right;
            } else {
                if (DUP == bst_duplicate_handling::STORE) {
                    go_left = false;
                    ptr = ptr->right;
                } else if (DUP == bst_duplicate_handling::COUNT) {
                    par = ptr->parent;
                    break;
                } else {
                    return ptr;
                }

            }
        }
        return insert_node(par, go_left, val, cnt);
        // Note for balancing:
        // ptr needs to be retraced up after insert
    }

    // important: this call is unchecked and the sort order may be violated
    iterator insert_before(const_iterator it, const T& val, int cnt = 1) {
        bool go_left = it.ptr->left->is_nil(); if (!go_left) --it;
        return insert_node(remove_const(it).ptr, go_left, val, cnt);
    }

    iterator erase(const K& key, int cnt = std::numeric_limits<int>::max()) {
        if (DUP == bst_duplicate_handling::STORE) {
            return erase(lower_bound(key), upper_bound(key));
        } else {
            return erase(find(key), cnt);
        }
    }

    iterator erase(const_iterator b, const_iterator e, int cnt = std::numeric_limits<int>::max()) {
        while (b != e) erase(b++, cnt);
        return remove_const(b);
    }

    iterator erase(const_iterator it, int cnt = std::numeric_limits<int>::max()) {
        return erase_node(remove_const(it).ptr, cnt);
        // Note for balancing:
        // for avl like trees: `par` needs to be retraced up after erase
        // for treap like trees: `ptr` needs to be retraced down before erase
    }

public: // const casting logic
    iterator remove_const(const_iterator it) {
        // the iterator itself is const, but this method is not so this is safe
        return const_cast<node_ptr>(it.ptr);
    }

    const binary_search_tree* const_this() const {
        // it is always safe to cast from non-const to const
        return static_cast<const binary_search_tree*>(this);
    }

protected: // insert & erase
    node_ptr insert_node(node_ptr par, bool go_left, const T& val, int cnt) {
        node_ptr ptr = go_left ? par->left : par->right;
        if (ptr->is_nil()) {
            ptr = _buy(val);
            ptr->parent = par;
            ptr->left = nil;
            ptr->right = nil;
            ptr->size = 0;
            make_link(par, ptr, go_left);
        }
        if (DUP != bst_duplicate_handling::COUNT) {
            cnt = 1;
        }
        propagate_size(ptr, nil, cnt);
        return ptr;
    }

    node_ptr erase_node(node_ptr ptr, int cnt) {
        if (ptr->is_nil()) return nil;
        if (DUP != bst_duplicate_handling::COUNT) {
            cnt = 1;
        }
        if (cnt < ptr->count()) {
            propagate_size(ptr, nil, -cnt);
            return ptr;
        }
        cnt = ptr->count();
        // physically erase from the tree
        if (!ptr->left->is_nil() && !ptr->right->is_nil()) {
            swap_with_descendant(ptr, bst_iterator_util<T>::inorder_next(ptr));
        }
        // there can be at most one child now
        node_ptr ch = (!ptr->left->is_nil()) ? ptr->left : ptr->right;
        node_ptr par = ptr->parent;
        make_link(par, ch, ptr);
        propagate_size(par, nil, -cnt);
        _free(ptr);
        return par;
    }

    node_ptr clone_subtree(const_node_ptr src) {
        // note this method gets called from copy-constructor;
        // this means that src may come from another tree
        if (src->is_nil()) return nil;
        auto ptr = _buy(src->val);
        ptr->parent = nil;
        make_link(ptr, clone_subtree(src->left), true);
        make_link(ptr, clone_subtree(src->right), false);
        ptr->balance = src->balance;
        ptr->size = src->size;
        return ptr;
    }

public: // pointer rewiring logic
    iterator rotate_left(const_iterator it) {
        node_ptr ptr = remove_const(it).ptr;
        node_ptr ch = ptr->right;
        int sz = ptr->size;
        make_link(ptr->parent, ch, ptr);
        make_link(ptr, ch->left, ch);
        ptr->size += ch->left->size - ch->size;
        ch->left = ptr, ptr->parent = ch;
        ch->size = sz;
        return ch;
    }

    iterator rotate_right(const_iterator it) {
        node_ptr ptr = remove_const(it).ptr;
        node_ptr ch = ptr->left;
        int sz = ptr->size;
        make_link(ptr->parent, ch, ptr);
        make_link(ptr, ch->right, ch);
        ptr->size += ch->right->size - ch->size;
        ch->right = ptr, ptr->parent = ch;
        ch->size = sz;
        return ch;
    }

protected:
    static void swap_with_descendant(node_ptr ptr, node_ptr des) {
        std::swap(ptr->parent, des->parent);
        if (ptr->parent == ptr) ptr->parent = des;
        make_link(ptr->parent, ptr, des);
        make_link(des->parent, des, ptr);
        std::swap(ptr->left, des->left);
        if (des->left == des) des->left = ptr;
        if (!ptr->left->is_nil()) ptr->left->parent = ptr;
        if (!des->left->is_nil()) des->left->parent = des;
        std::swap(ptr->right, des->right);
        if (des->right == des) des->right = ptr;
        if (!ptr->right->is_nil()) ptr->right->parent = ptr;
        if (!des->right->is_nil()) des->right->parent = des;
        std::swap(ptr->balance, des->balance);
        std::swap(ptr->size, des->size);
        propagate_size(ptr, des, des->count() - ptr->count());
    }

    static void make_link(node_ptr par, node_ptr ch, bool go_left) {
        if (par->is_nil() || go_left) par->left = ch;
        if (par->is_nil() || !go_left) par->right = ch;
        if (!ch->is_nil()) ch->parent = par;
    }

    static void make_link(node_ptr par, node_ptr ch, node_ptr old_ch) {
        // knowing what the old child was allows for simpler checks
        if (par->right == old_ch) par->right = ch;
        if (par->left == old_ch) par->left = ch;
        if (!ch->is_nil()) ch->parent = par;
    }

    static void propagate_size(node_ptr ptr, node_ptr end, int cnt) {
        if (cnt == 0) return;
        for (node_ptr tmp = ptr; tmp != end; tmp = tmp->parent) {
            tmp->size += cnt;
        }
    }

private: // allocation logic
    static const K& _key(const T& val) {
        return bst_key<K, T>::of(val);
    }

    void _init() {
        nil = alloc.allocate(1); // no construct
        nil->parent = nil; // self-loop indicates that this is the nil node
        nil->left = nil;   // always points to the root (which is now nil)
        nil->right = nil;  // always points to the root (which is now nil)
        nil->balance = 0;
        nil->size = 0;
    }

    void _destroy() {
        _free_all(root_ptr());
        alloc.deallocate(nil, 1); // no destruct
    }

    node_ptr _buy(const T& val) {
        node_ptr t = alloc.allocate(1);
        alloc.construct(t, val);
        return t;
    }

    void _free(node_ptr t) {
        alloc.destroy(t);
        alloc.deallocate(t, 1);
    }

    void _free_all(node_ptr t) {
        if (t->is_nil()) return;
        _free_all(t->left);
        _free_all(t->right);
        _free(t);
    }
};

} // container
} // altruct
