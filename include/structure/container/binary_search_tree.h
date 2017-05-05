#pragma once

#include <iterator>
#include <type_traits>
#include <algorithm>

namespace altruct {
namespace container {

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

private:
    static node_ptr remove_const(const_node_ptr ptr) { return const_cast<node_ptr>(ptr); }
};

/**
 * Bidirectional iterator.
 */
template<typename T>
struct bst_iterator : public std::iterator<std::bidirectional_iterator_tag, T> {
    typedef bst_node<T>* node_ptr;
    node_ptr ptr;
    bst_iterator(node_ptr ptr = nullptr) : ptr(ptr) { }
    int count() const { return ptr->count(); }
    T& operator * () { return ptr->val; }
    T* operator -> () { return &ptr->val; }
    bool operator == (const bst_iterator& rhs) const { return ptr == rhs.ptr; }
    bool operator != (const bst_iterator& rhs) const { return ptr != rhs.ptr; }
    bst_iterator& operator--() { ptr = bst_iterator_util<T>::inorder_prev(ptr); return *this; }
    bst_iterator operator--(int) { auto old = *this; --*this; return old; }
    bst_iterator& operator++() { ptr = bst_iterator_util<T>::inorder_next(ptr); return *this; }
    bst_iterator operator++(int) { auto old = *this; ++*this; return old; }
};

/**
 * Bidirectional const iterator.
 */
template<typename T>
struct bst_const_iterator : public std::iterator<std::bidirectional_iterator_tag, const T> {
    typedef const bst_node<T>* const_node_ptr;
    const_node_ptr ptr;
    bst_const_iterator(const_node_ptr ptr = nullptr) : ptr(ptr) { }
    bst_const_iterator(bst_iterator<T> it) : ptr(it.ptr) { }
    int count() const { return ptr->count(); }
    const T& operator * () const { return ptr->val; }
    const T* operator -> () const { return &ptr->val; }
    bool operator == (const bst_const_iterator& rhs) const { return ptr == rhs.ptr; }
    bool operator != (const bst_const_iterator& rhs) const { return ptr != rhs.ptr; }
    bst_const_iterator& operator--() { ptr = bst_iterator_util<T>::inorder_prev(ptr); return *this; }
    bst_const_iterator operator--(int) { auto old = *this; --*this; return old; }
    bst_const_iterator& operator++() { ptr = bst_iterator_util<T>::inorder_next(ptr); return *this; }
    bst_const_iterator operator++(int) { auto old = *this; ++*this; return old; }
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
    static const K& of(const std::pair<const K, V>& val) {
        return val.first;
    }
};

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
 * param T   - value type (for maps this is going to be pair<const K, V>)
 * param DUP - duplicate handling mode
 * param CMP - comparison functor type
 * param ALLOC - allocator type
 */
template<typename K, typename T = K, int DUP = bst_duplicate_handling::IGNORE, typename CMP = std::less<K>, typename ALLOC = allocator<bst_node<T>>>
class binary_search_tree {
public: // public types
    typedef K key_type;
    typedef T value_type;
    typedef typename std::conditional<std::is_same<key_type, value_type>::value,
        typename bst_const_iterator<T>,
        typename bst_iterator<T> > ::type iterator;
    typedef bst_const_iterator<T> const_iterator;

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

    node_ptr root() {
        return nil->left;
    }

    const_node_ptr root() const {
        return nil->left;
    }

public: // constructor & size
    binary_search_tree(const CMP& cmp = std::less<K>(), const ALLOC& alloc = ALLOC()) : cmp(cmp), alloc(alloc) {
        _init();
    }

    template<typename It>
    binary_search_tree(It begin, It end, const CMP& cmp = std::less<K>(), const ALLOC& alloc = ALLOC()) : cmp(cmp), alloc(alloc) {
        _init();
        for (It it = begin; it != end; ++it) {
            insert(*it);
        }
    }

    ~binary_search_tree() {
        _destroy();
    }

    void clear() {
        _free(root());
        nil->left = nil;
        nil->right = nil;
    }

    bool empty() const {
        return root()->size == 0;
    }

    int size() const {
        return root()->size;
    }

    int count(const T& val) const {
        return find(val).count();
    }

public: // iterators
    iterator begin() {
        return ++end();
    }

    const_iterator cbegin() const {
        return ++cend();
    }

    iterator end() {
        return nil;
    }

    const_iterator cend() const {
        return nil;
    }

public: // query & update
    int count_less(const K& key) const {
        int k = 0;
        for (const_node_ptr ptr = root(); ptr != nil;) {
            if (cmp(_key(ptr->val), key)) {
                k += ptr->size - ptr->right->size;
                ptr = ptr->right;
            } else {
                ptr = ptr->left;
            }
        }
        return k;
    }

    const_iterator find_kth(int k) const {
        for (const_node_ptr ptr = root(); ptr != nil;) {
            if (k < ptr->left->size) {
                ptr = ptr->left;
            } else if ((k -= ptr->size - ptr->right->size) >= 0) {
                ptr = ptr->right;
            } else {
                return ptr;
            }
        }
        return nil;
    }

    iterator find_kth(int k) {
        return remove_const(const_this()->find_kth(k));
    }

    const_iterator find(const K& key) const {
        for (const_node_ptr ptr = root(); ptr != nil;) {
            if (cmp(key, _key(ptr->val))) {
                ptr = ptr->left;
            } else if (cmp(_key(ptr->val), key)) {
                ptr = ptr->right;
            } else {
                return ptr;
            }
        }
        return nil;
    }

    iterator find(const K& key) {
        return remove_const(const_this()->find(key));
    }

    const_iterator lower_bound(const K& key) const {
        const_node_ptr res = nil;
        for (const_node_ptr ptr = root(); ptr != nil;) {
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
        for (const_node_ptr ptr = root(); ptr != nil;) {
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

    iterator insert(const T& val, int cnt = 1) {
        node_ptr ptr, par = nil; bool go_left = true;
        for (ptr = root(); ptr != nil;) {
            par = ptr;
            if (go_left = cmp(_key(val), _key(ptr->val))) {
                ptr = ptr->left;
            } else if (cmp(_key(ptr->val), _key(val))) {
                ptr = ptr->right;
            } else {
                if (DUP == bst_duplicate_handling::STORE) {
                    ptr = ptr->right;
                } else if (DUP == bst_duplicate_handling::COUNT) {
                    break;
                } else {
                    return ptr;
                }

            }
        }
        if (ptr == nil) {
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
        // Note for balancing:
        // ptr needs to be retraced up after insert
    }

    iterator erase(const K& key, int cnt = 1) {
        return erase(find(key), cnt);
    }

    iterator erase(const_iterator it, int cnt = 1) {
        node_ptr ptr = remove_const(it);
        if (ptr == nil) return nil;
        if (DUP != bst_duplicate_handling::COUNT) {
            cnt = 1;
        }
        if (cnt < it.count()) {
            propagate_size(ptr, nil, -cnt);
            return ptr;
        }
        cnt = it.count();
        // physically erase from the tree
        if (ptr->left != nil && ptr->right != nil) {
            ptr = swap_with_next(ptr);
        }
        // there can be at most one child now
        node_ptr ch = (ptr->left != nil) ? ptr->left : ptr->right;
        node_ptr par = ptr->parent;
        make_link(par, ch, ptr);
        propagate_size(par, nil, -cnt);
        return par;
        // Note for balancing:
        // for avl like trees: `par` needs to be retraced up after erase
        // for treap like trees: `ptr` needs to be retraced down before erase
    }

protected: // const casting logic
    const binary_search_tree* const_this() const {
        // it is always safe to cast from non-const to const
        return static_cast<const binary_search_tree*>(this);
    }

    node_ptr remove_const(const_iterator it) {
        // the iterator itself is const, but this method is not
        return const_cast<node_ptr>(it.ptr);
    }

protected: // pointer rewiring logic
    static node_ptr rotate_left(node_ptr ptr) {
        node_ptr ch = ptr->right;
        int sz = ptr->size;
        make_link(ptr->parent, ch, ptr);
        make_link(ptr, ch->left, ch);
        ptr->size += ch->left->size - ch->size;
        ch->left = ptr, ptr->parent = ch;
        ch->size = sz;
        return ch;
    }

    static node_ptr rotate_right(node_ptr ptr) {
        node_ptr ch = ptr->left;
        int sz = ptr->size;
        make_link(ptr->parent, ch, ptr);
        make_link(ptr, ch->right, ch);
        ptr->size += ch->right->size - ch->size;
        ch->right = ptr, ptr->parent = ch;
        ch->size = sz;
        return ch;
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

    static node_ptr swap_with_next(node_ptr ptr) {
        // we could rewire pointers, but it's simpler to just swap the values
        node_ptr ch = bst_iterator_util<T>::inorder_next(ptr);
        std::swap(ptr->val, ch->val);
        propagate_size(ch, ptr, ptr->count() - ch->count());
        return ch;
    }

    static void propagate_size(node_ptr ptr, node_ptr end, int cnt) {
        if (cnt == 0) return;
        for (node_ptr tmp = ptr; tmp != end; tmp = tmp->parent) {
            tmp->size += cnt;
        }
    }

private: // allocation logic
    void _init() {
        nil = alloc.allocate(1); // no construct
        nil->parent = nil; // self-loop indicates that this is the nil node
        nil->left = nil;   // always points to the root (which is now nil)
        nil->right = nil;  // always points to the root (which is now nil)
        nil->balance = 0;
        nil->size = 0;
    }

    void _destroy() {
        _free_all(root());
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
        if (t == nil) return;
        _free_all(t->left);
        _free_all(t->right);
        _free(t);
    }

    const K& _key(const T& val) const {
        return bst_key<K, T>::of(val);
    }
};

} // container
} // altruct
