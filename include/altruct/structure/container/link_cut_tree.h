#pragma once

#include <functional>
#include <vector>

namespace altruct {
namespace container {

/**
 * Link-Cut tree node.
 *
 * The original forest is decomposed into preferred paths.
 * Each preferred path is stored as a balanced Splay tree,
 * referred to as aux tree, ordered by depth in the forest.
 * 
 * Left/Right pointers are used exclusively by aux trees.
 * Parent pointer is used by aux trees, unles it is aux root,
 * in which case it is used as a parent pointer in the forest.
 * Since a parent won't have aux root set as its child, we can
 * use this information to determine whether node is aux root,
 * rather than looking at parent existence which is ambiguous.
 */
template <typename T>
struct link_cut_node {
    link_cut_node* parent;
    link_cut_node* left;
    link_cut_node* right;
    int size;
    int index;
    T val;

    link_cut_node(const T& val) : val(val) {}

    bool is_root() const {
        return (this != parent->left && this != parent->right);
    }
};

/**
 * Link-Cut tree.
 *
 * Note, Splay tree is used for balancing.
 *
 * Space complexity: `O(n)`.
 * Time complexities:
 *   find:   `O(h)`
 *   insert: `O(h)`
 *   erase:  `O(h)`
 * Where `h` is height that is proportional to `log(n)` with
 * very high probability.
 *
 * param T    - element type
 * param f_up - associative functor for upward propagation; commutativity is not required.
 * param id   - neutral element with respect to `f`; i.e. `f(e, id) = f(id, e) = e`.
 *              e.g. `0` for addition, `1` for multiplication, `+inf` for minimum, etc.
 * param ALLOC - allocator type
 */
template <
    typename T,
    typename F_UP = std::function<void(T& parent, const T& left, const T& right)>,
    typename ALLOC = std::allocator<link_cut_node<T>>>
class link_cut_tree {
    typedef link_cut_node<T>* node_ptr;

    F_UP f_up;
    ALLOC alloc;
    node_ptr nil;
    std::vector<node_ptr> nodes;

public: // constructor / destructor
    ~link_cut_tree() {
        _destroy();
    }

    link_cut_tree(const F_UP& f_up, T id = T(), const ALLOC& alloc = ALLOC())
        : f_up(f_up), alloc(alloc) {
        _init(id);
    }

    link_cut_tree(size_t sz, const F_UP& f_up, T id = T(), const ALLOC& alloc = ALLOC())
        : link_cut_tree(f_up, id, alloc) {
        for (int i = 0; i < sz; i++) make_tree(id);
    }

    template<typename It>
    link_cut_tree(It begin, It end, const F_UP& f_up, T id = T(), const ALLOC& alloc = ALLOC())
        : link_cut_tree(f_up, id, alloc) {
        for (auto it = begin; it != end; ++it) make_tree(*it);
    }

public: // index based methods
    // Returns the number of nodes.
    // Node indices are 1-based, 0 corresponds to nil.
    int size() { return nodes.size() - 1; }

    // Adds a new node as a separate tree (in the original forest).
    int add(const T& val) { return make_tree(val)->index; }

    // Links the given node to the parent (in the original forest).
    void link(int node, int parent) { return link(nodes[node], nodes[parent]); }

    // Cuts the given node from its parent (in the original forest).
    void cut(int node) { return cut(nodes[node]); }

    // Finds the root of the given node (in the original forest).
    int find_root(int node) { return find_root(nodes[node])->index; }

    // Finds the parent of the given node (in the original forest).
    int find_parent(int node) { return find_parent(nodes[node])->index; }

    // Finds the lowest common ancestor of the given nodes (in the original forest).
    int find_lca(int node1, int node2) { return find_lca(nodes[node1], nodes[node2])->index; }

    // Returns the depth of the node (in the original forest).
    int depth(int node) { return depth(nodes[node]); }

    // Performs path aggregation from the given node to its root (in the original forest).
    // (the given node becomes root in its aux tree)
    T& get(int node) { return get(nodes[node]); }

    // Traverses aux tree representing the path from the node to the root.
    // int f(const T& val, const T& leftVal, const T& rightVal, int depth)
    //   returns direction to go next (-1: left, +1: right, 0: stop)
    // Returns the last node traversed.
    template <typename F>
    int traverse(int node, F f) { return traverse(nodes[node], f)->index; }

private: // node_ptr based methods
    node_ptr make_tree(const T& val) {
        return _buy(val);
    }

    // Node and parent must not already be connected (in the original forest).
    // Node must be a root (in the original forest).
    void link(node_ptr node, node_ptr parent) {
        // if (find_root(node) == find_root(parent)) throw("already connected");
        access(node);
        // if (node->left != nil) throw("not a root");
        access(parent);
        node->parent = parent; // doesn't affect aux trees
    }

    // If the node is already a root (in the original forest), nothing happens.
    void cut(node_ptr node) {
        access(node);
        // if (node->left == nil) throw("already root");
        node->left->parent = nil;
        node->left = nil;
        update(node);
    }

    node_ptr find_root(node_ptr node) {
        access(node);
        while (node->left != nil) {
            node = node->left;
        }
        access(node); // splay
        return node;
    }

    node_ptr find_parent(node_ptr node) {
        access(node);
        node = node->left;
        while (node->right != nil) {
            node = node->right;
        }
        access(node); // splay
        return node;
    }

    node_ptr find_lca(node_ptr node1, node_ptr node2) {
        if (find_root(node1) != find_root(node2)) return nil;
        access(node1);
        return access(node2);
    }

    int depth(node_ptr node) {
        access(node);
        return node->left->size;
    }

    T& get(node_ptr node) {
        access(node);
        return node->val;
    }

    template <typename F>
    node_ptr traverse(node_ptr node, F f) {
        access(node);
        int depth = node->left->size;
        while (node != nil) {
            int dir = f(node->val, node->left->val, node->right->val, depth);
            if (!dir) break;
            node = ((dir < 0) ? node->left : node->right);
            depth += ((dir < 0) ? -node->size : 1) + node->left->size;
        }
        return node;
    }

private:
    // Makes the path from the root to the node preferred.
    // Also makes the node be the last node on the preferred path.
    // (node->right will be nil)
    // Also makes the node be the root of its aux tree.
    // Returns the lowest ancestor on the previous preferred path.
    node_ptr access(node_ptr node) {
        node_ptr last = nil;
        for (node_ptr temp = node; temp != nil; temp = temp->parent) {
            splay(temp);
            temp->right = last;
            update(temp);
            last = temp;
        }
        splay(node);
        return last;
    }

    // Splays the node within its aux tree.
    void splay(node_ptr node) {
        while (node != nil && !node->is_root()) {
            node_ptr par = node->parent;
            node_ptr gp = par->parent;
            bool zig_zig = ((node == par->left) == (par == gp->left));
            if (!par->is_root()) rotate(zig_zig ? par : node);
            rotate(node);
        }
    }

    // Rotates the node around its parent.
    void rotate(node_ptr node) {
        node_ptr par = node->parent;
        node_ptr gp = par->parent;
        bool is_left = (node == par->left);
        node_ptr ch = (is_left ? node->right : node->left);
        connect(ch, par, is_left, !is_left);
        connect(par, node, !is_left, is_left);
        connect(node, gp, par == gp->left, par == gp->right);
        update(par);
        update(node); // TODO: should we go all the way up?
    }

    // Connects the node to its new parent.
    void connect(node_ptr node, node_ptr parent, bool connect_left, bool connect_right) {
        if (node != nil) node->parent = parent;
        if (connect_left) parent->left = node;
        if (connect_right) parent->right = node;
    }

    // Updates the aux tree node.
    void update(node_ptr node) {
        if (node != nil) {
            node->size = node->left->size + node->right->size + 1;
            f_up(node->val, node->left->val, node->right->val);
        }
    }

private:
    void _init(const T& id) {
        nil = alloc.allocate(1);
        alloc.construct(nil, id);
        nil->parent = nil;
        nil->left = nil;
        nil->right = nil;
        nil->size = 0;
        nil->index = 0;
        nodes.push_back(nil);
    }

    void _destroy() {
        while (!nodes.empty()) {
            _free(nodes.back());
            nodes.pop_back();
        }
        nil = nullptr;
    }

    node_ptr _buy(const T& val) {
        node_ptr t = alloc.allocate(1);
        alloc.construct(t, val);
        t->parent = nil;
        t->left = nil;
        t->right = nil;
        t->size = 1;
        t->index = (int)nodes.size();
        nodes.push_back(t);
        return t;
    }

    void _free(node_ptr t) {
        alloc.destroy(t);
        alloc.deallocate(t, 1);
    }
};

} // container
} // altruct
