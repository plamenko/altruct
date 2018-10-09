#include "altruct/structure/container/binary_search_tree.h"

#include "altruct/algorithm/collections/collections.h"
#include "altruct/io/iostream_overloads.h"
#include "structure_test_util.h"

#include <functional>
#include <vector>

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::container;
using namespace altruct::collections;
using namespace altruct::test_util;

namespace {
    template<typename K, typename T = K, int DUP = bst_duplicate_handling::IGNORE, typename CMP = std::less<K>, typename ALLOC = allocator<bst_node<T>>>
    class binary_search_tree_dbg : public binary_search_tree<K, T, DUP, CMP, ALLOC> {
    public:
        typedef binary_search_tree<K, T, DUP, CMP, ALLOC> bst_t;
        using const_node_ptr = typename bst_t::const_node_ptr;
        using node_ptr = typename bst_t::node_ptr;

        binary_search_tree_dbg(const CMP& cmp = CMP(), const ALLOC& alloc = ALLOC()) :
            bst_t(cmp, alloc) {
        }

        template<typename It>
        binary_search_tree_dbg(It begin, It end, const CMP& cmp = CMP(), const ALLOC& alloc = ALLOC()) :
            bst_t(begin, end, cmp, alloc) {
        }

        binary_search_tree_dbg(std::initializer_list<T> list) :
            bst_t(list) {}

        binary_search_tree_dbg(binary_search_tree_dbg&& rhs) :
            bst_t(std::move(rhs)) {
        }

        binary_search_tree_dbg(const binary_search_tree_dbg& rhs) :
            bst_t(rhs) {
        }

        binary_search_tree_dbg& operator=(binary_search_tree_dbg&& rhs) {
            bst_t::operator=(std::move(rhs));
            return *this;
        }

        binary_search_tree_dbg& operator=(const binary_search_tree_dbg& rhs) {
            bst_t::operator=(rhs);
            return *this;
        }

        void debug_check() const {
            ASSERT_TRUE(bst_t::nil->parent == bst_t::nil) << "ERROR: nil not connected back to itself";
            ASSERT_TRUE(bst_t::nil->left == bst_t::nil->right) << "ERROR: nil left & right roots out of sync";
            debug_check(bst_t::root_ptr());
            for (auto it = this->begin(); it != this->end(); ++it) {
                auto itn = it; ++itn; if (itn == this->end()) break;
                ASSERT_FALSE(this->compare(*itn, *it)) << "ERROR: order violation";
            }
        }
        void debug_check(const_node_ptr ptr) const {
            if (ptr->is_nil()) return;
            if (!ptr->left->is_nil()) {
                ASSERT_FALSE(this->compare(ptr->val, ptr->left->val)) << "ERROR: parent < left";
                ASSERT_FALSE(ptr->left->parent != ptr) << "ERROR: left not connected back to parent";
                debug_check(ptr->left);
            }
            if (!ptr->right->is_nil()) {
                ASSERT_FALSE(this->compare(ptr->right->val, ptr->val)) << "ERROR: right < parent";
                ASSERT_FALSE(ptr->right->parent != ptr) << "ERROR: right not connected back to parent";
                debug_check(ptr->right);
            }
        }

        static void make_link(node_ptr par, node_ptr ch, bool go_left) {
            bst_t::make_link(par, ch, go_left);
        }
    };

    // Note: relies on iterator functioning properly
    template<typename K, typename T, int DUP, typename CMP, typename ALLOC, typename COLLECTION>
    void verify_structure(const binary_search_tree_dbg<K, T, DUP, CMP, ALLOC>& t, const COLLECTION& c) {
        t.debug_check();
        vector<T> va;
        for (auto it = t.cbegin(); it != t.cend(); ++it) {
            for (int i = 0; i < it.count(); i++) {
                va.push_back(*it);
            }
        }
        EXPECT_EQ(vector<T>(c.cbegin(), c.cend()), va);
        EXPECT_EQ(c.size(), t.size());
        EXPECT_EQ(c.empty(), t.empty());
    }

    template<typename T>
    bst_node<T> new_node(const T& val, bst_node<T>* nil) {
        bst_node<T> t(val);
        t.parent = nil;
        t.left = nil;
        t.right = nil;
        t.balance = 0;
        t.size = 0;
        return t;
    }

    template<typename It>
    struct iter_pair {
        std::pair<It, It> r;
        iter_pair(const std::pair<It, It>& r) : r(r) {}
        It begin() { return r.first; }
        It end() { return r.second; }
    };
    template<typename It>
    iter_pair<It> make_range(const std::pair<It, It> t) {
        return t;
    }
}

TEST(binary_search_tree_test, bst_node) {
    bst_node<int> nodes[4] {10, 20, 30, 40};
    (nodes + 0)->parent = (nodes + 0);
    EXPECT_TRUE((nodes + 0)->is_nil());
    (nodes + 1)->parent = (nodes + 0);
    EXPECT_FALSE((nodes + 1)->is_nil());

    (nodes + 1)->size = 25;
    (nodes + 1)->left = (nodes + 2);
    (nodes + 2)->size = 6;
    (nodes + 1)->right = (nodes + 3);
    (nodes + 3)->size = 8;
    EXPECT_EQ(11, (nodes + 1)->count());
}

TEST(binary_search_tree_test, bst_inorder) {
    typedef binary_search_tree_dbg<int> bst;
    bst_node<int> nil(-1); nil.parent = &nil;
    vector<bst_node<int>> nodes;
    for (int i = 0; i < 12; i++) {
        nodes.push_back(new_node(i, &nil));
    }
    bst::make_link(&nil, &nodes[5], true);
    bst::make_link(&nodes[5], &nodes[1], true);
    bst::make_link(&nodes[1], &nodes[0], true);
    bst::make_link(&nodes[1], &nodes[4], false);
    bst::make_link(&nodes[4], &nodes[3], true);
    bst::make_link(&nodes[3], &nodes[2], true);
    bst::make_link(&nodes[5], &nodes[9], false);
    bst::make_link(&nodes[9], &nodes[6], true);
    bst::make_link(&nodes[9], &nodes[11], false);
    bst::make_link(&nodes[6], &nodes[7], false);
    bst::make_link(&nodes[7], &nodes[8], false);
    bst::make_link(&nodes[11], &nodes[10], true);
    bst_node<int>* tmp = &nil;
    for (int i = -1; i < 12; i++) {
        EXPECT_EQ(i, tmp->val);
        tmp = bst_iterator_util<int>::inorder_next(tmp);
    }
    for (int i = 11; i >= -1; i--) {
        tmp = bst_iterator_util<int>::inorder_prev(tmp);
        EXPECT_EQ(i, tmp->val);
    }
    bst_iterator<int> it(tmp);
    for (int i = -1; i < 12; i++) {
        EXPECT_EQ(i, *it);
        ++it;
    }
    for (int i = 11; i >= -1; i--) {
        --it;
        EXPECT_EQ(i, *it);
    }
    for (int i = -1; i < 12; i++) {
        EXPECT_EQ(i, *it);
        it++;
    }
    for (int i = 11; i >= -1; i--) {
        it--;
        EXPECT_EQ(i, *it);
    }
    bst_const_iterator<int> cit(tmp);
    for (int i = -1; i < 12; i++) {
        EXPECT_EQ(i, *cit);
        ++cit;
    }
    for (int i = 11; i >= -1; i--) {
        --cit;
        EXPECT_EQ(i, *cit);
    }
    for (int i = -1; i < 12; i++) {
        EXPECT_EQ(i, *cit);
        cit++;
    }
    for (int i = 11; i >= -1; i--) {
        cit--;
        EXPECT_EQ(i, *cit);
    }
}

TEST(binary_search_tree_test, bst_iterator) {
    typedef pair<const int, string> ck_entry;
    ck_entry e{ 42, "abc" };
    bst_node<ck_entry> t(e), r({}), s({});
    t.size = 25;
    t.left = &r;
    r.size = 6;
    t.right = &s;
    s.size = 8;
    s.parent = &t;

    bst_iterator<ck_entry> it(&t), itr(&r), its(&s);
    EXPECT_EQ(e, *it);
    EXPECT_EQ(e.first, it->first);
    EXPECT_EQ(e.second, it->second);
    EXPECT_TRUE(it == &t);
    EXPECT_FALSE(it == &s);
    EXPECT_FALSE(it != &t);
    EXPECT_TRUE(it != &s);
    EXPECT_EQ(11, it.count());
    EXPECT_EQ(25, it.size());
    EXPECT_EQ(itr, it.left());
    EXPECT_EQ(its, it.right());
    EXPECT_EQ(it, its.parent());

    bst_const_iterator<ck_entry> cit(&t), citr(&r), cits(&s);
    EXPECT_EQ(e, *cit);
    EXPECT_EQ(e.first, cit->first);
    EXPECT_EQ(e.second, cit->second);
    EXPECT_TRUE(cit == &t);
    EXPECT_FALSE(cit == &s);
    EXPECT_FALSE(cit != &t);
    EXPECT_TRUE(cit != &s);
    EXPECT_TRUE(cit == it);
    EXPECT_EQ(11, cit.count());
    EXPECT_EQ(25, cit.size());
    EXPECT_EQ(citr, cit.left());
    EXPECT_EQ(cits, cit.right());
    EXPECT_EQ(cit, cits.parent());
}

TEST(binary_search_tree_test, constructor) {
    // default
    set<int> s0;
    binary_search_tree_dbg<int> t0;
    verify_structure(t0, s0);

    // range
    set<int> s1; for (int i = 0; i < 100; i++) s1.insert(rand() % 10);
    binary_search_tree_dbg<int> t1(s1.begin(), s1.end());
    verify_structure(t1, s1);

    // range + comparator
    set<int, std::greater<int>> s2; for (int i = 0; i < 110; i++) s2.insert(rand() % 1000000000);
    binary_search_tree_dbg<int, int, bst_duplicate_handling::IGNORE, std::greater<int>> t2(s2.begin(), s2.end(), std::greater<int>());
    verify_structure(t2, s2);

    // initializer list
    binary_search_tree_dbg<int> ti{ 42, 3, 15 };
    verify_structure(ti, (set<int>{42, 3, 15}));

    // move constructor
    binary_search_tree_dbg<int> t3(std::move(binary_search_tree_dbg<int>(s1.begin(), s1.end())));
    verify_structure(t3, s1);

    // copy constructor
    binary_search_tree_dbg<int> t4(t3);
    verify_structure(t3, s1);
    verify_structure(t4, s1);

    // move assignment
    t4 = std::move(binary_search_tree_dbg<int>(s1.begin(), s1.end()));
    verify_structure(t4, s1);

    // copy assignment
    t4 = t3;
    verify_structure(t4, s1);
    verify_structure(t3, s1);

    // clear
    t1.clear();
    verify_structure(t1, s0);
    // use after clear
    t1.insert(12), t1.insert(8), t1.insert(4);
    verify_structure(t1, (set<int>{12, 8, 4}));
}

TEST(binary_search_tree_test, swap) {
    set<int> s1; for (int i = 0; i < 100; i++) s1.insert(rand() % 1000000000);
    set<int> s2; for (int i = 0; i < 110; i++) s2.insert(rand() % 1000000000);
    binary_search_tree_dbg<int> t1(s1.begin(), s1.end());
    binary_search_tree_dbg<int> t2(s2.begin(), s2.end());
    verify_structure(t1, s1);
    verify_structure(t2, s2);
    t1.swap(t2);
    verify_structure(t2, s1);
    verify_structure(t1, s2);
    std::swap(t2, t1);
    verify_structure(t1, s1);
    verify_structure(t2, s2);
}

TEST(binary_search_tree_test, duplicate_handling) {
    set<int> s1; for (int i = 0; i < 110; i++) s1.insert(rand() % 1000000000);
    binary_search_tree_dbg<int, int, bst_duplicate_handling::IGNORE> t1(s1.begin(), s1.end());
    verify_structure(t1, s1);

    multiset<int> s2; for (int i = 0; i < 110; i++) s2.insert(rand() % 1000000000);
    binary_search_tree_dbg<int, int, bst_duplicate_handling::COUNT> t2(s2.begin(), s2.end());
    verify_structure(t2, s2);

    typedef pair<const int, string> ck_entry;
    multimap<int, string> s3; for (int i = 0; i < 110; i++) s3.insert({ rand() % 10, to_string(i) });
    binary_search_tree_dbg<int, ck_entry, bst_duplicate_handling::STORE> t3(s3.begin(), s3.end());
    verify_structure(t3, s3);
}

TEST(binary_search_tree_test, iterators) {
    set<int> s1; for (int i = 0; i < 110; i++) s1.insert(rand() % 1000000000);
    binary_search_tree_dbg<int> t1(s1.begin(), s1.end());
    const binary_search_tree_dbg<int> ct1(s1.begin(), s1.end());
    EXPECT_EQ((vector<int>(s1.begin(), s1.end())), (vector<int>(t1.begin(), t1.end())));
    EXPECT_EQ((vector<int>(s1.cbegin(), s1.cend())), (vector<int>(t1.cbegin(), t1.cend())));
    EXPECT_EQ((vector<int>(s1.rbegin(), s1.rend())), (vector<int>(t1.rbegin(), t1.rend())));
    EXPECT_EQ((vector<int>(s1.crbegin(), s1.crend())), (vector<int>(t1.crbegin(), t1.crend())));
    EXPECT_EQ((vector<int>(s1.begin(), s1.end())), (vector<int>(ct1.begin(), ct1.end())));
    EXPECT_EQ((vector<int>(s1.cbegin(), s1.cend())), (vector<int>(ct1.cbegin(), ct1.cend())));
    EXPECT_EQ((vector<int>(s1.rbegin(), s1.rend())), (vector<int>(ct1.rbegin(), ct1.rend())));
    EXPECT_EQ((vector<int>(s1.crbegin(), s1.crend())), (vector<int>(ct1.crbegin(), ct1.crend())));
}

TEST(binary_search_tree_test, root) {
    binary_search_tree_dbg<string, string> tc;
    tc.insert("cc");
    tc.insert("aaa");
    tc.insert("b");
    tc.insert("dddd");
    EXPECT_EQ(tc.find("cc"), tc.root());
    EXPECT_EQ(tc.find("aaa"), tc.root().left());
    EXPECT_EQ(tc.find("dddd"), tc.root().right());
    EXPECT_EQ(tc.end(), tc.root().parent());
}

TEST(binary_search_tree_test, relational_operators) {
    binary_search_tree_dbg<int> t{ 3, 8, 15, 16 };
    ASSERT_COMPARISON_OPERATORS(-1, (binary_search_tree_dbg<int>{ }), t);                  // empty
    ASSERT_COMPARISON_OPERATORS(0, (binary_search_tree_dbg<int>{ 3, 8, 15, 16 }), t);      // equal
    ASSERT_COMPARISON_OPERATORS(-1, (binary_search_tree_dbg<int>{ 3, 8, 15 }), t);         // shorter
    ASSERT_COMPARISON_OPERATORS(+1, (binary_search_tree_dbg<int>{ 3, 8, 15, 16, 17 }), t); // longer
    ASSERT_COMPARISON_OPERATORS(+1, (binary_search_tree_dbg<int>{ 3, 9, 15 }), t);         // shorter but larger
    ASSERT_COMPARISON_OPERATORS(-1, (binary_search_tree_dbg<int>{ 3, 7, 15, 16, 17 }), t); // longer but smaller

    typedef pair<const int, string> ck_entry;
    typedef binary_search_tree_dbg<int, ck_entry, bst_duplicate_handling::STORE> tree;
    //typedef multimap<int, string> tree;
    tree t2{ { 3, "abc" }, { 3, "d" }, { 15, "ef" }, { 16, "ghi" } };
    ASSERT_COMPARISON_OPERATORS(0, (tree{ { 3, "abc" }, { 3, "d" }, { 15, "ef" }, { 16, "ghi" } }), t2);    // equal
    ASSERT_COMPARISON_OPERATORS(+1, (tree{ { 3, "abc" }, { 4, "d" }, { 15, "ef" }, { 16, "ghi" } }), t2);   // key larger
}

TEST(binary_search_tree_test, query) {
    typedef pair<string, int> entry;
    typedef pair<const string, int> ck_entry;
    vector<int> c;
    vector<entry> d;
    vector<string> vn{ "b", "d", "n", "q" }; // smaller keys
    vector<string> vk{ "c", "e", "o", "r" };
    vector<entry> ve;
    vector<entry> vi{ { "c", 1 }, { "e", 1 }, { "o", 1 }, { "r", 1 } };
    vector<entry> vu{ { "c", 1 }, { "e", 3 }, { "o", 1 }, { "r", 2 } };
    vector<entry> vc{ { "c", 1 }, { "e", 1 }, { "e", 1 }, { "e", 1 }, { "o", 1 }, { "r", 1 }, { "r", 1 } };
    vector<entry> vs{ { "c", 1 }, { "e", 1 }, { "e", 2 }, { "e", 3 }, { "o", 1 }, { "r", 1 }, { "r", 2 } };
    // construct all with vs!
    binary_search_tree_dbg<string, ck_entry, bst_duplicate_handling::IGNORE> ti(vs.begin(), vs.end());
    binary_search_tree_dbg<string, ck_entry, bst_duplicate_handling::COUNT> tc(vs.begin(), vs.end());
    binary_search_tree_dbg<string, ck_entry, bst_duplicate_handling::STORE> ts(vs.begin(), vs.end());
    // count_less_or_equal
    c.clear(); for (auto& k : vk) c.push_back(ti.count_less_or_equal(k)); EXPECT_EQ((vector<int>{ 1, 2, 3, 4 }), c);
    c.clear(); for (auto& k : vk) c.push_back(tc.count_less_or_equal(k)); EXPECT_EQ((vector<int>{ 1, 4, 5, 7 }), c);
    c.clear(); for (auto& k : vk) c.push_back(ts.count_less_or_equal(k)); EXPECT_EQ((vector<int>{ 1, 4, 5, 7 }), c);
    // count_less
    c.clear(); for (auto& k : vk) c.push_back(ti.count_less(k)); EXPECT_EQ((vector<int>{ 0, 1, 2, 3 }), c);
    c.clear(); for (auto& k : vk) c.push_back(tc.count_less(k)); EXPECT_EQ((vector<int>{ 0, 1, 4, 5 }), c);
    c.clear(); for (auto& k : vk) c.push_back(ts.count_less(k)); EXPECT_EQ((vector<int>{ 0, 1, 4, 5 }), c);
    // count
    c.clear(); for (auto& k : vk) c.push_back(ti.count(k)); EXPECT_EQ((vector<int>{ 1, 1, 1, 1 }), c);
    c.clear(); for (auto& k : vk) c.push_back(tc.count(k)); EXPECT_EQ((vector<int>{ 1, 3, 1, 2 }), c);
    c.clear(); for (auto& k : vk) c.push_back(ts.count(k)); EXPECT_EQ((vector<int>{ 1, 3, 1, 2 }), c);
    // find_kth
    d.clear(); for (int k = 0; k < ti.size(); k++) d.push_back(*ti.find_kth(k)); EXPECT_EQ(vi, d);
    d.clear(); for (int k = 0; k < tc.size(); k++) d.push_back(*tc.find_kth(k)); EXPECT_EQ(vc, d);
    d.clear(); for (int k = 0; k < ts.size(); k++) d.push_back(*ts.find_kth(k)); EXPECT_EQ(vs, d);
    // find
    d.clear(); for (auto& k : vk) d.push_back(*ti.find(k)); EXPECT_EQ(vi, d);
    d.clear(); for (auto& k : vk) d.push_back(*tc.find(k)); EXPECT_EQ(vi, d);
    d.clear(); for (auto& k : vk) d.push_back(*ts.find(k)); EXPECT_EQ(vi, d);
    // lower_bound
    d.clear(); for (auto& k : vk) d.push_back(*ti.lower_bound(k)); EXPECT_EQ(vi, d);
    d.clear(); for (auto& k : vk) d.push_back(*tc.lower_bound(k)); EXPECT_EQ(vi, d);
    d.clear(); for (auto& k : vk) d.push_back(*ts.lower_bound(k)); EXPECT_EQ(vi, d);
    // lower_bound (for non-existing)
    d.clear(); for (auto& k : vn) d.push_back(*ti.lower_bound(k)); EXPECT_EQ(vi, d);
    d.clear(); for (auto& k : vn) d.push_back(*tc.lower_bound(k)); EXPECT_EQ(vi, d);
    d.clear(); for (auto& k : vn) d.push_back(*ts.lower_bound(k)); EXPECT_EQ(vi, d);
    // upper_bound
    d.clear(); for (auto& k : vk) d.push_back(*--ti.upper_bound(k)); EXPECT_EQ(vi, d);
    d.clear(); for (auto& k : vk) d.push_back(*--tc.upper_bound(k)); EXPECT_EQ(vi, d);
    d.clear(); for (auto& k : vk) d.push_back(*--ts.upper_bound(k)); EXPECT_EQ(vu, d);
    // upper_bound (for non-existing)
    d.clear(); for (auto& k : vn) d.push_back(*ti.upper_bound(k)); EXPECT_EQ(vi, d);
    d.clear(); for (auto& k : vn) d.push_back(*tc.upper_bound(k)); EXPECT_EQ(vi, d);
    d.clear(); for (auto& k : vn) d.push_back(*ts.upper_bound(k)); EXPECT_EQ(vi, d);
    // equal_range
    d.clear(); for (auto& k : vk) for (auto e : make_range(ti.equal_range(k))) d.push_back(e); EXPECT_EQ(vi, d);
    d.clear(); for (auto& k : vk) for (auto e : make_range(tc.equal_range(k))) d.push_back(e); EXPECT_EQ(vi, d);
    d.clear(); for (auto& k : vk) for (auto e : make_range(ts.equal_range(k))) d.push_back(e); EXPECT_EQ(vs, d);
    // equal_range (for non-existing)
    d.clear(); for (auto& k : vn) for (auto e : make_range(ti.equal_range(k))) d.push_back(e); EXPECT_EQ(ve, d);
    d.clear(); for (auto& k : vn) for (auto e : make_range(tc.equal_range(k))) d.push_back(e); EXPECT_EQ(ve, d);
    d.clear(); for (auto& k : vn) for (auto e : make_range(ts.equal_range(k))) d.push_back(e); EXPECT_EQ(ve, d);
}

TEST(binary_search_tree_test, insert) {
    typedef pair<string, int> entry;
    typedef pair<const string, int> ck_entry;
    vector<entry> d;
    vector<entry> vi{ { "c", 1 }, { "e", 1 }, { "o", 1 }, { "r", 1 } };
    vector<entry> vc{ { "c", 1 }, { "e", 1 }, { "e", 1 }, { "e", 1 }, { "o", 1 }, { "r", 1 }, { "r", 1 } };
    vector<entry> vs{ { "c", 1 }, { "e", 1 }, { "e", 2 }, { "e", 3 }, { "o", 1 }, { "r", 1 }, { "r", 2 } };
    vector<entry> vp = vs;
    do {
        // test all 420 key permutations, but keep entries with the same key in the same order
        vector<entry> vd; map<string, int> m; for (auto& e : vp) e.second = ++m[e.first], vd.push_back({ e.first, 1 });
        // construct empty
        binary_search_tree_dbg<string, ck_entry, bst_duplicate_handling::IGNORE> ti;
        binary_search_tree_dbg<string, ck_entry, bst_duplicate_handling::COUNT> tc;
        binary_search_tree_dbg<string, ck_entry, bst_duplicate_handling::STORE> ts;
        // feed all with vp!
        d.clear(); for (const auto& e : vp) d.push_back(*ti.insert(e)); EXPECT_EQ(vd, d); verify_structure(ti, vi);
        d.clear(); for (const auto& e : vp) d.push_back(*tc.insert(e)); EXPECT_EQ(vd, d); verify_structure(tc, vc);
        d.clear(); for (const auto& e : vp) d.push_back(*ts.insert(e)); EXPECT_EQ(vp, d); verify_structure(ts, vs);
    } while (next_permutation(vp.begin(), vp.end(), [](const entry& e1, const entry& e2){ return e1.first < e2.first; }));
}

TEST(binary_search_tree_test, erase) {
    typedef pair<string, int> entry;
    typedef pair<const string, int> ck_entry;
    vector<entry> vs{ { "c", 1 }, { "e", 1 }, { "e", 2 }, { "e", 3 }, { "o", 1 }, { "r", 1 }, { "r", 2 }, { "r", 3 } };
    do {
        // test all 1120 key permutations, but keep entries with the same key in the same order
        map<string, int> m; for (auto& e : vs) e.second = ++m[e.first];
        // construct all with vs!
        binary_search_tree_dbg<string, ck_entry, bst_duplicate_handling::IGNORE> ti(vs.begin(), vs.end());
        binary_search_tree_dbg<string, ck_entry, bst_duplicate_handling::COUNT> tc(vs.begin(), vs.end());
        binary_search_tree_dbg<string, ck_entry, bst_duplicate_handling::STORE> ts(vs.begin(), vs.end());
        // erase by key
        ti.erase("e"); verify_structure(ti, vector<entry>{{ "c", 1 }, { "o", 1 }, { "r", 1 } });
        tc.erase("e", 1); verify_structure(tc, vector<entry>{{ "c", 1 }, { "e", 1 }, { "e", 1 }, { "o", 1 }, { "r", 1 }, { "r", 1 }, { "r", 1 } });
        tc.erase("e"); verify_structure(tc, vector<entry>{{ "c", 1 }, { "o", 1 }, { "r", 1 }, { "r", 1 }, { "r", 1 } });
        ts.erase("e"); verify_structure(ts, vector<entry>{{ "c", 1 }, { "o", 1 }, { "r", 1 }, { "r", 2 }, { "r", 3 } });
        // erase by position
        ti.erase(ti.find_kth(2)); verify_structure(ti, vector<entry>{{ "c", 1 }, { "o", 1 } });
        tc.erase(tc.find_kth(3), 1); verify_structure(tc, vector<entry>{{ "c", 1 }, { "o", 1 }, { "r", 1 }, { "r", 1 } });
        tc.erase(tc.find_kth(3)); verify_structure(tc, vector<entry>{{ "c", 1 }, { "o", 1 } });
        ts.erase(ts.find_kth(3)); verify_structure(ts, vector<entry>{{ "c", 1 }, { "o", 1 }, { "r", 1 }, { "r", 3 } });
    } while (next_permutation(vs.begin(), vs.end(), [](const entry& e1, const entry& e2){ return e1.first < e2.first; }));
}

TEST(binary_search_tree_test, insert_erase_with_count) {
    binary_search_tree_dbg<string, string, bst_duplicate_handling::COUNT> tc;
    tc.insert("aaa", 5);
    tc.insert("b", 2);
    tc.insert("cc", 4);
    verify_structure(tc, vector<string>{ "aaa", "aaa", "aaa", "aaa", "aaa", "b", "b", "cc", "cc", "cc", "cc" });
    tc.erase("d", 5);
    verify_structure(tc, vector<string>{ "aaa", "aaa", "aaa", "aaa", "aaa", "b", "b", "cc", "cc", "cc", "cc" });
    tc.erase("aaa", 3);
    tc.erase("cc", 1);
    verify_structure(tc, vector<string>{ "aaa", "aaa", "b", "b", "cc", "cc", "cc" });
    tc.insert("b", 1);
    tc.insert("e", 2);
    verify_structure(tc, vector<string>{ "aaa", "aaa", "b", "b", "b", "cc", "cc", "cc", "e", "e" });
}

TEST(binary_search_tree_test, erase_range) {
    binary_search_tree_dbg<string, string, bst_duplicate_handling::COUNT> tc;
    tc.insert("dddd", 2);
    tc.insert("b", 3);
    tc.insert("aaa", 5);
    tc.insert("cc", 4);
    auto b = tc.find("b");
    auto e = tc.find("dddd");
    tc.erase(b, e, 2);
    verify_structure(tc, vector<string>{ "aaa", "aaa", "aaa", "aaa", "aaa", "b", "cc", "cc", "dddd", "dddd" });
    tc.erase(tc.begin(), tc.end(), 1);
    verify_structure(tc, vector<string>{ "aaa", "aaa", "aaa", "aaa", "cc", "dddd" });
}

TEST(binary_search_tree_test, insert_before) {
    binary_search_tree_dbg<string, string, bst_duplicate_handling::COUNT> tc;
    tc.insert("dddd", 2);
    tc.insert("b", 3);
    tc.insert("aaa", 5);
    tc.insert("cc", 4);
    tc.insert_before(tc.find("b"), "abc", 2);
    verify_structure(tc, vector<string>{ "aaa", "aaa", "aaa", "aaa", "aaa", "abc", "abc", "b", "b", "b", "cc", "cc", "cc", "cc", "dddd", "dddd" });
}

TEST(binary_search_tree_test, iterator_add_pos) {
    binary_search_tree_dbg<string, string> tc;
    tc.insert("dddd", 2);
    tc.insert("b", 3);
    tc.insert("aaa", 5);
    tc.insert("cc", 4);
    EXPECT_EQ(0, tc.find("aaa").pos());
    EXPECT_EQ(1, tc.find("b").pos());
    EXPECT_EQ(2, tc.find("cc").pos());
    EXPECT_EQ(3, tc.find("dddd").pos());
    EXPECT_EQ(4, tc.find("c").pos());
    auto it = tc.find("b");
    EXPECT_EQ("b", *it.add(0));
    EXPECT_EQ("aaa", *it.add(-1));
    EXPECT_EQ("cc", *it.add(+1));
    EXPECT_EQ("dddd", *it.add(+2));
    EXPECT_EQ(tc.end(), it.add(+3));
}
