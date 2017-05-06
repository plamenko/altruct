#include "structure/container/treap.h"

#include "algorithm/collections/collections.h"
#include "algorithm/random/xorshift.h"
#include "io/iostream_overloads.h"
#include "structure_test_util.h"

#include <functional>
#include <list>
#include <set>
#include <map>
#include <chrono>

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::container;
using namespace altruct::collections;
using namespace altruct::test_util;

namespace {
    template<typename K, typename T = K, int DUP = bst_duplicate_handling::IGNORE, typename CMP = std::less<K>, typename RAND = std::function<int()>, typename ALLOC = allocator<bst_node<T>>>
    class treap_dbg : public treap<K, T, DUP, CMP, RAND, ALLOC> {
    public:
        treap_dbg(const CMP& cmp = CMP(), const RAND& rnd = rand, const ALLOC& alloc = ALLOC()) :
            treap(cmp, rnd, alloc) {
        }

        template<typename It>
        treap_dbg(It begin, It end, const CMP& cmp = CMP(), const RAND& rnd = rand, const ALLOC& alloc = ALLOC()) :
            treap(begin, end, cmp, rnd, alloc) {
        }

        treap_dbg(std::initializer_list<T> list) :
            treap(list) {}

        treap_dbg(treap_dbg&& rhs) :
            treap(std::move(rhs)) {
        }

        treap_dbg(const treap_dbg& rhs) :
            treap(rhs) {
        }

        treap_dbg& operator=(treap_dbg&& rhs) {
            treap::operator=(std::move(rhs));
            return *this;
        }

        treap_dbg& operator=(const treap_dbg& rhs) {
            treap::operator=(rhs);
            return *this;
        }

        void debug_check(const_node_ptr ptr = nullptr) const {
            if (ptr == nullptr) {
                ptr = root();
                ASSERT_TRUE(nil->parent == nil) << "ERROR: nil not connected back to itself";
                ASSERT_TRUE(nil->left == nil->right) << "ERROR: nil left & right roots out of sync";
            }
            if (ptr == nil) {
                return;
            }
            if (!ptr->left->is_nil()) {
                ASSERT_FALSE(cmp(_key(ptr->val), _key(ptr->left->val))) << "ERROR: parent < left";
                ASSERT_FALSE(ptr->left->parent != ptr) << "ERROR: left not connected back to parent";
                debug_check(ptr->left);
            }
            if (!ptr->right->is_nil()) {
                ASSERT_FALSE(cmp(_key(ptr->right->val), _key(ptr->val))) << "ERROR: right < parent";
                ASSERT_FALSE(ptr->right->parent != ptr) << "ERROR: right not connected back to parent";
                debug_check(ptr->right);
            }
        }

        static void make_link(node_ptr par, node_ptr ch, bool go_left) {
            treap::make_link(par, ch, go_left);
        }
    };

    // Note: relies on iterator functioning properly
    template<typename K, typename T, int DUP, typename CMP, typename ALLOC, typename COLLECTION>
    void verify_structure(const treap_dbg<K, T, DUP, CMP, ALLOC>& t, const COLLECTION& c) {
        t.debug_check();
        list<T> va;
        for (auto it = t.cbegin(); it != t.cend(); ++it) {
            for (int i = 0; i < it.count(); i++) {
                va.push_back(*it);
            }
        }
        EXPECT_EQ(list<T>(c.cbegin(), c.cend()), va);
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

TEST(treap_test, bst_entry) {
    typedef bst_entry<int, string> entry;
    entry e0;
    EXPECT_EQ(0, e0.key());
    EXPECT_EQ("", e0.val());

    entry e1{ 42, "aaa" };
    EXPECT_EQ(42, e1.key());
    EXPECT_EQ("aaa", e1.val());

    entry e2(42, "aaa");
    EXPECT_EQ(42, e2.key());
    EXPECT_EQ("aaa", e2.val());

    entry e3(pair<int, string>(42, "aaa"));
    EXPECT_EQ(42, e3.key());
    EXPECT_EQ("aaa", e3.val());

    entry e4(pair<const int, string>(42, "aaa"));
    EXPECT_EQ(42, e4.key());
    EXPECT_EQ("aaa", e4.val());

    const entry e5{ 100, "b" };
    EXPECT_EQ(100, e5.key());
    EXPECT_EQ("b", e5.val());

    e2.val() = "x";
    EXPECT_EQ("x", e2.val());

    EXPECT_EQ("abc", (bst_key<string, string>::of("abc")));
    EXPECT_EQ(42, (bst_key<int, entry>::of({ 42, "def" })));

    ASSERT_COMPARISON_OPERATORS(-1, e0, e3);
    ASSERT_COMPARISON_OPERATORS(0, e1, e3);
    ASSERT_COMPARISON_OPERATORS(+1, e5, e3);
}

TEST(treap_test, bst_node) {
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

TEST(treap_test, bst_inorder) {
    typedef treap_dbg<int> bst;
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

TEST(treap_test, bst_iterator) {
    typedef bst_entry<int, string> entry;
    entry e{ 42, "abc" };
    bst_node<entry> t(e), r({}), s({});
    t.size = 25;
    t.left = &r;
    r.size = 6;
    t.right = &s;
    s.size = 8;

    bst_iterator<entry> it(&t);
    EXPECT_EQ(e, *it);
    EXPECT_EQ(e.key(), it->key());
    EXPECT_EQ(e.val(), it->val());
    EXPECT_TRUE(it == &t);
    EXPECT_FALSE(it == &s);
    EXPECT_EQ(11, it.count());

    bst_const_iterator<entry> cit(&t);
    EXPECT_EQ(e, *cit);
    EXPECT_EQ(e.key(), cit->key());
    EXPECT_EQ(e.val(), cit->val());
    EXPECT_TRUE(cit == &t);
    EXPECT_FALSE(cit == &s);
    EXPECT_TRUE(cit == it);
    EXPECT_EQ(11, cit.count());
}

TEST(treap_test, constructor) {
    // default
    set<int> s0;
    treap_dbg<int> t0;
    verify_structure(t0, s0);

    // range
    set<int> s1; for (int i = 0; i < 100; i++) s1.insert(rand() % 10);
    treap_dbg<int> t1(s1.begin(), s1.end());
    verify_structure(t1, s1);

    // range + comparator
    set<int, std::greater<int>> s2; for (int i = 0; i < 110; i++) s2.insert(rand() % 1000000000);
    treap_dbg<int, int, bst_duplicate_handling::IGNORE, std::greater<int>> t2(s2.begin(), s2.end(), std::greater<int>());
    verify_structure(t2, s2);

    // initializer list
    treap_dbg<int> ti{ 42, 3, 15 };
    verify_structure(ti, (set<int>{42, 3, 15}));

    // move constructor
    treap_dbg<int> t3(std::move(treap_dbg<int>(s1.begin(), s1.end())));
    verify_structure(t3, s1);

    // copy constructor
    treap_dbg<int> t4(t3);
    verify_structure(t3, s1);
    verify_structure(t4, s1);

    // move assignment
    t4 = std::move(treap_dbg<int>(s1.begin(), s1.end()));
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

TEST(treap_test, swap) {
    set<int> s1; for (int i = 0; i < 100; i++) s1.insert(rand() % 1000000000);
    set<int> s2; for (int i = 0; i < 110; i++) s2.insert(rand() % 1000000000);
    treap_dbg<int> t1(s1.begin(), s1.end());
    treap_dbg<int> t2(s2.begin(), s2.end());
    verify_structure(t1, s1);
    verify_structure(t2, s2);
    t1.swap(t2);
    verify_structure(t2, s1);
    verify_structure(t1, s2);
    std::swap(t2, t1);
    verify_structure(t1, s1);
    verify_structure(t2, s2);
}

TEST(treap_test, duplicate_handling) {
    set<int> s1; for (int i = 0; i < 110; i++) s1.insert(rand() % 1000000000);
    treap_dbg<int, int, bst_duplicate_handling::IGNORE> t1(s1.begin(), s1.end());
    verify_structure(t1, s1);

    multiset<int> s2; for (int i = 0; i < 110; i++) s2.insert(rand() % 1000000000);
    treap_dbg<int, int, bst_duplicate_handling::COUNT> t2(s2.begin(), s2.end());
    verify_structure(t2, s2);

    typedef bst_entry<int, string> entry;
    multimap<int, string> s3; for (int i = 0; i < 110; i++) s3.insert({ rand() % 10, to_string(i) });
    treap_dbg<int, entry, bst_duplicate_handling::STORE> t3(s3.begin(), s3.end());
    verify_structure(t3, s3);
}

TEST(treap_test, iterators) {
    set<int> s1; for (int i = 0; i < 110; i++) s1.insert(rand() % 1000000000);
    treap_dbg<int> t1(s1.begin(), s1.end());
    EXPECT_EQ((vector<int>(s1.begin(), s1.end())), (vector<int>(t1.begin(), t1.end())));
    EXPECT_EQ((vector<int>(s1.cbegin(), s1.cend())), (vector<int>(t1.cbegin(), t1.cend())));
    EXPECT_EQ((vector<int>(s1.rbegin(), s1.rend())), (vector<int>(t1.rbegin(), t1.rend())));
    EXPECT_EQ((vector<int>(s1.crbegin(), s1.crend())), (vector<int>(t1.crbegin(), t1.crend())));
}

TEST(treap_test, relational_operators) {
    treap_dbg<int> t{ 3, 8, 15, 16 };
    ASSERT_COMPARISON_OPERATORS(-1, (treap_dbg<int>{ }), t);                  // empty
    ASSERT_COMPARISON_OPERATORS(0, (treap_dbg<int>{ 3, 8, 15, 16 }), t);      // equal
    ASSERT_COMPARISON_OPERATORS(-1, (treap_dbg<int>{ 3, 8, 15 }), t);         // shorter
    ASSERT_COMPARISON_OPERATORS(+1, (treap_dbg<int>{ 3, 8, 15, 16, 17 }), t); // longer
    ASSERT_COMPARISON_OPERATORS(+1, (treap_dbg<int>{ 3, 9, 15 }), t);         // shorter but larger
    ASSERT_COMPARISON_OPERATORS(-1, (treap_dbg<int>{ 3, 7, 15, 16, 17 }), t); // longer but smaller

    typedef bst_entry<int, string> entry;
    typedef treap_dbg<int, entry, bst_duplicate_handling::STORE> tree;
    tree t2{ { 3, "abc" }, { 3, "d" }, { 15, "ef" }, { 16, "ghi" } };
    ASSERT_COMPARISON_OPERATORS(0, (tree{ { 3, "abc" }, { 3, "d" }, { 15, "ef" }, { 16, "ghi" } }), t2);    // equal
    ASSERT_COMPARISON_OPERATORS(+1, (tree{ { 3, "abc" }, { 3, "dx" }, { 15, "ef" }, { 16, "ghi" } }), t2);  // keys equal, but value larger
}

TEST(treap_test, query) {
    typedef bst_entry<string, int> entry;
    vector<int> c;
    list<entry> d;
    vector<string> vn{ "b", "d", "n", "q" }; // smaller keys
    vector<string> vk{ "c", "e", "o", "r" };
    list<entry> ve;
    list<entry> vi{ { "c", 1 }, { "e", 1 }, { "o", 1 }, { "r", 1 } };
    list<entry> vu{ { "c", 1 }, { "e", 3 }, { "o", 1 }, { "r", 2 } };
    list<entry> vc{ { "c", 1 }, { "e", 1 }, { "e", 1 }, { "e", 1 }, { "o", 1 }, { "r", 1 }, { "r", 1 } };
    list<entry> vs{ { "c", 1 }, { "e", 1 }, { "e", 2 }, { "e", 3 }, { "o", 1 }, { "r", 1 }, { "r", 2 } };
    // construct all with vs!
    treap_dbg<string, entry, bst_duplicate_handling::IGNORE> ti(vs.begin(), vs.end());
    treap_dbg<string, entry, bst_duplicate_handling::COUNT> tc(vs.begin(), vs.end());
    treap_dbg<string, entry, bst_duplicate_handling::STORE> ts(vs.begin(), vs.end());
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

TEST(treap_test, insert) {
    typedef bst_entry<string, int> entry;
    list<entry> d;
    list<entry> vi{ { "c", 1 }, { "e", 1 }, { "o", 1 }, { "r", 1 } };
    list<entry> vc{ { "c", 1 }, { "e", 1 }, { "e", 1 }, { "e", 1 }, { "o", 1 }, { "r", 1 }, { "r", 1 } };
    list<entry> vs{ { "c", 1 }, { "e", 1 }, { "e", 2 }, { "e", 3 }, { "o", 1 }, { "r", 1 }, { "r", 2 } };
    list<entry> vp = vs;
    do {
        // test all 420 key permutations, but keep entries with the same key in the same order
        list<entry> vd; map<string, int> m; for (auto& e : vp) e.val() = ++m[e.key()], vd.push_back({ e.key(), 1 });
        // construct empty
        treap_dbg<string, entry, bst_duplicate_handling::IGNORE> ti;
        treap_dbg<string, entry, bst_duplicate_handling::COUNT> tc;
        treap_dbg<string, entry, bst_duplicate_handling::STORE> ts;
        // feed all with vp!
        d.clear(); for (const auto& e : vp) d.push_back(*ti.insert(e)); EXPECT_EQ(vd, d); verify_structure(ti, vi);
        d.clear(); for (const auto& e : vp) d.push_back(*tc.insert(e)); EXPECT_EQ(vd, d); verify_structure(tc, vc);
        d.clear(); for (const auto& e : vp) d.push_back(*ts.insert(e)); EXPECT_EQ(vp, d); verify_structure(ts, vs);
    } while (next_permutation(vp.begin(), vp.end(), [](const entry& e1, const entry& e2){ return e1.key() < e2.key(); }));
}

TEST(treap_test, erase) {
    typedef bst_entry<string, int> entry;
    //list<entry> vs{ { "e", 1 }, { "r", 1 }, { "r", 2 }, { "e", 2 }, { "c", 1 }, { "e", 3 }, { "o", 1 }, { "r", 3 } };
    list<entry> vs{ { "c", 1 }, { "e", 1 }, { "e", 2 }, { "e", 3 }, { "o", 1 }, { "r", 1 }, { "r", 2 }, { "r", 3 } };
    do {
        // test all 1120 key permutations, but keep entries with the same key in the same order
        map<string, int> m; for (auto& e : vs) e.val() = ++m[e.key()];
        // construct all with vs!
        treap_dbg<string, entry, bst_duplicate_handling::IGNORE> ti(vs.begin(), vs.end());
        treap_dbg<string, entry, bst_duplicate_handling::COUNT> tc(vs.begin(), vs.end());
        treap_dbg<string, entry, bst_duplicate_handling::STORE> ts(vs.begin(), vs.end());
        // erase by key
        ti.erase("e"); verify_structure(ti, list<entry>{{ "c", 1 }, { "o", 1 }, { "r", 1 } });
        tc.erase("e", 1); verify_structure(tc, list<entry>{{ "c", 1 }, { "e", 1 }, { "e", 1 }, { "o", 1 }, { "r", 1 }, { "r", 1 }, { "r", 1 } });
        tc.erase("e"); verify_structure(tc, list<entry>{{ "c", 1 }, { "o", 1 }, { "r", 1 }, { "r", 1 }, { "r", 1 } });
        ts.erase("e"); verify_structure(ts, list<entry>{{ "c", 1 }, { "o", 1 }, { "r", 1 }, { "r", 2 }, { "r", 3 } });
        // erase by position
        ti.erase(ti.find_kth(2)); verify_structure(ti, list<entry>{{ "c", 1 }, { "o", 1 } });
        tc.erase(tc.find_kth(3), 1); verify_structure(tc, list<entry>{{ "c", 1 }, { "o", 1 }, { "r", 1 }, { "r", 1 } });
        tc.erase(tc.find_kth(3)); verify_structure(tc, list<entry>{{ "c", 1 }, { "o", 1 } });
        ts.erase(ts.find_kth(3)); verify_structure(ts, list<entry>{{ "c", 1 }, { "o", 1 }, { "r", 1 }, { "r", 3 } });
    } while (next_permutation(vs.begin(), vs.end(), [](const entry& e1, const entry& e2){ return e1.key() < e2.key(); }));
}

TEST(treap_test, insert_erase_with_count) {
    treap_dbg<string, string, bst_duplicate_handling::COUNT> tc;
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

namespace {
namespace x{
    struct rdtsc_clock {
        typedef unsigned long long rep;
        typedef std::ratio<1, 2666666666> period; // My machine is 2.67 GHz
        typedef std::chrono::duration<rep, period> duration;
        typedef std::chrono::time_point<rdtsc_clock> time_point;
        static const bool is_steady = true;
        static time_point now() { return time_point(duration(__rdtsc())); }
    };
}
namespace x {
    typedef rdtsc_clock clock;
    //typedef std::chrono::high_resolution_clock clock;
    typedef std::chrono::duration<double, typename clock::period> duration;
    duration since(clock::time_point t0) { return duration(clock::now() - t0); }
} // x
}

template<typename S, typename T, typename K>
void test_perf(const std::function<int()>& rnd, const string& title) {
    S ms;
    T mt(std::less<K>(), rnd);
    using namespace std::chrono;
    x::duration ds_i, ds_e, ds_c, ds_t; uint32_t cs_c = 0, cs_t = 0;
    x::duration dt_i, dt_e, dt_c, dt_t; uint32_t ct_c = 0, ct_t = 0;
    uint32_t iter = 0; double dur = 0, max_dur = 5.0; auto T0 = clock();
    while (true) {
        if (++iter % 10000 == 0 && (dur = double(clock() - T0) / CLOCKS_PER_SEC) > max_dur) break;
        // 20% erase, 40% insert, 40% find
        int prob = rnd() % 100;
        if ((prob -= 20) < 0 && mt.size() != 0) {
            int val = *mt.find_kth(rnd() % mt.size());
            auto ts = x::clock::now(); ms.erase(ms.find(val)); ds_e += x::since(ts);
            auto tt = x::clock::now(); mt.erase(mt.find(val), 1); dt_e += x::since(tt);
        } else if ((prob -= 40) < 0) {
            int val = rnd() % 1000;
            auto ts = x::clock::now(); ms.insert(val); ds_i += x::since(ts);
            auto tt = x::clock::now(); mt.insert(val); dt_i += x::since(tt);
        } else if ((prob -= 40) < 0) {
            int val = rnd() % 1000;
            auto ts = x::clock::now(); cs_c += (int)ms.count(val); ds_c += x::since(ts);
            auto tt = x::clock::now(); ct_c += (int)mt.count(val); dt_c += x::since(tt);
        }
    }
    if (true) {
        auto ts = x::clock::now(); for (const auto& e : ms) cs_t++; ds_t += x::since(ts);
        auto tt = x::clock::now(); for (const auto& e : mt) ct_t++; dt_t += x::since(tt);
    }
    verify_structure(mt, ms);
    auto dd = [](const x::duration& d) { return duration_cast<duration<double>>(d).count(); };
    printf("impl cont rand insert erase count iter hits size iter dur\n");
    printf("std %s %lf %lf %lf %lf %d %d %d %lf\n", title.c_str(), dd(ds_i), dd(ds_e), dd(ds_c), dd(ds_t), cs_c, cs_t, iter, dur);
    printf("treap %s %lf %lf %lf %lf %d %d %d %lf\n", title.c_str(), dd(dt_i), dd(dt_e), dd(dt_c), dd(dt_t), ct_c, ct_t, iter, dur);
}

TEST(treap_test, perf) {
    return; // skip perf tests by default
    altruct::random::xorshift_64star xrnd;
    auto crnd_func = [](){ return rand(); };
    auto xrnd_func = [&](){ return int(xrnd.next() % (1 << 30)); };
    auto init_crand_func = [&]() { srand(12345); return crnd_func; };
    auto init_xrand_func = [&]() { xrnd.seed(12345); return xrnd_func; };
    typedef treap_dbg<int, int, bst_duplicate_handling::IGNORE> treap_set_int;
    test_perf<set<int>, treap_set_int, int>(init_crand_func(), "set_int crand");
    test_perf<set<int>, treap_set_int, int>(init_xrand_func(), "set_int xrand");
    typedef treap_dbg<int, int, bst_duplicate_handling::COUNT> treap_countset_int;
    test_perf<multiset<int>, treap_countset_int, int>(init_crand_func(), "countset_int crand");
    test_perf<multiset<int>, treap_countset_int, int>(init_xrand_func(), "countset_int xrand");
    typedef treap_dbg<int, int, bst_duplicate_handling::STORE> treap_multiset_int;
    test_perf<multiset<int>, treap_multiset_int, int>(init_crand_func(), "multiset_int crand");
    test_perf<multiset<int>, treap_multiset_int, int>(init_xrand_func(), "multiset_int xrand");
}
