#include "altruct/structure/container/treap.h"

#include "altruct/algorithm/collections/collections.h"
#include "altruct/algorithm/random/xorshift.h"
#include "altruct/io/iostream_overloads.h"
#include "structure_test_util.h"

#include <functional>
#include <vector>
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
        typedef treap<K, T, DUP, CMP, RAND, ALLOC> treap_t;

        treap_dbg(const CMP& cmp = CMP(), const RAND& rnd = rand, const ALLOC& alloc = ALLOC()) :
            treap_t(cmp, rnd, alloc) {
        }

        template<typename It>
        treap_dbg(It begin, It end, const CMP& cmp = CMP(), const RAND& rnd = rand, const ALLOC& alloc = ALLOC()) :
            treap_t(begin, end, cmp, rnd, alloc) {
        }

        treap_dbg(std::initializer_list<T> list) :
            treap_t(list) {}

        treap_dbg(treap_dbg&& rhs) :
            treap_t(std::move(rhs)) {
        }

        treap_dbg(const treap_dbg& rhs) :
            treap_t(rhs) {
        }

        treap_dbg& operator=(treap_dbg&& rhs) {
            treap_t::operator=(std::move(rhs));
            return *this;
        }

        treap_dbg& operator=(const treap_dbg& rhs) {
            treap_t::operator=(rhs);
            return *this;
        }

        void debug_check() const {
            ASSERT_TRUE(end().parent() == end()) << "ERROR: nil not connected back to itself";
            ASSERT_TRUE(end().left() == end().right()) << "ERROR: nil left & right roots out of sync";
            debug_check(root());
            for (auto it = begin(); it != end(); ++it) {
                auto itn = it; ++itn; if (itn == end()) break;
                ASSERT_FALSE(tree.compare(*itn, *it)) << "ERROR: order violation";
            }
        }
        void debug_check(const_iterator it) const {
            if (it == end()) return;
            if (it.left() != end()) {
                ASSERT_FALSE(tree.compare(*it, *it.left())) << "ERROR: parent < left";
                ASSERT_FALSE(it.left().parent() != it) << "ERROR: left not connected back to parent";
                debug_check(it.left());
            }
            if (it.right() != end()) {
                ASSERT_FALSE(tree.compare(*it.right(), *it)) << "ERROR: right < parent";
                ASSERT_FALSE(it.right().parent() != it) << "ERROR: right not connected back to parent";
                debug_check(it.right());
            }
        }
    };

    // Note: relies on iterator functioning properly
    template<typename K, typename T, int DUP, typename CMP, typename RAND, typename ALLOC, typename COLLECTION>
    void verify_structure(const treap_dbg<K, T, DUP, CMP, RAND, ALLOC>& t, const COLLECTION& c) {
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

    typedef pair<const int, string> ck_entry;
    multimap<int, string> s3; for (int i = 0; i < 110; i++) s3.insert({ rand() % 10, to_string(i) });
    treap_dbg<int, ck_entry, bst_duplicate_handling::STORE> t3(s3.begin(), s3.end());
    verify_structure(t3, s3);
}

TEST(treap_test, treap_iterator) {
    typedef pair<const int, string> ck_entry;
    treap_dbg<int, ck_entry, bst_duplicate_handling::COUNT> tc;
    ck_entry e{ 42, "abc" };
    ck_entry e2{ 13, "de" };
    tc.insert(e, 11);
    tc.insert(e2, 14);

    auto it = tc.find(42);
    auto it2 = tc.find(13);
    EXPECT_EQ(e, *it);
    EXPECT_EQ(e.first, it->first);
    EXPECT_EQ(e.second, it->second);
    EXPECT_TRUE(it == tc.find(42));
    EXPECT_FALSE(it == it2);
    EXPECT_FALSE(it != tc.find(42));
    EXPECT_TRUE(it != it2);
    EXPECT_EQ(11, it.count());
    EXPECT_EQ(25, it.size());

    const auto& cc = tc;
    auto cit = cc.find(42);
    auto cit2 = cc.find(13);
    EXPECT_EQ(e, *cit);
    EXPECT_EQ(e.first, cit->first);
    EXPECT_EQ(e.second, cit->second);
    EXPECT_TRUE(cit == it);
    EXPECT_FALSE(cit == it2);
    EXPECT_FALSE(cit != it);
    EXPECT_TRUE(cit != it2);
    EXPECT_TRUE(cit == it);
    EXPECT_EQ(11, cit.count());
    EXPECT_EQ(25, cit.size());
}

TEST(treap_test, iterators) {
    set<int> s1; for (int i = 0; i < 110; i++) s1.insert(rand() % 1000000000);
    treap_dbg<int> t1(s1.begin(), s1.end());
    const treap_dbg<int> ct1(s1.begin(), s1.end());
    EXPECT_EQ((vector<int>(s1.begin(), s1.end())), (vector<int>(t1.begin(), t1.end())));
    EXPECT_EQ((vector<int>(s1.cbegin(), s1.cend())), (vector<int>(t1.cbegin(), t1.cend())));
    EXPECT_EQ((vector<int>(s1.rbegin(), s1.rend())), (vector<int>(t1.rbegin(), t1.rend())));
    EXPECT_EQ((vector<int>(s1.crbegin(), s1.crend())), (vector<int>(t1.crbegin(), t1.crend())));
    EXPECT_EQ((vector<int>(s1.begin(), s1.end())), (vector<int>(ct1.begin(), ct1.end())));
    EXPECT_EQ((vector<int>(s1.cbegin(), s1.cend())), (vector<int>(ct1.cbegin(), ct1.cend())));
    EXPECT_EQ((vector<int>(s1.rbegin(), s1.rend())), (vector<int>(ct1.rbegin(), ct1.rend())));
    EXPECT_EQ((vector<int>(s1.crbegin(), s1.crend())), (vector<int>(ct1.crbegin(), ct1.crend())));
}

template<typename T, typename It>
void collect_subtree(vector<T>& v, It it) {
    if (it == it.parent()) return;
    collect_subtree(v, it.left());
    v.push_back(*it);
    collect_subtree(v, it.right());
}

TEST(treap_test, root) {
    treap_dbg<string, string> tc;
    tc.insert("cc");
    tc.insert("aaa");
    tc.insert("b");
    tc.insert("dddd");
    vector<string> va; collect_subtree(va, tc.root());
    EXPECT_EQ((vector<string>{"aaa", "b", "cc", "dddd"}), va);
    EXPECT_EQ(tc.end(), tc.root().parent());
}

TEST(treap_test, relational_operators) {
    treap_dbg<int> t{ 3, 8, 15, 16 };
    ASSERT_COMPARISON_OPERATORS(-1, (treap_dbg<int>{ }), t);                  // empty
    ASSERT_COMPARISON_OPERATORS(0, (treap_dbg<int>{ 3, 8, 15, 16 }), t);      // equal
    ASSERT_COMPARISON_OPERATORS(-1, (treap_dbg<int>{ 3, 8, 15 }), t);         // shorter
    ASSERT_COMPARISON_OPERATORS(+1, (treap_dbg<int>{ 3, 8, 15, 16, 17 }), t); // longer
    ASSERT_COMPARISON_OPERATORS(+1, (treap_dbg<int>{ 3, 9, 15 }), t);         // shorter but larger
    ASSERT_COMPARISON_OPERATORS(-1, (treap_dbg<int>{ 3, 7, 15, 16, 17 }), t); // longer but smaller

    typedef pair<const int, string> ck_entry;
    typedef treap_dbg<int, ck_entry, bst_duplicate_handling::STORE> tree;
    tree t2{ { 3, "abc" }, { 3, "d" }, { 15, "ef" }, { 16, "ghi" } };
    ASSERT_COMPARISON_OPERATORS(0, (tree{ { 3, "abc" }, { 3, "d" }, { 15, "ef" }, { 16, "ghi" } }), t2);    // equal
    ASSERT_COMPARISON_OPERATORS(+1, (tree{ { 3, "abc" }, { 4, "d" }, { 15, "ef" }, { 16, "ghi" } }), t2);  // key larger
}

TEST(treap_test, query) {
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
    treap_dbg<string, ck_entry, bst_duplicate_handling::IGNORE> ti(vs.begin(), vs.end());
    treap_dbg<string, ck_entry, bst_duplicate_handling::COUNT> tc(vs.begin(), vs.end());
    treap_dbg<string, ck_entry, bst_duplicate_handling::STORE> ts(vs.begin(), vs.end());
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
        treap_dbg<string, ck_entry, bst_duplicate_handling::IGNORE> ti;
        treap_dbg<string, ck_entry, bst_duplicate_handling::COUNT> tc;
        treap_dbg<string, ck_entry, bst_duplicate_handling::STORE> ts;
        // feed all with vp!
        d.clear(); for (const auto& e : vp) d.push_back(*ti.insert(e)); EXPECT_EQ(vd, d); verify_structure(ti, vi);
        d.clear(); for (const auto& e : vp) d.push_back(*tc.insert(e)); EXPECT_EQ(vd, d); verify_structure(tc, vc);
        d.clear(); for (const auto& e : vp) d.push_back(*ts.insert(e)); EXPECT_EQ(vp, d); verify_structure(ts, vs);
    } while (next_permutation(vp.begin(), vp.end(), [](const entry& e1, const entry& e2){ return e1.first < e2.first; }));
}

TEST(treap_test, erase) {
    typedef pair<string, int> entry;
    typedef pair<const string, int> ck_entry;
    vector<entry> vs{ { "c", 1 }, { "e", 1 }, { "e", 2 }, { "e", 3 }, { "o", 1 }, { "r", 1 }, { "r", 2 }, { "r", 3 } };
    do {
        // test all 1120 key permutations, but keep entries with the same key in the same order
        map<string, int> m; for (auto& e : vs) e.second = ++m[e.first];
        // construct all with vs!
        treap_dbg<string, ck_entry, bst_duplicate_handling::IGNORE> ti(vs.begin(), vs.end());
        treap_dbg<string, ck_entry, bst_duplicate_handling::COUNT> tc(vs.begin(), vs.end());
        treap_dbg<string, ck_entry, bst_duplicate_handling::STORE> ts(vs.begin(), vs.end());
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

TEST(treap_test, erase_range) {
    treap_dbg<string, string, bst_duplicate_handling::COUNT> tc;
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

TEST(treap_test, insert_before) {
    treap_dbg<string, string, bst_duplicate_handling::COUNT> tc;
    tc.insert("dddd", 2);
    tc.insert("b", 3);
    tc.insert("aaa", 5);
    tc.insert("cc", 4);
    tc.insert_before(tc.find("b"), "abc", 2);
    verify_structure(tc, vector<string>{ "aaa", "aaa", "aaa", "aaa", "aaa", "abc", "abc", "b", "b", "b", "cc", "cc", "cc", "cc", "dddd", "dddd" });
}

TEST(treap_test, iterator_add_pos) {
    treap_dbg<string, string> tc;
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

namespace {
namespace x{
    struct rdtsc_clock {
        typedef unsigned long long rep;
        typedef std::ratio<1, 2666666666> period; // My machine is 2.67 GHz
        typedef std::chrono::duration<rep, period> duration;
        typedef std::chrono::time_point<rdtsc_clock> time_point;
        static const bool is_steady = true;
        // TODO: extract this to altruct::chrono and use rtdsc intrinsics on each compiler
        //static time_point now() { return time_point(duration(__rdtsc())); }
        static time_point now() { return time_point(duration(0)); }
    };
}
namespace x {
    //typedef rdtsc_clock clock;
    typedef std::chrono::high_resolution_clock clock;
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
