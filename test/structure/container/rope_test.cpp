#include "structure/container/rope.h"

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
    template<typename T, typename RAND, typename ALLOC, typename COLLECTION>
    void verify_structure(const rope<T, RAND, ALLOC>& t, const COLLECTION& c) {
        EXPECT_EQ(list<T>(c.cbegin(), c.cend()), list<T>(t.cbegin(), t.cend()));
        EXPECT_EQ(c.size(), t.size());
        EXPECT_EQ(c.empty(), t.empty());
    }
}

TEST(rope_test, constructor) {
    // default
    vector<int> v0;
    rope<int> t0;
    verify_structure(t0, v0);

    // range
    vector<int> v1; for (int i = 0; i < 100; i++) v1.push_back(rand() % 10);
    rope<int> t1(v1.begin(), v1.end());
    verify_structure(t1, v1);

    // initializer list
    rope<int> ti{ 42, 3, 15 };
    verify_structure(ti, (vector<int>{42, 3, 15}));

    // move constructor
    rope<int> t3(std::move(rope<int>(v1.begin(), v1.end())));
    verify_structure(t3, v1);

    // copy constructor
    rope<int> t4(t3);
    verify_structure(t3, v1);
    verify_structure(t4, v1);

    // move assignment
    t4 = std::move(rope<int>(v1.begin(), v1.end()));
    verify_structure(t4, v1);

    // copy assignment
    t4 = t3;
    verify_structure(t4, v1);
    verify_structure(t3, v1);

    // clear
    t1.clear();
    verify_structure(t1, v0);
    // use after clear
    t1.push_back(12), t1.push_back(8), t1.push_back(4);
    verify_structure(t1, (vector<int>{12, 8, 4}));
}

TEST(rope_test, swap) {
    vector<int> v1; for (int i = 0; i < 100; i++) v1.push_back(rand() % 1000000000);
    vector<int> v2; for (int i = 0; i < 110; i++) v2.push_back(rand() % 1000000000);
    rope<int> t1(v1.begin(), v1.end());
    rope<int> t2(v2.begin(), v2.end());
    verify_structure(t1, v1);
    verify_structure(t2, v2);
    t1.swap(t2);
    verify_structure(t2, v1);
    verify_structure(t1, v2);
    std::swap(t2, t1);
    verify_structure(t1, v1);
    verify_structure(t2, v2);
}

TEST(rope_test, iterators) {
    vector<int> v1; for (int i = 0; i < 110; i++) v1.push_back(rand() % 1000000000);
    rope<int> t1(v1.begin(), v1.end());
    EXPECT_EQ((vector<int>(v1.begin(), v1.end())), (vector<int>(t1.begin(), t1.end())));
    EXPECT_EQ((vector<int>(v1.cbegin(), v1.cend())), (vector<int>(t1.cbegin(), t1.cend())));
    EXPECT_EQ((vector<int>(v1.rbegin(), v1.rend())), (vector<int>(t1.rbegin(), t1.rend())));
    EXPECT_EQ((vector<int>(v1.crbegin(), v1.crend())), (vector<int>(t1.crbegin(), t1.crend())));
    // TODO: test all the functionality mandated by random_access_iterator_tag
}

TEST(rope_test, relational_operators) {
    rope<int> t{ 3, 8, 15, 16 };
    ASSERT_COMPARISON_OPERATORS(-1, (rope<int>{ }), t);                  // empty
    ASSERT_COMPARISON_OPERATORS(0, (rope<int>{ 3, 8, 15, 16 }), t);      // equal
    ASSERT_COMPARISON_OPERATORS(-1, (rope<int>{ 3, 8, 15 }), t);         // shorter
    ASSERT_COMPARISON_OPERATORS(+1, (rope<int>{ 3, 8, 15, 16, 17 }), t); // longer
    ASSERT_COMPARISON_OPERATORS(+1, (rope<int>{ 3, 9, 15 }), t);         // shorter but larger
    ASSERT_COMPARISON_OPERATORS(-1, (rope<int>{ 3, 7, 15, 16, 17 }), t); // longer but smaller
}

TEST(rope_test, insert_and_query) {
    vector<pair<int, string>> ve{ { 0, "s" }, { 1, "t" }, { 2, "g" }, { 2, "n" }, { 2, "i" }, { 0, " " }, { 0, "t" }, { 1, "t" }, { 1, "e" }, { 6, "r" }, { 2, "s" } };
    rope<string> t;
    string se;
    for (auto& e : ve) {
        t.insert(e.first, e.second);
        se.insert(e.first, e.second);
    }
    string sa;
    for (int i = 0; i < t.size(); i++) {
        sa += t[i];
    }
    EXPECT_EQ(se, sa);
}

TEST(rope_test, update_and_query) {
    rope<char> t;
    string se, sa;
    string qe, qa;
    for (int iter = 0; iter < 10000; iter++) {
        int prob = rand() % 100;
        if ((prob -= 20) < 0 && !se.empty()) {
            int pos = rand() % se.size();
            t.erase(pos);
            se.erase(se.begin() + pos);
        } else if ((prob -= 40)) {
            int pos = rand() % (se.size() + 1);
            char val = 'a' + rand() % 26;
            t.insert(pos, val);
            se.insert(se.begin() + pos, val);
        } else if ((prob -= 20) < 0 && !se.empty()) {
            int pos = rand() % se.size();
            char val = 'a' + rand() % 26;
            se[pos] = val;
            t[pos] = val;
        } else if ((prob -= 20) < 0 && !se.empty()) {
            int pos = rand() % se.size();
            qe.push_back(se[pos]);
            qa.push_back(t[pos]);
        }
    }
    EXPECT_EQ(se, string(t.cbegin(), t.cend()));
    EXPECT_EQ(qe, qa);
}
