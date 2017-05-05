#include "structure/container/binary_heap.h"

#include "algorithm/collections/collections.h"
#include "io/iostream_overloads.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::container;
using namespace altruct::collections;

namespace {
    template<typename T>
    void verify_structure(const binary_heap<T>& bh, const vector<T>& v) {
        for (size_t i = bh.size() - 1; i > 0; i--) {
            EXPECT_FALSE(bh.cmp(bh.v[i], bh.v[(i - 1) / 2])) << "incorrect parent-child order";
        }
        EXPECT_EQ(sorted(v), sorted(bh.v));
    }
}

TEST(binary_heap_test, constructor) {
    vector<int> v1(100); for (int& a : v1) a = rand() % 10;
    vector<int> v2(110); for (int& a : v2) a = rand() % 1000000000;
    binary_heap<int>bh1(v1); verify_structure(bh1, v1);
    binary_heap<int>bh2(v2); verify_structure(bh2, v2);
    EXPECT_EQ(100, bh1.size());
    EXPECT_EQ(110, bh2.size());
}

TEST(binary_heap_test, insert) {
    vector<int> v1;
    binary_heap<int>bh1;
    for (int i = 0; i < 100; i++) {
        int a = rand() % 10;
        v1.push_back(a);
        bh1.insert(a);
    }
    verify_structure(bh1, v1);
}

TEST(binary_heap_test, pop_front) {
    vector<int> v1(100); for (int& a : v1) a = rand() % 10;
    binary_heap<int>bh1(v1);
    vector<int> s1;
    for (; bh1.size() > 0; bh1.pop_front()) {
        s1.push_back(bh1.front());
    }
    EXPECT_EQ(sorted(v1), s1);
}

TEST(binary_heap_test, sort) {
    vector<int> v1(100); for (int& a : v1) a = rand() % 10;
    binary_heap<int> bh1(v1);
    bh1.sort();
    EXPECT_EQ(sorted(v1), bh1.v);
}

TEST(binary_heap_test, test_perf) {
    return; // disable by default

    vector<int> va(10000000);
    for (int& a : va) a = rand() % 10;

    auto T0 = clock();
    binary_heap<int> bh(va);
    bh.sort();
    auto dT0 = clock() - T0;
    printf("%d ms\n", dT0);
    // mod10  1448 ms
    // any    3723 ms

    auto T1 = clock();
    sort(va.begin(), va.end());
    auto dT1 = clock() - T1;
    printf("%d ms\n", dT1);
    // mod10  162 ms
    // any    815 ms

    if (bh.v != va) {
        printf("ERROR\n");
        //cout << va << endl;
        //cout << bh.v << endl;
    }
}
