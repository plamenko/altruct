#include "altruct/algorithm/search/binary_search.h"

#include "gtest/gtest.h"

#include <vector>

using namespace std;
using namespace altruct::search;

TEST(binary_search_test, poly_zero) {
    auto p = [](double x) { return ((1 * x - 2) * x + 0) * x + 1; };
    EXPECT_NEAR(-0.618033988749895, lower_bound(-1.0, 0.0, 0.0, p, 0.0, false), 1e-9);
    EXPECT_NEAR(-0.618033988749895, upper_bound(-1.0, 0.0, 0.0, p, 0.0, false), 1e-9);
    EXPECT_NEAR(+1.000000000000000, lower_bound(0.0, 4/3.0, 0.0, p, 0.0, true), 1e-9);
    EXPECT_NEAR(+1.000000000000000, upper_bound(0.0, 4/3.0, 0.0, p, 0.0, true), 1e-9);
    EXPECT_NEAR(+1.618033988749895, lower_bound(4/3.0, 2.0, 0.0, p, 0.0, false), 1e-9);
    EXPECT_NEAR(+1.618033988749895, upper_bound(4/3.0, 2.0, 0.0, p, 0.0, false), 1e-9);
}

TEST(binary_search_test, index) {
    int a0[] = { 0 };
    int a1[] = { 5 };
    int a2[] = { 2, 2, 5, 5, 5, 5, 8 };
    auto deref_idx0 = [&](int idx){ return a0[idx]; };
    auto deref_idx1 = [&](int idx){ return a1[idx]; };
    auto deref_idx2 = [&](int idx){ return a2[idx]; };

    EXPECT_EQ(0, lower_bound(0, 0, 1, deref_idx0, 10));
    EXPECT_EQ(0, upper_bound(0, 0, 1, deref_idx0, 10));

    EXPECT_EQ(0, lower_bound(0, 1, 1, deref_idx1, 4));
    EXPECT_EQ(0, lower_bound(0, 1, 1, deref_idx1, 5));
    EXPECT_EQ(1, lower_bound(0, 1, 1, deref_idx1, 6));

    EXPECT_EQ(0, upper_bound(0, 1, 1, deref_idx1, 4));
    EXPECT_EQ(1, upper_bound(0, 1, 1, deref_idx1, 5));
    EXPECT_EQ(1, upper_bound(0, 1, 1, deref_idx1, 6));

    EXPECT_EQ(0, lower_bound(0, 7, 1, deref_idx2, 1));
    EXPECT_EQ(0, lower_bound(0, 7, 1, deref_idx2, 2));
    EXPECT_EQ(2, lower_bound(0, 7, 1, deref_idx2, 3));
    EXPECT_EQ(2, lower_bound(0, 7, 1, deref_idx2, 4));
    EXPECT_EQ(2, lower_bound(0, 7, 1, deref_idx2, 5));
    EXPECT_EQ(6, lower_bound(0, 7, 1, deref_idx2, 6));
    EXPECT_EQ(6, lower_bound(0, 7, 1, deref_idx2, 7));
    EXPECT_EQ(6, lower_bound(0, 7, 1, deref_idx2, 8));
    EXPECT_EQ(7, lower_bound(0, 7, 1, deref_idx2, 9));

    EXPECT_EQ(0, upper_bound(0, 7, 1, deref_idx2, 1));
    EXPECT_EQ(2, upper_bound(0, 7, 1, deref_idx2, 2));
    EXPECT_EQ(2, upper_bound(0, 7, 1, deref_idx2, 3));
    EXPECT_EQ(2, upper_bound(0, 7, 1, deref_idx2, 4));
    EXPECT_EQ(6, upper_bound(0, 7, 1, deref_idx2, 5));
    EXPECT_EQ(6, upper_bound(0, 7, 1, deref_idx2, 6));
    EXPECT_EQ(6, upper_bound(0, 7, 1, deref_idx2, 7));
    EXPECT_EQ(7, upper_bound(0, 7, 1, deref_idx2, 8));
    EXPECT_EQ(7, upper_bound(0, 7, 1, deref_idx2, 9));
}

TEST(binary_search_test, pointer) {
    int a0[] = { 0 };
    int a1[] = { 5 };
    int a[] = { 2, 2, 5, 5, 5, 5, 8 };
    auto deref_ptr = [](int* ptr){ return *ptr; };

    EXPECT_EQ(a0, lower_bound(a0, a0, 1, deref_ptr, 10));
    EXPECT_EQ(a0, upper_bound(a0, a0, 1, deref_ptr, 10));

    EXPECT_EQ(a1, lower_bound(a1, a1 + 1, 1, deref_ptr, 4));
    EXPECT_EQ(a1, lower_bound(a1, a1 + 1, 1, deref_ptr, 5));
    EXPECT_EQ(a1 + 1, lower_bound(a1, a1 + 1, 1, deref_ptr, 6));

    EXPECT_EQ(a1, upper_bound(a1, a1 + 1, 1, deref_ptr, 4));
    EXPECT_EQ(a1 + 1, upper_bound(a1, a1 + 1, 1, deref_ptr, 5));
    EXPECT_EQ(a1 + 1, upper_bound(a1, a1 + 1, 1, deref_ptr, 6));

    EXPECT_EQ(a + 0, lower_bound(a, a + 7, 1, deref_ptr, 1));
    EXPECT_EQ(a + 0, lower_bound(a, a + 7, 1, deref_ptr, 2));
    EXPECT_EQ(a + 2, lower_bound(a, a + 7, 1, deref_ptr, 3));
    EXPECT_EQ(a + 2, lower_bound(a, a + 7, 1, deref_ptr, 4));
    EXPECT_EQ(a + 2, lower_bound(a, a + 7, 1, deref_ptr, 5));
    EXPECT_EQ(a + 6, lower_bound(a, a + 7, 1, deref_ptr, 6));
    EXPECT_EQ(a + 6, lower_bound(a, a + 7, 1, deref_ptr, 7));
    EXPECT_EQ(a + 6, lower_bound(a, a + 7, 1, deref_ptr, 8));
    EXPECT_EQ(a + 7, lower_bound(a, a + 7, 1, deref_ptr, 9));

    EXPECT_EQ(a + 0, upper_bound(a, a + 7, 1, deref_ptr, 1));
    EXPECT_EQ(a + 2, upper_bound(a, a + 7, 1, deref_ptr, 2));
    EXPECT_EQ(a + 2, upper_bound(a, a + 7, 1, deref_ptr, 3));
    EXPECT_EQ(a + 2, upper_bound(a, a + 7, 1, deref_ptr, 4));
    EXPECT_EQ(a + 6, upper_bound(a, a + 7, 1, deref_ptr, 5));
    EXPECT_EQ(a + 6, upper_bound(a, a + 7, 1, deref_ptr, 6));
    EXPECT_EQ(a + 6, upper_bound(a, a + 7, 1, deref_ptr, 7));
    EXPECT_EQ(a + 7, upper_bound(a, a + 7, 1, deref_ptr, 8));
    EXPECT_EQ(a + 7, upper_bound(a, a + 7, 1, deref_ptr, 9));
}

TEST(binary_search_test, iterator) {
    vector<int> v0{};
    vector<int> v1{ 5 };
    vector<int> v{ 2, 2, 5, 5, 5, 5, 8 };
    auto deref_it = [](vector<int>::iterator it){ return *it; };

    EXPECT_EQ(v0.end(), lower_bound(v0.begin(), v0.end(), 1, deref_it, 10));
    EXPECT_EQ(v0.end(), upper_bound(v0.begin(), v0.end(), 1, deref_it, 10));

    EXPECT_EQ(v1.begin(), lower_bound(v1.begin(), v1.end(), 1, deref_it, 4));
    EXPECT_EQ(v1.begin(), lower_bound(v1.begin(), v1.end(), 1, deref_it, 5));
    EXPECT_EQ(v1.end(), lower_bound(v1.begin(), v1.end(), 1, deref_it, 6));

    EXPECT_EQ(v1.begin(), upper_bound(v1.begin(), v1.end(), 1, deref_it, 4));
    EXPECT_EQ(v1.end(), upper_bound(v1.begin(), v1.end(), 1, deref_it, 5));
    EXPECT_EQ(v1.end(), upper_bound(v1.begin(), v1.end(), 1, deref_it, 6));

    EXPECT_EQ(v.begin() + 0, lower_bound(v.begin(), v.end(), 1, deref_it, 1));
    EXPECT_EQ(v.begin() + 0, lower_bound(v.begin(), v.end(), 1, deref_it, 2));
    EXPECT_EQ(v.begin() + 2, lower_bound(v.begin(), v.end(), 1, deref_it, 3));
    EXPECT_EQ(v.begin() + 2, lower_bound(v.begin(), v.end(), 1, deref_it, 4));
    EXPECT_EQ(v.begin() + 2, lower_bound(v.begin(), v.end(), 1, deref_it, 5));
    EXPECT_EQ(v.begin() + 6, lower_bound(v.begin(), v.end(), 1, deref_it, 6));
    EXPECT_EQ(v.begin() + 6, lower_bound(v.begin(), v.end(), 1, deref_it, 7));
    EXPECT_EQ(v.begin() + 6, lower_bound(v.begin(), v.end(), 1, deref_it, 8));
    EXPECT_EQ(v.begin() + 7, lower_bound(v.begin(), v.end(), 1, deref_it, 9));
    EXPECT_EQ(v.end(), lower_bound(v.begin(), v.end(), 1, deref_it, 9));

    EXPECT_EQ(v.begin() + 0, upper_bound(v.begin(), v.end(), 1, deref_it, 1));
    EXPECT_EQ(v.begin() + 2, upper_bound(v.begin(), v.end(), 1, deref_it, 2));
    EXPECT_EQ(v.begin() + 2, upper_bound(v.begin(), v.end(), 1, deref_it, 3));
    EXPECT_EQ(v.begin() + 2, upper_bound(v.begin(), v.end(), 1, deref_it, 4));
    EXPECT_EQ(v.begin() + 6, upper_bound(v.begin(), v.end(), 1, deref_it, 5));
    EXPECT_EQ(v.begin() + 6, upper_bound(v.begin(), v.end(), 1, deref_it, 6));
    EXPECT_EQ(v.begin() + 6, upper_bound(v.begin(), v.end(), 1, deref_it, 7));
    EXPECT_EQ(v.begin() + 7, upper_bound(v.begin(), v.end(), 1, deref_it, 8));
    EXPECT_EQ(v.begin() + 7, upper_bound(v.begin(), v.end(), 1, deref_it, 9));
    EXPECT_EQ(v.end(), upper_bound(v.begin(), v.end(), 1, deref_it, 9));
}
