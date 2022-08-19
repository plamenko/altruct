#include "altruct/algorithm/search/binary_search.h"

#include "gtest/gtest.h"

#include <vector>

using namespace std;
using namespace altruct::search;

TEST(binary_search_test, poly_zero_double) {
    auto p = [](double x) { return ((1 * x - 2) * x + 0) * x + 1; };
    double eps = 0.0;
    EXPECT_NEAR(-0.618033988749895, lower_bound_num(-1.0, 0.0, eps, 0.0, p, false), 1e-9);
    EXPECT_NEAR(-0.618033988749895, upper_bound_num(-1.0, 0.0, eps, 0.0, p, false), 1e-9);
    EXPECT_NEAR(+1.000000000000000, lower_bound_num(0.0, 4/3.0, eps, 0.0, p, true), 1e-9);
    EXPECT_NEAR(+1.000000000000000, upper_bound_num(0.0, 4/3.0, eps, 0.0, p, true), 1e-9);
    EXPECT_NEAR(+1.618033988749895, lower_bound_num(4/3.0, 2.0, eps, 0.0, p, false), 1e-9);
    EXPECT_NEAR(+1.618033988749895, upper_bound_num(4/3.0, 2.0, eps, 0.0, p, false), 1e-9);
}

TEST(binary_search_test, poly_zero_double_eps) {
    auto p = [](double x) { return ((1 * x - 2) * x + 0) * x + 1; };
    double eps = 0.001;
    EXPECT_NEAR(-0.618, lower_bound_num(-1.0, 0.0, eps, 0.0, p, false), eps);
    EXPECT_NEAR(-0.618, upper_bound_num(-1.0, 0.0, eps, 0.0, p, false), eps);
    EXPECT_NEAR(+1.000, lower_bound_num(0.0, 4 / 3.0, eps, 0.0, p, true), eps);
    EXPECT_NEAR(+1.000, upper_bound_num(0.0, 4 / 3.0, eps, 0.0, p, true), eps);
    EXPECT_NEAR(+1.618, lower_bound_num(4 / 3.0, 2.0, eps, 0.0, p, false), eps);
    EXPECT_NEAR(+1.618, upper_bound_num(4 / 3.0, 2.0, eps, 0.0, p, false), eps);
}

TEST(binary_search_test, index) {
    int eps = 1;

    int a0[] = { 0 };
    auto deref_idx0 = [&](int idx) { return a0[idx]; };
    EXPECT_EQ(0, binary_search_pred(0, 0, [&](int idx) { return a0[idx] >= 10; }));
    EXPECT_EQ(0, lower_bound_num(0, 0, eps, 10, deref_idx0));
    EXPECT_EQ(0, upper_bound_num(0, 0, eps, 10, deref_idx0));

    int a1[] = { 5 };
    auto deref_idx1 = [&](int idx) { return a1[idx]; };
    EXPECT_EQ(0, binary_search_pred(0, 1, [&](int idx) { return a1[idx] >= 4; }));
    EXPECT_EQ(0, binary_search_pred(0, 1, [&](int idx) { return a1[idx] >= 5; }));
    EXPECT_EQ(1, binary_search_pred(0, 1, [&](int idx) { return a1[idx] >= 6; }));
    EXPECT_EQ(0, lower_bound_num(0, 1, eps, 4, deref_idx1));
    EXPECT_EQ(0, lower_bound_num(0, 1, eps, 5, deref_idx1));
    EXPECT_EQ(1, lower_bound_num(0, 1, eps, 6, deref_idx1));
    EXPECT_EQ(0, binary_search_pred(0, 1, [&](int idx) { return a1[idx] > 4; }));
    EXPECT_EQ(1, binary_search_pred(0, 1, [&](int idx) { return a1[idx] > 5; }));
    EXPECT_EQ(1, binary_search_pred(0, 1, [&](int idx) { return a1[idx] > 6; }));
    EXPECT_EQ(0, upper_bound_num(0, 1, eps, 4, deref_idx1));
    EXPECT_EQ(1, upper_bound_num(0, 1, eps, 5, deref_idx1));
    EXPECT_EQ(1, upper_bound_num(0, 1, eps, 6, deref_idx1));

    int a2[] = { 2, 2, 5, 5, 5, 5, 8 };
    auto deref_idx2 = [&](int idx) { return a2[idx]; };
    EXPECT_EQ(0, binary_search_pred(0, 7, [&](int idx) { return a2[idx] >= 1; }));
    EXPECT_EQ(0, binary_search_pred(0, 7, [&](int idx) { return a2[idx] >= 2; }));
    EXPECT_EQ(2, binary_search_pred(0, 7, [&](int idx) { return a2[idx] >= 3; }));
    EXPECT_EQ(2, binary_search_pred(0, 7, [&](int idx) { return a2[idx] >= 4; }));
    EXPECT_EQ(2, binary_search_pred(0, 7, [&](int idx) { return a2[idx] >= 5; }));
    EXPECT_EQ(6, binary_search_pred(0, 7, [&](int idx) { return a2[idx] >= 6; }));
    EXPECT_EQ(6, binary_search_pred(0, 7, [&](int idx) { return a2[idx] >= 7; }));
    EXPECT_EQ(6, binary_search_pred(0, 7, [&](int idx) { return a2[idx] >= 8; }));
    EXPECT_EQ(7, binary_search_pred(0, 7, [&](int idx) { return a2[idx] >= 9; }));
    EXPECT_EQ(0, lower_bound_num(0, 7, eps, 1, deref_idx2));
    EXPECT_EQ(0, lower_bound_num(0, 7, eps, 2, deref_idx2));
    EXPECT_EQ(2, lower_bound_num(0, 7, eps, 3, deref_idx2));
    EXPECT_EQ(2, lower_bound_num(0, 7, eps, 4, deref_idx2));
    EXPECT_EQ(2, lower_bound_num(0, 7, eps, 5, deref_idx2));
    EXPECT_EQ(6, lower_bound_num(0, 7, eps, 6, deref_idx2));
    EXPECT_EQ(6, lower_bound_num(0, 7, eps, 7, deref_idx2));
    EXPECT_EQ(6, lower_bound_num(0, 7, eps, 8, deref_idx2));
    EXPECT_EQ(7, lower_bound_num(0, 7, eps, 9, deref_idx2));
    EXPECT_EQ(0, binary_search_pred(0, 7, [&](int idx) { return a2[idx] > 1; }));
    EXPECT_EQ(2, binary_search_pred(0, 7, [&](int idx) { return a2[idx] > 2; }));
    EXPECT_EQ(2, binary_search_pred(0, 7, [&](int idx) { return a2[idx] > 3; }));
    EXPECT_EQ(2, binary_search_pred(0, 7, [&](int idx) { return a2[idx] > 4; }));
    EXPECT_EQ(6, binary_search_pred(0, 7, [&](int idx) { return a2[idx] > 5; }));
    EXPECT_EQ(6, binary_search_pred(0, 7, [&](int idx) { return a2[idx] > 6; }));
    EXPECT_EQ(6, binary_search_pred(0, 7, [&](int idx) { return a2[idx] > 7; }));
    EXPECT_EQ(7, binary_search_pred(0, 7, [&](int idx) { return a2[idx] > 8; }));
    EXPECT_EQ(7, binary_search_pred(0, 7, [&](int idx) { return a2[idx] > 9; }));
    EXPECT_EQ(0, upper_bound_num(0, 7, eps, 1, deref_idx2));
    EXPECT_EQ(2, upper_bound_num(0, 7, eps, 2, deref_idx2));
    EXPECT_EQ(2, upper_bound_num(0, 7, eps, 3, deref_idx2));
    EXPECT_EQ(2, upper_bound_num(0, 7, eps, 4, deref_idx2));
    EXPECT_EQ(6, upper_bound_num(0, 7, eps, 5, deref_idx2));
    EXPECT_EQ(6, upper_bound_num(0, 7, eps, 6, deref_idx2));
    EXPECT_EQ(6, upper_bound_num(0, 7, eps, 7, deref_idx2));
    EXPECT_EQ(7, upper_bound_num(0, 7, eps, 8, deref_idx2));
    EXPECT_EQ(7, upper_bound_num(0, 7, eps, 9, deref_idx2));

    int a3[] = { 8, 5, 5, 5, 5, 2, 2 };
    auto deref_idx3 = [&](int idx) { return a3[idx]; };
    EXPECT_EQ(7, lower_bound_num(0, 7, eps, 1, deref_idx3, true));
    EXPECT_EQ(5, lower_bound_num(0, 7, eps, 2, deref_idx3, true));
    EXPECT_EQ(5, lower_bound_num(0, 7, eps, 3, deref_idx3, true));
    EXPECT_EQ(5, lower_bound_num(0, 7, eps, 4, deref_idx3, true));
    EXPECT_EQ(1, lower_bound_num(0, 7, eps, 5, deref_idx3, true));
    EXPECT_EQ(1, lower_bound_num(0, 7, eps, 6, deref_idx3, true));
    EXPECT_EQ(1, lower_bound_num(0, 7, eps, 7, deref_idx3, true));
    EXPECT_EQ(0, lower_bound_num(0, 7, eps, 8, deref_idx3, true));
    EXPECT_EQ(0, lower_bound_num(0, 7, eps, 9, deref_idx3, true));
    EXPECT_EQ(7, upper_bound_num(0, 7, eps, 1, deref_idx3, true));
    EXPECT_EQ(7, upper_bound_num(0, 7, eps, 2, deref_idx3, true));
    EXPECT_EQ(5, upper_bound_num(0, 7, eps, 3, deref_idx3, true));
    EXPECT_EQ(5, upper_bound_num(0, 7, eps, 4, deref_idx3, true));
    EXPECT_EQ(5, upper_bound_num(0, 7, eps, 5, deref_idx3, true));
    EXPECT_EQ(1, upper_bound_num(0, 7, eps, 6, deref_idx3, true));
    EXPECT_EQ(1, upper_bound_num(0, 7, eps, 7, deref_idx3, true));
    EXPECT_EQ(1, upper_bound_num(0, 7, eps, 8, deref_idx3, true));
    EXPECT_EQ(0, upper_bound_num(0, 7, eps, 9, deref_idx3, true));
}

TEST(binary_search_test, pointer) {
    int a0[] = { 0 };
    int a1[] = { 5 };
    int a[] = { 2, 2, 5, 5, 5, 5, 8 };
    auto deref_ptr = [](int* ptr){ return *ptr; };

    EXPECT_EQ(a0, binary_search_pred(a0, a0, [](int* ptr) { return *ptr >= 10; }));
    EXPECT_EQ(a0, binary_search_pred(a0, a0, [](int* ptr) { return *ptr > 10; }));

    EXPECT_EQ(a1 + 0, binary_search_pred(a1, a1 + 1, [](int* ptr) { return *ptr >= 4; }));
    EXPECT_EQ(a1 + 0, binary_search_pred(a1, a1 + 1, [](int* ptr) { return *ptr >= 5; }));
    EXPECT_EQ(a1 + 1, binary_search_pred(a1, a1 + 1, [](int* ptr) { return *ptr >= 6; }));

    EXPECT_EQ(a1 + 0, binary_search_pred(a1, a1 + 1, [](int* ptr) { return *ptr > 4; }));
    EXPECT_EQ(a1 + 1, binary_search_pred(a1, a1 + 1, [](int* ptr) { return *ptr > 5; }));
    EXPECT_EQ(a1 + 1, binary_search_pred(a1, a1 + 1, [](int* ptr) { return *ptr > 6; }));

    EXPECT_EQ(a + 0, binary_search_pred(a, a + 7, [](int* ptr) { return *ptr >= 1; }));
    EXPECT_EQ(a + 0, binary_search_pred(a, a + 7, [](int* ptr) { return *ptr >= 2; }));
    EXPECT_EQ(a + 2, binary_search_pred(a, a + 7, [](int* ptr) { return *ptr >= 3; }));
    EXPECT_EQ(a + 2, binary_search_pred(a, a + 7, [](int* ptr) { return *ptr >= 4; }));
    EXPECT_EQ(a + 2, binary_search_pred(a, a + 7, [](int* ptr) { return *ptr >= 5; }));
    EXPECT_EQ(a + 6, binary_search_pred(a, a + 7, [](int* ptr) { return *ptr >= 6; }));
    EXPECT_EQ(a + 6, binary_search_pred(a, a + 7, [](int* ptr) { return *ptr >= 7; }));
    EXPECT_EQ(a + 6, binary_search_pred(a, a + 7, [](int* ptr) { return *ptr >= 8; }));
    EXPECT_EQ(a + 7, binary_search_pred(a, a + 7, [](int* ptr) { return *ptr >= 9; }));

    EXPECT_EQ(a + 0, binary_search_pred(a, a + 7, [](int* ptr) { return *ptr > 1; }));
    EXPECT_EQ(a + 2, binary_search_pred(a, a + 7, [](int* ptr) { return *ptr > 2; }));
    EXPECT_EQ(a + 2, binary_search_pred(a, a + 7, [](int* ptr) { return *ptr > 3; }));
    EXPECT_EQ(a + 2, binary_search_pred(a, a + 7, [](int* ptr) { return *ptr > 4; }));
    EXPECT_EQ(a + 6, binary_search_pred(a, a + 7, [](int* ptr) { return *ptr > 5; }));
    EXPECT_EQ(a + 6, binary_search_pred(a, a + 7, [](int* ptr) { return *ptr > 6; }));
    EXPECT_EQ(a + 6, binary_search_pred(a, a + 7, [](int* ptr) { return *ptr > 7; }));
    EXPECT_EQ(a + 7, binary_search_pred(a, a + 7, [](int* ptr) { return *ptr > 8; }));
    EXPECT_EQ(a + 7, binary_search_pred(a, a + 7, [](int* ptr) { return *ptr > 9; }));
}

TEST(binary_search_test, iterator) {
    vector<int> v0{};
    vector<int> v1{ 5 };
    vector<int> v{ 2, 2, 5, 5, 5, 5, 8 };

    EXPECT_EQ(v0.end(), binary_search_pred(v0.begin(), v0.end(), [](const vector<int>::iterator& it) { return *it >= 10; }));
    EXPECT_EQ(v0.end(), binary_search_pred(v0.begin(), v0.end(), [](const vector<int>::iterator& it) { return *it > 10; }));

    EXPECT_EQ(v1.begin(), binary_search_pred(v1.begin(), v1.end(), [](const vector<int>::iterator& it) { return *it >= 4; }));
    EXPECT_EQ(v1.begin(), binary_search_pred(v1.begin(), v1.end(), [](const vector<int>::iterator& it) { return *it >= 5; }));
    EXPECT_EQ(v1.end(), binary_search_pred(v1.begin(), v1.end(), [](const vector<int>::iterator& it) { return *it >= 6; }));

    EXPECT_EQ(v1.begin(), binary_search_pred(v1.begin(), v1.end(), [](const vector<int>::iterator& it) { return *it > 4; }));
    EXPECT_EQ(v1.end(), binary_search_pred(v1.begin(), v1.end(), [](const vector<int>::iterator& it) { return *it > 5; }));
    EXPECT_EQ(v1.end(), binary_search_pred(v1.begin(), v1.end(), [](const vector<int>::iterator& it) { return *it >= 6; }));

    EXPECT_EQ(v.begin() + 0, binary_search_pred(v.begin(), v.end(), [](const vector<int>::iterator& it) { return *it >= 1; }));
    EXPECT_EQ(v.begin() + 0, binary_search_pred(v.begin(), v.end(), [](const vector<int>::iterator& it) { return *it >= 2; }));
    EXPECT_EQ(v.begin() + 2, binary_search_pred(v.begin(), v.end(), [](const vector<int>::iterator& it) { return *it >= 3; }));
    EXPECT_EQ(v.begin() + 2, binary_search_pred(v.begin(), v.end(), [](const vector<int>::iterator& it) { return *it >= 4; }));
    EXPECT_EQ(v.begin() + 2, binary_search_pred(v.begin(), v.end(), [](const vector<int>::iterator& it) { return *it >= 5; }));
    EXPECT_EQ(v.begin() + 6, binary_search_pred(v.begin(), v.end(), [](const vector<int>::iterator& it) { return *it >= 6; }));
    EXPECT_EQ(v.begin() + 6, binary_search_pred(v.begin(), v.end(), [](const vector<int>::iterator& it) { return *it >= 7; }));
    EXPECT_EQ(v.begin() + 6, binary_search_pred(v.begin(), v.end(), [](const vector<int>::iterator& it) { return *it >= 8; }));
    EXPECT_EQ(v.begin() + 7, binary_search_pred(v.begin(), v.end(), [](const vector<int>::iterator& it) { return *it >= 9; }));
    EXPECT_EQ(v.end(), binary_search_pred(v.begin(), v.end(), [](const vector<int>::iterator& it) { return *it >= 9; }));

    EXPECT_EQ(v.begin() + 0, binary_search_pred(v.begin(), v.end(), [](const vector<int>::iterator& it) { return *it > 1; }));
    EXPECT_EQ(v.begin() + 2, binary_search_pred(v.begin(), v.end(), [](const vector<int>::iterator& it) { return *it > 2; }));
    EXPECT_EQ(v.begin() + 2, binary_search_pred(v.begin(), v.end(), [](const vector<int>::iterator& it) { return *it > 3; }));
    EXPECT_EQ(v.begin() + 2, binary_search_pred(v.begin(), v.end(), [](const vector<int>::iterator& it) { return *it > 4; }));
    EXPECT_EQ(v.begin() + 6, binary_search_pred(v.begin(), v.end(), [](const vector<int>::iterator& it) { return *it > 5; }));
    EXPECT_EQ(v.begin() + 6, binary_search_pred(v.begin(), v.end(), [](const vector<int>::iterator& it) { return *it > 6; }));
    EXPECT_EQ(v.begin() + 6, binary_search_pred(v.begin(), v.end(), [](const vector<int>::iterator& it) { return *it > 7; }));
    EXPECT_EQ(v.begin() + 7, binary_search_pred(v.begin(), v.end(), [](const vector<int>::iterator& it) { return *it > 8; }));
    EXPECT_EQ(v.begin() + 7, binary_search_pred(v.begin(), v.end(), [](const vector<int>::iterator& it) { return *it > 9; }));
    EXPECT_EQ(v.end(), binary_search_pred(v.begin(), v.end(), [](const vector<int>::iterator& it) { return *it > 9; }));
}
