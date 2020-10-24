#include "altruct/algorithm/math/squares_r.h"
#include "altruct/algorithm/math/factorization.h"

#include <algorithm>
#include <vector>
#include <map>

#include "gtest/gtest.h"

#include "altruct/io/iostream_overloads.h"

using namespace std;
using namespace altruct::math;

TEST(squares_r_test, cornacchia) {
    EXPECT_EQ(make_pair(7, 3), cornacchia(6, 103));
}

TEST(squares_r_test, squares_r) {
    vector<int> va, vu;
    for (int i = 1; i <= 30; i++) {
        auto vf = factor_integer_slow(i);
        va.push_back(squares_r(vf, false));
        vu.push_back(squares_r(vf, true));
    }
    EXPECT_EQ((vector<int>{4, 4, 0, 4, 8, 0, 0, 4, 4, 8, 0, 0, 8, 0, 0, 4, 8, 4, 0, 8, 0, 0, 0, 0, 12, 8, 0, 0, 8, 0}), va);
    EXPECT_EQ((vector<int>{1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 2, 1, 0, 0, 1, 0}), vu);
}

TEST(squares_r_test, squares_r_prime) {
    ASSERT_EQ(make_pair(1, 1), squares_r_prime(2));

    ASSERT_EQ(make_pair(1, 2), squares_r_prime(5));  // 1 (mod 4)
    ASSERT_EQ(make_pair(2, 3), squares_r_prime(13)); // 1 (mod 4)
    ASSERT_EQ(make_pair(1, 4), squares_r_prime(17)); // 1 (mod 4)
    ASSERT_EQ(make_pair(2, 5), squares_r_prime(29)); // 1 (mod 4)

    ASSERT_EQ(make_pair(0, 4), squares_r_prime(23)); // 3 (mod 4), no solution
}

TEST(squares_r_test, squares_r_prime_table) {
    unordered_map<int, pair<int, int>> expected{ { 2, {1, 1} }, { 5, {1, 2} }, { 8, {2, 2} }, { 10, {1, 3} }, { 13, {2, 3} }, { 17, {1, 4} }, {18, { 3, 3 }}, { 20, {2, 4} } };
    EXPECT_EQ(expected, squares_r_prime_table(20));
}

TEST(squares_r_test, squares_r_list) {
    for (int i = 1; i < 1000; i++) {
        auto vf = factor_integer_slow(i);
        for (int u = 0; u < 2; u++) {
            auto vr = squares_r_list(vf, u == 1);
            int r = squares_r(vf, u == 1);
            ASSERT_EQ(r, vr.size()) << "ERROR: " << i << " " << r << " " << testing::PrintToString(vr) << endl;
            for (const auto& t : vr) {
                ASSERT_EQ(i, sqT(t.first) + sqT(t.second)) << "ERROR: " << i << " " << testing::PrintToString(t) << endl;
            }
        }
    }
}

TEST(squares_r_test, squares_r_list_tbl_max) {
    auto tbl = squares_r_prime_table(100);
    for (int i = 1; i < 1000; i++) {
        int b_max = i / 4 + 1;
        auto vf = factor_integer_slow(i);
        for (int u = 0; u < 2; u++) {
            auto vr = squares_r_list(vf, u == 1, tbl);
            int r = squares_r(vf, u == 1);
            ASSERT_EQ(r, vr.size()) << "ERROR: " << i << " " << r << " " << testing::PrintToString(vr) << endl;
            vector<pair<int, int>> evr4;
            for (const auto& t : vr) {
                if (t.second <= b_max) evr4.push_back(t);
                ASSERT_EQ(i, sqT(t.first) + sqT(t.second)) << "ERROR: " << i << " " << testing::PrintToString(t) << endl;
            }
            EXPECT_EQ(evr4, squares_r_list(vf, u == 1, tbl, b_max));
        }
    }
}
