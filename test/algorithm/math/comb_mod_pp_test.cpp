#include "altruct/algorithm/math/comb_mod_pp.h"
#include "altruct/algorithm/math/ranges.h"
#include "altruct/algorithm/collections/collections.h"

#include "gtest/gtest.h"

#include <ostream>
#include <vector>

using namespace altruct::math;

template<typename T, int ID, int STORAGE_TYPE>
std::ostream& operator << (std::ostream& os, const modulo<T, ID, STORAGE_TYPE>& rhs) {
    return os << rhs.v;
}

TEST(comb_mod_pp_test, factorials_mod_pp_skipped) {
    EXPECT_EQ((std::vector<int64_t>{
        1, 1, 1, 3, 3, 15, 15, 105, 105, 945, 945, 10395, 10395, 135135, \
        135135, 2027025, 2027025, 34459425, 34459425, 654729075, 654729075, \
        13749310575, 13749310575, 316234143225, 316234143225, 7905853580625, \
        7905853580625, 213458046676875, 213458046676875, 560783819416255, \
        560783819416255, 495799799264545, 495799799264545, 598794679933249, \
        598794679933249, 691615474496483, 691615474496483, 819974605832143, \
        819974605832143, 453812235860105 }), factorials_mod_pp_skipped<int64_t>(39, 2, 50));

    EXPECT_EQ((std::vector<int64_t>{
        1, 1, 2, 2, 8, 40, 40, 280, 2240, 2240, 22400, 246400, 246400, \
        3203200, 44844800, 44844800, 717516800, 12197785600, 12197785600, \
        231757926400, 4635158528000, 4635158528000, 101973487616000, \
        80587762126861, 80587762126861, 161673864319684, 85697830418804, \
        85697830418804, 134736798685373, 201326784172135, 201326784172135, \
        64396346496715, 1771766948390, 1771766948390, 60240076245260, \
        49491347637610, 49491347637610, 184050805834378, 199523262582947, \
        199523262582947}), factorials_mod_pp_skipped<int64_t>(39, 3, 30));

    EXPECT_EQ((std::vector<int64_t>{
        1, 1, 2, 6, 24, 24, 144, 1008, 8064, 72576, 72576, 798336, 9580032, \
        124540416, 1743565824, 1743565824, 27897053184, 474249904128, \
        8536498274304, 66826035571151, 66826035571151, 68202704025421, \
        69948013949887, 82925414597401, 82861317525124, 82861317525124, \
        56310759559474, 89879033496423, 37059715243594, 25689994017351, \
        25689994017351, 33450361412881, 21369817165317, 37631944971086, \
        39709517688799, 39709517688799, 94398593828014, 59520432574018, \
        68325510078309, 89774238757176}), factorials_mod_pp_skipped<int64_t>(39, 5, 20));
}

TEST(comb_mod_pp_test, factorial_mod_pp_skipped_all) {
    int e_max = 5, n_max = 100;
    for (auto p : { 2, 3, 5, 7 }) {
        auto tbl = factorials_mod_pp_skipped<int64_t>(std::max(powT(p, e_max), n_max), p, e_max);
        for (int e = 1; e <= e_max; e++) {
            for (int n = 0; n <= n_max; n++) {
                auto r0 = std::pair<int64_t, int64_t>{ tbl[n] % powT(p, e), n / p };
                auto r1 = factorial_mod_pp_skipped_slow<int64_t>(int64_t(n), p, e);
                auto r3 = factorial_mod_pp_skipped<int64_t>(int64_t(n), p, e, tbl.data());
                EXPECT_EQ(r0, r1) << "factorial_mod_pp_skipped_slow(" << n << ", " << p << ", " << e << ")";
                EXPECT_EQ(r0, r3) << "factorial_mod_pp_skipped(" << n << ", " << p << ", " << e << ")";
            }
        }
    }
}

TEST(comb_mod_pp_test, factorial_mod_pp_reduced_all) {
    int e_max = 5, n_max = 100;
    for (auto p : { 2, 3, 5, 7 }) {
        auto tbl = factorials_mod_pp_skipped<int64_t>(std::max(powT(p, e_max), n_max), p, e_max);
        for (int e = 1; e <= e_max; e++) {
            for (int n = 0; n <= n_max; n++) {
                int64_t r = 1, a = 0, m = powT(p, e);
                for (int t = n; t > 0; t /= p) r = modulo_mul(r, tbl[t], m), a += t / p;
                auto r0 = std::make_pair(r, a);
                auto r1 = factorial_mod_pp_reduced_slow<int64_t>(int64_t(n), p, e);
                auto r2 = factorial_mod_pp_reduced_2<int64_t>(int64_t(n), p, e, tbl.data());
                auto r3 = factorial_mod_pp_reduced<int64_t>(int64_t(n), p, e, tbl.data());
                EXPECT_EQ(r0, r1) << "factorial_mod_pp_reduced_slow(" << n << ", " << p << ", " << e << ")";
                EXPECT_EQ(r0, r2) << "factorial_mod_pp_reduced_2(" << n << ", " << p << ", " << e << ")";
                EXPECT_EQ(r0, r3) << "factorial_mod_pp_reduced(" << n << ", " << p << ", " << e << ")";
            }
        }
    }
}

TEST(comb_mod_pp_test, binomial_mod_pp_reduced_all) {
    int e_max = 5, n_max = 50;
    for (auto p : { 2, 3, 5, 7 }) {
        auto tbl = factorials_mod_pp_skipped<int64_t>(std::max(powT(p, e_max), n_max), p, e_max);
        for (int e = 1; e <= e_max; e++) {
            for (int n = 0; n <= n_max; n++) {
                for (int k = 0; k <= (n + 1) / 2; k++) {
                    int64_t r = 1, a = 0, m = powT(p, e);
                    for (int t = n; t > 0; t /= p) r = modulo_mul(r, tbl[t], m), a += t / p;
                    for (int t = k; t > 0; t /= p) r = modulo_div(r, tbl[t], m), a -= t / p;
                    for (int t = n - k; t > 0; t /= p) r = modulo_div(r, tbl[t], m), a -= t / p;
                    auto r0 = std::make_pair(r, a);
                    auto r1 = binomial_mod_pp_reduced_slow<int64_t>(int64_t(n), int64_t(k), p, e);
                    auto r2 = binomial_mod_pp_reduced_2<int64_t>(int64_t(n), int64_t(k), p, e, tbl.data());
                    auto r3 = binomial_mod_pp_reduced<int64_t>(int64_t(n), int64_t(k), p, e, tbl.data());
                    EXPECT_EQ(r0, r1) << "binomial_mod_pp_reduced_slow(" << n << ", " << k << ", " << p << ", " << e << ")";
                    EXPECT_EQ(r0, r2) << "binomial_mod_pp_reduced_2(" << n << ", " << k << ", " << p << ", " << e << ")";
                    EXPECT_EQ(r0, r3) << "binomial_mod_pp_reduced(" << n << ", " << k << ", " << p << ", " << e << ")";
                }
            }
        }
    }
}

TEST(comb_mod_pp_test, reduced_2_p37e1) {
    using std::make_pair;
    const int P = 37, E = 1, PE = powT(P, E);
    auto tbl = factorials_mod_pp_skipped<int>(PE, P, E);
    EXPECT_EQ(make_pair(25, 0LL), factorial_mod_pp_reduced_2(10LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(31, 2LL), factorial_mod_pp_reduced_2(100LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(7, 27LL), factorial_mod_pp_reduced_2(1000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(19, 277LL), factorial_mod_pp_reduced_2(10000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(3, 2776LL), factorial_mod_pp_reduced_2(100000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(30, 27776LL), factorial_mod_pp_reduced_2(1000000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(25, 2777777774LL), factorial_mod_pp_reduced_2(100000000000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(26, 2777777777777770LL), factorial_mod_pp_reduced_2(100000000000000000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(20, 2LL), binomial_mod_pp_reduced_2(1000000LL, 1234LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(32, 4LL), binomial_mod_pp_reduced_2(100000000000LL, 12345678LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(4, 3LL), binomial_mod_pp_reduced_2(100000000000000000LL, 1234567891235LL, P, E, tbl.data()));
}

TEST(comb_mod_pp_test, reduced_2_p3e4) {
    using std::make_pair;
    const int P = 3, E = 4, PE = 81;
    auto tbl = factorials_mod_pp_skipped<int>(PE, P, E);
    EXPECT_EQ(make_pair(7, 4LL), factorial_mod_pp_reduced_2(10LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(65, 48LL), factorial_mod_pp_reduced_2(100LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(31, 498LL), factorial_mod_pp_reduced_2(1000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(5, 4996LL), factorial_mod_pp_reduced_2(10000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(22, 49995LL), factorial_mod_pp_reduced_2(100000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(58, 499993LL), factorial_mod_pp_reduced_2(1000000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(64, 49999999991LL), factorial_mod_pp_reduced_2(100000000000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(2, 49999999999999978LL), factorial_mod_pp_reduced_2(100000000000000000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(17, 2LL), binomial_mod_pp_reduced_2(1000000LL, 1234LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(59, 6LL), binomial_mod_pp_reduced_2(100000000000LL, 12345678LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(38, 10LL), binomial_mod_pp_reduced_2(100000000000000000LL, 1234567891235LL, P, E, tbl.data()));
}

TEST(comb_mod_pp_test, reduced_p37e1) {
    using std::make_pair;
    const int P = 37, E = 1;
    auto tbl = factorials_mod_pp_skipped<int>(P * E, P, E);
    EXPECT_EQ(make_pair(25, 0LL), factorial_mod_pp_reduced(10LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(31, 2LL), factorial_mod_pp_reduced(100LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(7, 27LL), factorial_mod_pp_reduced(1000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(19, 277LL), factorial_mod_pp_reduced(10000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(3, 2776LL), factorial_mod_pp_reduced(100000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(30, 27776LL), factorial_mod_pp_reduced(1000000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(25, 2777777774LL), factorial_mod_pp_reduced(100000000000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(26, 2777777777777770LL), factorial_mod_pp_reduced(100000000000000000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(20, 2LL), binomial_mod_pp_reduced(1000000LL, 1234LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(32, 4LL), binomial_mod_pp_reduced(100000000000LL, 12345678LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(4, 3LL), binomial_mod_pp_reduced(100000000000000000LL, 1234567891235LL, P, E, tbl.data()));
}

TEST(comb_mod_pp_test, reduced_p3e4) {
    using std::make_pair;
    const int P = 3, E = 4;
    auto tbl = factorials_mod_pp_skipped<int>(P * E, P, E);
    EXPECT_EQ(make_pair(7, 4LL), factorial_mod_pp_reduced(10LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(65, 48LL), factorial_mod_pp_reduced(100LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(31, 498LL), factorial_mod_pp_reduced(1000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(5, 4996LL), factorial_mod_pp_reduced(10000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(22, 49995LL), factorial_mod_pp_reduced(100000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(58, 499993LL), factorial_mod_pp_reduced(1000000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(64, 49999999991LL), factorial_mod_pp_reduced(100000000000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(2, 49999999999999978LL), factorial_mod_pp_reduced(100000000000000000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(17, 2LL), binomial_mod_pp_reduced(1000000LL, 1234LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(59, 6LL), binomial_mod_pp_reduced(100000000000LL, 12345678LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(38, 10LL), binomial_mod_pp_reduced(100000000000000000LL, 1234567891235LL, P, E, tbl.data()));
}

TEST(comb_mod_pp_test, reduced_p37e4) {
    using std::make_pair;
    const int P = 37, E = 4;
    auto tbl = factorials_mod_pp_skipped<int>(P * E, P, E);
    EXPECT_EQ(make_pair(1754639, 0LL), factorial_mod_pp_reduced(10LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(1404181, 2LL), factorial_mod_pp_reduced(100LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(832174, 27LL), factorial_mod_pp_reduced(1000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(900155, 277LL), factorial_mod_pp_reduced(10000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(1176936, 2776LL), factorial_mod_pp_reduced(100000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(1579560, 27776LL), factorial_mod_pp_reduced(1000000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(333765, 2777777774LL), factorial_mod_pp_reduced(100000000000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(410504, 2777777777777770LL), factorial_mod_pp_reduced(100000000000000000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(1067285, 2LL), binomial_mod_pp_reduced(1000000LL, 1234LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(628588, 4LL), binomial_mod_pp_reduced(100000000000LL, 12345678LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(600514, 3LL), binomial_mod_pp_reduced(100000000000000000LL, 1234567891235LL, P, E, tbl.data()));
}
