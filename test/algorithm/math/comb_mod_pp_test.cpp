#include "algorithm/math/comb_mod_pp.h"
#include "algorithm/math/ranges.h"
#include "algorithm/collections/collections.h"

#include "gtest/gtest.h"

#include <ostream>
#include <vector>

using namespace std;
using namespace altruct::math;

template<typename T, int ID, int STORAGE_TYPE>
std::ostream& operator << (std::ostream& os, const modulo<T, ID, STORAGE_TYPE>& rhs) {
    return os << rhs.v;
}

TEST(comb_mod_pp_test, factorials_mod_pp_skipped) {
    EXPECT_EQ((vector<int64_t>{
        1, 1, 1, 3, 3, 15, 15, 105, 105, 945, 945, 10395, 10395, 135135, \
        135135, 2027025, 2027025, 34459425, 34459425, 654729075, 654729075, \
        13749310575, 13749310575, 316234143225, 316234143225, 7905853580625, \
        7905853580625, 213458046676875, 213458046676875, 560783819416255, \
        560783819416255, 495799799264545, 495799799264545, 598794679933249, \
        598794679933249, 691615474496483, 691615474496483, 819974605832143, \
        819974605832143, 453812235860105 }), factorials_mod_pp_skipped<int64_t>(39, 2, 50));

    EXPECT_EQ((vector<int64_t>{
        1, 1, 2, 2, 8, 40, 40, 280, 2240, 2240, 22400, 246400, 246400, \
        3203200, 44844800, 44844800, 717516800, 12197785600, 12197785600, \
        231757926400, 4635158528000, 4635158528000, 101973487616000, \
        80587762126861, 80587762126861, 161673864319684, 85697830418804, \
        85697830418804, 134736798685373, 201326784172135, 201326784172135, \
        64396346496715, 1771766948390, 1771766948390, 60240076245260, \
        49491347637610, 49491347637610, 184050805834378, 199523262582947, \
        199523262582947}), factorials_mod_pp_skipped<int64_t>(39, 3, 30));

    EXPECT_EQ((vector<int64_t>{
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
        auto tbl = factorials_mod_pp_skipped<int64_t>(max(powT(p, e_max), n_max), p, e_max);
        for (int e = 1; e <= e_max; e++) {
            for (int n = 0; n <= n_max; n++) {
                int64_t r0 = tbl[n] % powT(p, e);
                int64_t r1 = factorial_mod_pp_skipped_slow<int64_t>(int64_t(n), p, e);
                int64_t r3 = factorial_mod_pp_skipped<int64_t>(int64_t(n), p, e, tbl.data());
                EXPECT_EQ(r0, r1) << "factorial_mod_pp_skipped_slow(" << n << ", " << p << ", " << e << ")";
                EXPECT_EQ(r0, r3) << "factorial_mod_pp_skipped(" << n << ", " << p << ", " << e << ")";
            }
        }
    }
}

TEST(comb_mod_pp_test, factorial_mod_pp_reduced_all) {
    int e_max = 5, n_max = 100;
    for (auto p : { 2, 3, 5, 7 }) {
        auto tbl = factorials_mod_pp_skipped<int64_t>(max(powT(p, e_max), n_max), p, e_max);
        for (int e = 1; e <= e_max; e++) {
            for (int n = 0; n <= n_max; n++) {
                int64_t r0 = 1; for (int t = n; t > 0; t /= p) r0 = r0 * tbl[t] % powT(p, e);
                int64_t r1 = factorial_mod_pp_reduced_slow<int64_t>(int64_t(n), p, e);
                int64_t r2 = factorial_mod_pp_reduced_2<int64_t>(int64_t(n), p, e, tbl.data()).first.v;
                int64_t r3 = factorial_mod_pp_reduced<int64_t>(int64_t(n), p, e, tbl.data());
                EXPECT_EQ(r0, r1) << "factorial_mod_pp_reduced_slow(" << n << ", " << p << ", " << e << ")";
                EXPECT_EQ(r0, r3) << "factorial_mod_pp_reduced_2(" << n << ", " << p << ", " << e << ")";
                EXPECT_EQ(r0, r3) << "factorial_mod_pp_reduced(" << n << ", " << p << ", " << e << ")";
            }
        }
    }
}

TEST(comb_mod_pp_test, binomial_mod_pp_reduced_all) {
    int e_max = 5, n_max = 50;
    for (auto p : { 2, 3, 5, 7 }) {
        auto tbl = factorials_mod_pp_skipped<int64_t>(max(powT(p, e_max), n_max), p, e_max);
        for (int e = 1; e <= e_max; e++) {
            for (int n = 0; n <= n_max; n++) {
                for (int k = 0; k <= (n + 1) / 2; k++) {
                    int64_t r0 = 1; for (int t = n; t > 0; t /= p) r0 = r0 * tbl[t] % powT(p, e);
                    int64_t s0 = 1; for (int t = k; t > 0; t /= p) s0 = s0 * tbl[t] % powT(p, e);
                    int l = n - k; for (int t = l; t > 0; t /= p) s0 = s0 * tbl[t] % powT(p, e);
                    r0 = modulo_div<int64_t>(r0, s0, powT(p, e));
                    int64_t r1 = binomial_mod_pp_reduced_slow<int64_t>(int64_t(n), int64_t(k), p, e);
                    int64_t r2 = binomial_mod_pp_reduced_2<int64_t>(int64_t(n), int64_t(k), p, e, tbl.data()).first.v;
                    int64_t r3 = binomial_mod_pp_reduced<int64_t>(int64_t(n), int64_t(k), p, e, tbl.data());
                    EXPECT_EQ(r0, r1) << "binomial_mod_pp_reduced_slow(" << n << ", " << k << ", " << p << ", " << e << ")";
                    EXPECT_EQ(r0, r3) << "binomial_mod_pp_reduced_2(" << n << ", " << k << ", " << p << ", " << e << ")";
                    EXPECT_EQ(r0, r3) << "binomial_mod_pp_reduced(" << n << ", " << k << ", " << p << ", " << e << ")";
                }
            }
        }
    }
}

TEST(comb_mod_pp_test, reduced_2_p37e1) {
    typedef moduloX<int> modx;
    const int P = 37, E = 1, PE = powT(P, E);
    auto tbl = factorials_mod_pp_skipped<int>(PE, P, E);
    EXPECT_EQ(make_pair(modx(25, PE), 0LL), factorial_mod_pp_reduced_2(10LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(modx(31, PE), 2LL), factorial_mod_pp_reduced_2(100LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(modx(7, PE), 27LL), factorial_mod_pp_reduced_2(1000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(modx(19, PE), 277LL), factorial_mod_pp_reduced_2(10000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(modx(3, PE), 2776LL), factorial_mod_pp_reduced_2(100000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(modx(30, PE), 27776LL), factorial_mod_pp_reduced_2(1000000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(modx(25, PE), 2777777774LL), factorial_mod_pp_reduced_2(100000000000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(modx(26, PE), 2777777777777770LL), factorial_mod_pp_reduced_2(100000000000000000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(modx(20, PE), 2LL), binomial_mod_pp_reduced_2(1000000LL, 1234LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(modx(32, PE), 4LL), binomial_mod_pp_reduced_2(100000000000LL, 12345678LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(modx(4, PE), 3LL), binomial_mod_pp_reduced_2(100000000000000000LL, 1234567891235LL, P, E, tbl.data()));
}

TEST(comb_mod_pp_test, reduced_2_p3e4) {
    typedef moduloX<int> modx;
    const int P = 3, E = 4, PE = 81;
    auto tbl = factorials_mod_pp_skipped<int>(PE, P, E);
    EXPECT_EQ(make_pair(modx(7, PE), 4LL), factorial_mod_pp_reduced_2(10LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(modx(65, PE), 48LL), factorial_mod_pp_reduced_2(100LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(modx(31, PE), 498LL), factorial_mod_pp_reduced_2(1000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(modx(5, PE), 4996LL), factorial_mod_pp_reduced_2(10000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(modx(22, PE), 49995LL), factorial_mod_pp_reduced_2(100000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(modx(58, PE), 499993LL), factorial_mod_pp_reduced_2(1000000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(modx(64, PE), 49999999991LL), factorial_mod_pp_reduced_2(100000000000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(modx(2, PE), 49999999999999978LL), factorial_mod_pp_reduced_2(100000000000000000LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(modx(17, PE), 2LL), binomial_mod_pp_reduced_2(1000000LL, 1234LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(modx(59, PE), 6LL), binomial_mod_pp_reduced_2(100000000000LL, 12345678LL, P, E, tbl.data()));
    EXPECT_EQ(make_pair(modx(38, PE), 10LL), binomial_mod_pp_reduced_2(100000000000000000LL, 1234567891235LL, P, E, tbl.data()));
}

TEST(comb_mod_pp_test, reduced_p37e1) {
    const int P = 37, E = 1;
    auto tbl = factorials_mod_pp_skipped<int>(P * E, P, E);
    EXPECT_EQ(25, factorial_mod_pp_reduced(10LL, P, E, tbl.data()));
    EXPECT_EQ(31, factorial_mod_pp_reduced(100LL, P, E, tbl.data()));
    EXPECT_EQ(7, factorial_mod_pp_reduced(1000LL, P, E, tbl.data()));
    EXPECT_EQ(19, factorial_mod_pp_reduced(10000LL, P, E, tbl.data()));
    EXPECT_EQ(3, factorial_mod_pp_reduced(100000LL, P, E, tbl.data()));
    EXPECT_EQ(30, factorial_mod_pp_reduced(1000000LL, P, E, tbl.data()));
    EXPECT_EQ(25, factorial_mod_pp_reduced(100000000000LL, P, E, tbl.data()));
    EXPECT_EQ(26, factorial_mod_pp_reduced(100000000000000000LL, P, E, tbl.data()));
    EXPECT_EQ(20, binomial_mod_pp_reduced(1000000LL, 1234LL, P, E, tbl.data()));
    EXPECT_EQ(32, binomial_mod_pp_reduced(100000000000LL, 12345678LL, P, E, tbl.data()));
    EXPECT_EQ(4, binomial_mod_pp_reduced(100000000000000000LL, 1234567891235LL, P, E, tbl.data()));
}

TEST(comb_mod_pp_test, reduced_p3e4) {
    const int P = 3, E = 4;
    auto tbl = factorials_mod_pp_skipped<int>(P * E, P, E);
    EXPECT_EQ(7, factorial_mod_pp_reduced(10LL, P, E, tbl.data()));
    EXPECT_EQ(65, factorial_mod_pp_reduced(100LL, P, E, tbl.data()));
    EXPECT_EQ(31, factorial_mod_pp_reduced(1000LL, P, E, tbl.data()));
    EXPECT_EQ(5, factorial_mod_pp_reduced(10000LL, P, E, tbl.data()));
    EXPECT_EQ(22, factorial_mod_pp_reduced(100000LL, P, E, tbl.data()));
    EXPECT_EQ(58, factorial_mod_pp_reduced(1000000LL, P, E, tbl.data()));
    EXPECT_EQ(64, factorial_mod_pp_reduced(100000000000LL, P, E, tbl.data()));
    EXPECT_EQ(2, factorial_mod_pp_reduced(100000000000000000LL, P, E, tbl.data()));
    EXPECT_EQ(17, binomial_mod_pp_reduced(1000000LL, 1234LL, P, E, tbl.data()));
    EXPECT_EQ(59, binomial_mod_pp_reduced(100000000000LL, 12345678LL, P, E, tbl.data()));
    EXPECT_EQ(38, binomial_mod_pp_reduced(100000000000000000LL, 1234567891235LL, P, E, tbl.data()));
}

TEST(comb_mod_pp_test, reduced_p37e4) {
    const int P = 37, E = 4;
    auto tbl = factorials_mod_pp_skipped<int>(P * E, P, E);
    EXPECT_EQ(1754639, factorial_mod_pp_reduced(10LL, P, E, tbl.data()));
    EXPECT_EQ(1404181, factorial_mod_pp_reduced(100LL, P, E, tbl.data()));
    EXPECT_EQ(832174, factorial_mod_pp_reduced(1000LL, P, E, tbl.data()));
    EXPECT_EQ(900155, factorial_mod_pp_reduced(10000LL, P, E, tbl.data()));
    EXPECT_EQ(1176936, factorial_mod_pp_reduced(100000LL, P, E, tbl.data()));
    EXPECT_EQ(1579560, factorial_mod_pp_reduced(1000000LL, P, E, tbl.data()));
    EXPECT_EQ(333765, factorial_mod_pp_reduced(100000000000LL, P, E, tbl.data()));
    EXPECT_EQ(410504, factorial_mod_pp_reduced(100000000000000000LL, P, E, tbl.data()));
    EXPECT_EQ(1067285, binomial_mod_pp_reduced(1000000LL, 1234LL, P, E, tbl.data()));
    EXPECT_EQ(628588, binomial_mod_pp_reduced(100000000000LL, 12345678LL, P, E, tbl.data()));
    EXPECT_EQ(600514, binomial_mod_pp_reduced(100000000000000000LL, 1234567891235LL, P, E, tbl.data()));
}
