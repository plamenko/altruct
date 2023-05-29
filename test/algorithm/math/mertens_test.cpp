#include "altruct/algorithm/collections/collections.h"
#include "altruct/algorithm/math/mertens.h"
#include "altruct/algorithm/math/primes.h"
#include "altruct/algorithm/math/ranges.h"
#include "altruct/structure/math/modulo.h"

#include "gtest/gtest.h"

#include <functional>

using namespace std;
using namespace altruct::math;
using namespace altruct::container;
using namespace altruct::collections;

namespace {
typedef moduloX<int> modx;

vector<int> primes_table(int n) {
    vector<int> p(n);
    int m = primes(p.data(), nullptr, n);
    p.resize(m);
    return p;
}

auto prime_pi_table(int n) {
    auto pa = primes_table(n + 1);
    return make_sqrt_map<int64_t, int>([&](int64_t n) {
        return int(std::upper_bound(pa.begin(), pa.end(), n) - pa.begin());
    }, n);
}

vector<moduloX<int>> to_modx(int M, const vector<int>& v) {
    return transform(v, [&](int a) {
        return moduloX<int>(a, M);
    });
}
} // namespace

TEST(mertens_test, sieve_mertens) {
    int n = 31;
    auto pa = primes_table(n);
    auto expected1 = vector<int>{ 0, 1, 0, -1, -1, -2, -1, -2, -2, -2, -1, -2, -2, -3, -2, -1, -1, -2, -2, -3, -3, -2, -1, -2, -2, -2, -1, -1, -1, -2, -3 };
    vector<int> actual1(n); sieve_mertens(actual1, n, pa.data(), (int)pa.size(), int(1));
    EXPECT_EQ(expected1, actual1);

    typedef moduloX<int> modx;
    auto expected2 = to_modx(1009, expected1);
    vector<modx> actual2(n); sieve_mertens(actual2, n, pa.data(), (int)pa.size(), modx(1, 1009));
    EXPECT_EQ(expected2, actual2);
}

TEST(mertens_test, sieve_mertens_odd) {
    int n = 31;
    auto pa = primes_table(n);
    auto expected1 = vector<int>{ 0, 1, 1, 0, 0, -1, -1, -2, -2, -2, -2, -3, -3, -4, -4, -3, -3, -4, -4, -5, -5, -4, -4, -5, -5, -5, -5, -5, -5, -6, -6 };
    vector<int> actual1(n); sieve_mertens_odd(actual1, n, pa.data(), (int)pa.size(), int(1));
    EXPECT_EQ(expected1, actual1);

    typedef moduloX<int> modx;
    auto expected2 = to_modx(1009, expected1);
    vector<modx> actual2(n); sieve_mertens_odd(actual2, n, pa.data(), (int)pa.size(), modx(1, 1009));
    EXPECT_EQ(expected2, actual2);
}

TEST(mertens_test, sieve_mertens_even) {
    int n = 31;
    auto expected1 = vector<int>{ 0, 0, -1, -1, -1, -1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 3 };
    vector<int> actual1(n); sieve_mertens_even(actual1, n, int(1));
    EXPECT_EQ(expected1, actual1);

    typedef moduloX<int> modx;
    auto expected2 = to_modx(1009, expected1);
    vector<modx> actual2(n); sieve_mertens_even(actual2, n, modx(1, 1009));
    EXPECT_EQ(expected2, actual2);
}

TEST(mertens_test, sieve_mertens_even_odd) {
    int n = 31;
    auto pa = primes_table(n);
    auto expected_even1 = vector<int>{ 0, 0, -1, -1, -1, -1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 3 };
    auto expected_odd1 = vector<int>{ 0, 1, 1, 0, 0, -1, -1, -2, -2, -2, -2, -3, -3, -4, -4, -3, -3, -4, -4, -5, -5, -4, -4, -5, -5, -5, -5, -5, -5, -6, -6 };
    vector<int> actual_even1(n), actual_odd1(n);
    sieve_mertens_even_odd(actual_even1, actual_odd1, n, pa.data(), (int)pa.size(), int(1));
    EXPECT_EQ(expected_even1, actual_even1);
    EXPECT_EQ(expected_odd1, actual_odd1);

    typedef moduloX<int> modx;
    auto expected_even2 = to_modx(1009, expected_even1);
    auto expected_odd2 = to_modx(1009, expected_odd1);
    vector<modx> actual_even2(n), actual_odd2(n);
    sieve_mertens_even_odd(actual_even2, actual_odd2, n, pa.data(), (int)pa.size(), modx(1, 1009));
    EXPECT_EQ(expected_even2, actual_even2);
    EXPECT_EQ(expected_odd2, actual_odd2);
}

TEST(mertens_test, mertens) {
    auto v_M = to_modx(1009, { 0, 1, 0, -1, -1, -2, -1, -2, -2, -2, -1, -2, -2, -3, -2, -1, -1, -2, -2, -3, -3, -2, -1, -2, -2, -2, -1, -1, -1, -2, -3 });
    for (int n = 1; n < v_M.size(); n++) {
        // preprocess `U = n^(2/3)` values
        int U = int(isq(icbrt(n)));
        sqrt_map<int, modx> M(U, n);
        for (int i = 0; i < U; i++) M[i] = v_M[i];
        // call mertens
        mertens(n, M, modx(1, 1009));
        for (int i = 1; i <= n; i++) {
            EXPECT_EQ(v_M[n / i], M[n / i]) << "n: " << n << ", i: " << i;
        }
    }
}

TEST(mertens_test, mertens_pi) {
    auto v_M = to_modx(1009, { 0, 1, 0, -1, -1, -2, -1, -2, -2, -2, -1, -2, -2, -3, -2, -1, -1, -2, -2, -3, -3, -2, -1, -2, -2, -2, -1, -1, -1, -2, -3 });
    for (int n = 1; n < v_M.size(); n++) {
        // preprocess
        int u = max(3, int(sqrt(n * log(n))));
        auto pa = primes_table(u);
        auto pi_tbl = prime_pi_table(n);
        // call mertens
        auto M = mertens(n, pi_tbl, pa.data(), (int)pa.size(), modx(1, 1009));
        for (int i = 1; i <= n; i++) {
            EXPECT_EQ(v_M[n / i], M[n / i]) << "n: " << n << ", i: " << i;
        }
    }
}

TEST(mertens_test, mertens_odd) {
    auto v_M1 = to_modx(1009, { 0, 1, 1, 0, 0, -1, -1, -2, -2, -2, -2, -3, -3, -4, -4, -3, -3, -4, -4, -5, -5, -4, -4, -5, -5, -5, -5, -5, -5, -6, -6 });
    for (int n = 1; n < v_M1.size(); n++) {
        // preprocess `U = n^(2/3)` values
        int U = int(isq(icbrt(n)));
        sqrt_map<int, modx> M1(U, n);
        for (int i = 0; i < U; i++) M1[i] = v_M1[i];
        // call mertens_odd
        mertens_odd(n, M1, modx(1, 1009));
        for (int i = 1; i <= n; i++) {
            EXPECT_EQ(v_M1[n / i], M1[n / i]) << "n: " << n << ", i: " << i;
        }
    }
}

TEST(mertens_test, mertens_odd_pi) {
    auto v_M1 = to_modx(1009, { 0, 1, 1, 0, 0, -1, -1, -2, -2, -2, -2, -3, -3, -4, -4, -3, -3, -4, -4, -5, -5, -4, -4, -5, -5, -5, -5, -5, -5, -6, -6 });
    for (int n = 1; n < v_M1.size(); n++) {
        // preprocess
        int u = max(3, int(sqrt(n * log(n))));
        auto pa = primes_table(u);
        auto pi_tbl = prime_pi_table(n);
        // call mertens_odd
        auto M1 = mertens_odd(n, pi_tbl, pa.data(), (int)pa.size(), modx(1, 1009));
        for (int i = 1; i <= n; i++) {
            EXPECT_EQ(v_M1[n / i], M1[n / i]) << "n: " << n << ", i: " << i;
        }
    }
}

TEST(mertens_test, mertens_even) {
    auto v_M0 = to_modx(1009, { 0, 0, -1, -1, -1, -1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 3 });
    for (int n = 1; n < v_M0.size(); n++) {
        // preprocess `U = n^(2/3)` values
        int U = int(isq(icbrt(n)));
        sqrt_map<int, modx> M0(U, n);
        for (int i = 0; i < U; i++) M0[i] = v_M0[i];
        // call mertens_even
        mertens_even(n, M0, modx(1, 1009));
        for (int i = 1; i <= n; i++) {
            EXPECT_EQ(v_M0[n / i], M0[n / i]) << "n: " << n << ", i: " << i;
        }
    }
}

TEST(mertens_test, mertens_even_pi) {
    auto v_M0 = to_modx(1009, { 0, 0, -1, -1, -1, -1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 3 });
    for (int n = 1; n < v_M0.size(); n++) {
        // preprocess
        int u = max(3, int(sqrt(n * log(n))));
        auto pa = primes_table(u);
        auto pi_tbl = prime_pi_table(n);
        // call mertens_even
        auto M0 = mertens_even(n, pi_tbl, pa.data(), (int)pa.size(), modx(1, 1009));
        for (int i = 1; i <= n; i++) {
            EXPECT_EQ(v_M0[n / i], M0[n / i]) << "n: " << n << ", i: " << i;
        }
    }
}
