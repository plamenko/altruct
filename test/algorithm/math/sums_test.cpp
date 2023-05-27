#include "altruct/algorithm/math/sums.h"
#include "altruct/algorithm/math/primes.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

typedef modulo<int, 1000000007> field;

TEST(sums_test, sum_ratio) {
    auto f = [](int a, int b, int q, int n) {
        int s = 0;
        for (int k = 0; k < n; k++) {
            s += div_floor(a * k + b, q);
        }
        return s;
    };
    int U = 20;
    for (int a = -U; a < U; a++) {
        for (int b = -U; b < U; b++) {
            for (int q = -U; q < U; q++) {
                if (q == 0) continue;
                for (int n = -3; n < U; n++) {
                    EXPECT_EQ(f(a, b, q, n), (sum_ratio<int>(a, b, q, n)));
                }
            }
        }
    }
}

TEST(sums_test, sum_ratio_modx) {
    using modx = moduloX<int>;
    auto f = [](int a, int b, int q, int n) {
        int s = 0;
        for (int k = 0; k < n; k++) {
            s += div_floor(a * k + b, q);
        }
        return modx(s, 101);
    };
    int U = 20;
    for (int a = -U; a < U; a++) {
        for (int b = -U; b < U; b++) {
            for (int q = -U; q < U; q++) {
                if (q == 0) continue;
                for (int n = -3; n < U; n++) {
                    EXPECT_EQ(f(a, b, q, n), (sum_ratio<modx>(a, b, q, n, modx(0, 101))));
                }
            }
        }
    }
}

TEST(sums_test, sum) {
    auto f = [](int k) { return k*k; };
    EXPECT_EQ(0, (sum<int>(f, 1, 0)));
    EXPECT_EQ(1, (sum<int>(f, 1, 1)));
    EXPECT_EQ(5, (sum<int>(f, 1, 2)));
    EXPECT_EQ(385, (sum<int>(f, 1, 10)));
    EXPECT_EQ(355, (sum<int>(f, 5, 10)));
}

vector<field> calc_sum_pow(int p, int n) {
    vector<field> v;
    for (int k = 0; k <= n; k++) {
        v.push_back(sum_pow<field>(p, k));
    }
    return v;

}
TEST(sums_test, sum_pow) {
    EXPECT_EQ((vector<field>{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }), calc_sum_pow(0, 10));
    EXPECT_EQ((vector<field>{ 0, 1, 3, 6, 10, 15, 21, 28, 36, 45, 55 }), calc_sum_pow(1, 10));
    EXPECT_EQ((vector<field>{ 0, 1, 5, 14, 30, 55, 91, 140, 204, 285, 385 }), calc_sum_pow(2, 10));
    EXPECT_EQ((vector<field>{ 0, 1, 9, 36, 100, 225, 441, 784, 1296, 2025, 3025 }), calc_sum_pow(3, 10));
    EXPECT_EQ((vector<field>{ 0, 1, 129, 2316, 18700, 96825, 376761, 1200304, 3297456, 8080425, 18080425 }), calc_sum_pow(7, 10));
}

vector<field> calc_sum_powx(int p, field x, int n) {
    vector<field> v;
    for (int k = 0; k <= n; k++) {
        v.push_back(sum_powx<field>(p, x, k));
    }
    return v;

}
TEST(sums_test, sum_powx) {
    EXPECT_EQ((vector<field>{0, 2, 6, 14, 30, 62, 126, 254, 510, 1022, 2046}), calc_sum_powx(0, 2, 10));
    EXPECT_EQ((vector<field>{0, 3, 21, 102, 426, 1641, 6015, 21324, 73812, 250959, 841449}), calc_sum_powx(1, 3, 10));
    EXPECT_EQ((vector<field>{0, 4, 68, 644, 4740, 30340, 177796, 980612, 5174916, 26408580, 131266180}), calc_sum_powx(2, 4, 10));
    EXPECT_EQ((vector<field>{0, 5, 205, 3580, 43580, 434205, 3809205, 30606080, 230606080, 654434198, 420059128}), calc_sum_powx(3, 5, 10));
    EXPECT_EQ((vector<field>{0, 6, 4614, 477006, 21710670, 629210670, 689904595, 229236226, 639265204, 946487221, 702254587}), calc_sum_powx(7, 6, 10));
}

TEST(sums_test, sum_sqrt) {
    // Sum[[n/k], {k,1,n}]
    auto f0 = [](int m) { return m; };
    vector<int> ve0, va0;
    for (int n = 0; n < 100; n++) {
        ve0.push_back(sum<int>([&](int k) { return f0(n / k); }, 1, n));
        va0.push_back(sum_sqrt<int>(f0, n));
    }
    EXPECT_EQ(ve0, va0);

    // Sum[k*[n/k], {k,1,n}]
    auto f1 = [](int k, int m) { return k * m; };
    auto sf1 = [](int k, int m) { return sum_pow<int>(1, k) * m; };
    vector<int> ve1, va1;
    for (int n = 0; n < 100; n++) {
        ve1.push_back(sum<int>([&](int k) { return f1(k, n / k); }, 1, n));
        va1.push_back(sum_sqrt2<int>(sf1, n));
    }
    EXPECT_EQ(ve1, va1);

    // Sum[k^2[n/k]^2, {k,1,n}]
    auto f2 = [](int k, int m) { return k * k * m * m; };
    auto sf2 = [](int k, int m) { return sum_pow<int>(2, k) * m * m; };
    vector<int> ve2, va2;
    for (int n = 0; n < 100; n++) {
        ve2.push_back(sum<int>([&](int k) { return f2(k, n / k); }, 1, n));
        va2.push_back(sum_sqrt2<int>(sf2, n));
    }
    EXPECT_EQ(ve2, va2);

    // Sum[(k+3)*[n/k+2]^2, {k,1,n}]
    auto f3 = [](int k) { return k + 3; };
    auto sf3 = [](int k) { return sum_pow<int>(1, k) + 3 * k; };
    auto g3 = [](int m) { return m + 2; };
    vector<int> ve3, va3a, va3b;
    for (int n = 0; n < 100; n++) {
        ve3.push_back(sum<int>([&](int k) { return f3(k) * g3(n / k); }, 1, n));
        va3a.push_back(sum_sqrt2m<int>(sf3, g3, n));
        va3b.push_back(sum_sqrt2m<int>(f3, sf3, g3, n, 0));
    }
    EXPECT_EQ(ve3, va3a);
    EXPECT_EQ(ve3, va3b);
}
