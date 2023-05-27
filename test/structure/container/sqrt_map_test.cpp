#include "altruct/algorithm/math/base.h"
#include "altruct/structure/container/sqrt_map.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;
using namespace altruct::container;

namespace {
template<typename F>
void assert_oor(F f, const char* expected_msg = nullptr) {
    try {
        f();
        FAIL();
    } catch (const std::out_of_range& e) {
        if (expected_msg) {
            EXPECT_EQ(e.what(), std::string(expected_msg));
        }
    } catch (...) {
        FAIL() << "Expected std::out_of_range";
    }
}
}

TEST(sqrt_map_test, sqrt_map) {
    int n = 100;
    int q = isqrt(n);
    sqrt_map<int, double> m(q, n);
    for (int i = 1; i <= q; i++) {
        EXPECT_EQ(0, m.count(i)) << "unexpected element at " << i;
        m[i] = 1.0 / i;
    }
    for (int i = 1; i <= q; i++) {
        int k = n / i;
        EXPECT_EQ(k <= q ? 1 : 0, m.count(k)) << "unexpected element at " << k;
        m[k] = 1.0 / k;
    }

    for (int i = 1; i <= q; i++) {
        EXPECT_EQ(1, m.count(i)) << "no element at " << i;
        EXPECT_EQ(1.0 / i, m.at(i)) << "unexpected element at " << i;
    }
    for (int i = 1; i <= q; i++) {
        int k = n / i;
        EXPECT_EQ(1, m.count(k)) << "no element at " << k;
        EXPECT_EQ(1.0 / k, m.at(k)) << "unexpected element at " << k;
    }

    n = 200;
    m.reset_max(n);
    for (int i = 1; i <= q; i++) {
        EXPECT_EQ(1, m.count(i)) << "no element at " << i;
        EXPECT_EQ(1.0 / i, m.at(i)) << "unexpected element at " << i;
    }
    for (int i = 1; i <= q; i++) {
        int k = n / i;
        EXPECT_EQ(k <= q ? 1 : 0, m.count(k)) << "unexpected element at " << k;
        m[k] = 1.0 / k;
    }

    n = 50;
    m.reset_max(n);
    for (int i = 1; i <= q; i++) {
        EXPECT_EQ(1, m.count(i)) << "no element at " << i;
        EXPECT_EQ(1.0 / i, m.at(i)) << "unexpected element at " << i;
    }
    for (int i = 1; i <= q; i++) {
        int k = n / i;
        EXPECT_EQ(k <= q ? 1 : 0, m.count(k)) << "unexpected element at " << k;
        m[k] = 1.0 / k;
    }
}

TEST(sqrt_map_test, accessors) {
    int n = 1000;
    int q = isqrt(n);
    sqrt_map<int, double> mm(q, n);
    const auto& cm = mm;
    for (int k = 1; k <= n; k++) {
        mm[n / k] = n / k;
    }
    for (int k = 1; k <= n; k++) {
        EXPECT_EQ(n / k, mm.at(n / k)) << "unexpected element at " << (n / k);
        EXPECT_EQ(n / k, cm.at(n / k)) << "unexpected element at " << (n / k);
        EXPECT_EQ(n / k, mm.el(n / k)) << "unexpected element el " << (n / k);
        EXPECT_EQ(n / k, cm.el(n / k)) << "unexpected element el " << (n / k);
        EXPECT_EQ(n / k, mm[n / k]) << "unexpected element [] " << (n / k);
        EXPECT_EQ(n / k, cm[n / k]) << "unexpected element [] " << (n / k);
        EXPECT_EQ(n / k, mm(n / k)) << "unexpected element [] " << (n / k);
        EXPECT_EQ(n / k, cm(n / k)) << "unexpected element [] " << (n / k);
    }
    for (int i = 1; i <= q; i++) {
        EXPECT_EQ(i, mm.lo(i)) << "unexpected element lo " << (i);
        EXPECT_EQ(i, cm.lo(i)) << "unexpected element lo " << (i);
    }
    for (int k = 1; k <= (n / (q + 1)); k++) {
        int i = n / k;
        EXPECT_EQ(i, mm.hi(k)) << "unexpected element hi " << (k);
        EXPECT_EQ(i, cm.hi(k)) << "unexpected element hi " << (k);
    }
}

TEST(sqrt_map_test, out_of_range) {
    int n = 1000;
    int q = isqrt(n);
    sqrt_map<int, double> m(q, n);
    assert_oor([&](){ m.at(-1); });
    assert_oor([&](){ m.at(n + 1); });
    assert_oor([&](){ m.at(n / 123); }, "invalid sqrt_map<I, T> key");
    m[n / 123] = 42;
    ASSERT_EQ(42, m.at(n / 123));
    assert_oor([&](){ m.at(n / 3); }, "invalid sqrt_map<I, T> key");
    m[n / 3] = 51;
    ASSERT_EQ(51, m.at(n / 3));
}

void test_insert_erase(sqrt_map<int, double>& m, int key, double val, double val2) {
    EXPECT_EQ(0, m.count(key)) << "unexpected element at " << key;
    EXPECT_EQ(make_pair(key, true), m.insert({ key, val })) << "insert failed at " << key;
    EXPECT_EQ(1, m.count(key)) << "no element at " << key;
    EXPECT_EQ(val, m.at(key)) << "unexpected element at " << key;
    EXPECT_EQ(make_pair(key, false), m.insert({ key, val2 })) << "no element at " << key;
    EXPECT_EQ(1, m.count(key)) << "no element at " << key;
    EXPECT_EQ(val, m.at(key)) << "unexpected element at " << key;
    EXPECT_EQ(make_pair(key, false), m.insert({ key, val2 })) << "no element at " << key;
    EXPECT_EQ(1, m.count(key)) << "no element at " << key;
    EXPECT_EQ(val, m.at(key)) << "unexpected element at " << key;

    EXPECT_EQ(1, m.erase(key)) << "erase failed at " << key;
    EXPECT_EQ(0, m.erase(key)) << "unexpected element at " << key;
    EXPECT_EQ(0, m.erase(key)) << "unexpected element at " << key;

    EXPECT_EQ(make_pair(key, true), m.insert({ key, val2 })) << "insert failed at " << key;
    EXPECT_EQ(1, m.count(key)) << "no element at " << key;
    EXPECT_EQ(val2, m.at(key)) << "unexpected element at " << key;
}

TEST(sqrt_map_test, insert_erase) {
    int n = 100;
    int q = isqrt(n);
    sqrt_map<int, double> m(q, n);
    test_insert_erase(m, 17, 3.14, 2.71);
    test_insert_erase(m, 5, 0.61, 1.61);
}
