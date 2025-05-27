#include "altruct/algorithm/math/ranges.h"
#include "altruct/structure/math/modulo.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

const int P = 1000000007;
typedef moduloX<int> modx;

namespace {
vector<modx> make_vector_modx(std::initializer_list<int> l, int M) {
    vector<modx> r(l.begin(), l.end());
    for (modx& e : r) e.M() = M;
    return r;
}
} // namespace

TEST(ranges_test, range) {
    vector<modx> expected{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
    vector<modx> expected5{ 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75 };
    vector<modx> table(16); range(table.begin(), table.end(), modx(1, P));
    EXPECT_EQ(expected, table);
    EXPECT_EQ(expected, range(16, modx(1, P)));
    EXPECT_EQ(expected5, range(16, modx(5, P)));
}

TEST(ranges_test, powers) {
    vector<modx> expected{ 1, 5, 25, 125, 625, 3125, 15625, 78125, 390625, 1953125, 9765625, 48828125, 244140625, 220703118, 103515583, 517577915 };
    vector<modx> table(16); powers(table.begin(), table.end(), modx(5, P));
    EXPECT_EQ(expected, table);
    EXPECT_EQ(expected, powers(16, modx(5, P)));
}

TEST(ranges_test, factorials) {
    vector<modx> expected{ 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 227020758, 178290591, 674358851 };
    vector<modx> table(16); factorials(table.begin(), table.end(), modx(1, P));
    EXPECT_EQ(expected, table);
    EXPECT_EQ(expected, make_factorials(16, modx(1, P)));
}

TEST(ranges_test, inv_factorials) {
    vector<modx> expected{ 1, 1, 500000004, 166666668, 41666667, 808333339, 301388891, 900198419, 487524805, 831947206, 283194722, 571199524, 380933296, 490841026, 320774361, 821384963 };
    vector<modx> table(16); inv_factorials(table.begin(), table.end(), modx(1, P));
    EXPECT_EQ(expected, table);
    EXPECT_EQ(expected, make_inv_factorials(16, modx(1, P)));
    vector<modx> expected17{ 1, 1, 9, 3, 5, 1, 3, 15, 4, 8, 11, 1, 10, 6, 15, 1, 16 };
    
    EXPECT_EQ(expected17, make_inv_factorials(17, modx(1, 17)));
    EXPECT_EQ(expected17, make_inv_factorials(17, modx(1, 17), 0));
    EXPECT_EQ(expected17, make_inv_factorials(17, modx(1, 17), 1));
    EXPECT_EQ(expected17, make_inv_factorials(17, modx(2, 17), 2));
    EXPECT_EQ(expected17, make_inv_factorials(17, modx(6, 17), 3));
    EXPECT_EQ(expected17, make_inv_factorials(17, modx(7, 17), 4));
    EXPECT_EQ(expected17, make_inv_factorials(17, modx(16, 17), 16));
}

TEST(ranges_test, inverses) {
    vector<modx> expected{ 0, 1, 500000004, 333333336, 250000002, 400000003, 166666668, 142857144, 125000001, 111111112, 700000005, 818181824, 83333334, 153846155, 71428572, 466666670 };
    vector<modx> actual(16); inverses(actual.begin(), actual.end(), modx(1, P));
    EXPECT_EQ(expected, actual);
    EXPECT_EQ(expected, make_inverses(16, modx(1, P)));

    vector<modx> expected17{ 0, 1, 9, 6, 13, 7, 3, 5, 15, 2, 12, 14, 10, 4, 11, 8, 16 };
    const auto ifact17 = make_inv_factorials(17, modx(1, 17));
    vector<modx> actual17 = ifact17; inverses_from_ifact(actual17.begin(), actual17.end(), modx(1, 17));
    EXPECT_EQ(expected17, actual17);
}

TEST(ranges_test, power) {
    vector<modx> expected{ 0, 1, 1024, 59049, 1048576, 9765625, 60466176, 282475249, 73741817, 486784380, 999999937, 937424426, 917363797, 858490890, 254652953, 650386593 };
    auto table = range(16, modx(1, P));
    power(table.begin(), table.end(), 10);
    EXPECT_EQ(expected, table);
}

TEST(ranges_test, invert) {
    vector<modx> expected{ 0, 1, 500000004, 333333336, 250000002, 400000003, 166666668, 142857144, 125000001, 111111112, 700000005, 818181824, 83333334, 153846155, 71428572, 466666670 };
    auto table = range(16, modx(1, P));
    invert(table.begin(), table.end(), modx(1, P));
    EXPECT_EQ(expected, table);
}

TEST(ranges_test, invert_field) {
    vector<modx> expected = make_vector_modx({ 3, 7, 2, 8, 13, 16, 11, 4, 10, 9, 12, 6, 5, 1, 15, 14 }, 17);
    vector<modx> table = make_vector_modx({ 6, 5, 9, 15, 4, 16, 14, 13, 12, 2, 10, 3, 7, 1, 8, 11 }, 17);
    invert_field(table.begin(), table.end(), modx(1, 17));
    EXPECT_EQ(expected, table);
}

TEST(ranges_test, negate) {
    vector<modx> expected{ 0, 1000000006, 1000000005, 1000000004, 1000000003, 1000000002, 1000000001, 1000000000, 999999999, 999999998, 999999997, 999999996, 999999995, 999999994, 999999993, 999999992 };
    auto table = range(16, modx(1, P));
    altruct::math::negate(table.begin(), table.end());
    EXPECT_EQ(expected, table);
}

TEST(ranges_test, alternate) {
    vector<modx> expected{ 0, 1000000006, 2, 1000000004, 4, 1000000002, 6, 1000000000, 8, 999999998, 10, 999999996, 12, 999999994, 14, 999999992 };
    auto table = range(16, modx(1, P));
    alternate(table.begin(), table.end());
    EXPECT_EQ(expected, table);
}

TEST(ranges_test, accumulate) {
    vector<modx> expected{ 0, 1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66, 78, 91, 105, 120 };
    auto table = range(16, modx(1, P));
    accumulate(table.begin(), table.end());
    EXPECT_EQ(expected, table);
}

TEST(ranges_test, differentiate) {
    vector<modx> expected{ 1, 0, 1, 4, 18, 96, 600, 4320, 35280, 322560, 3265920, 36288000, 439084800, 748019165, 951269840, 496068260 };
    auto table = make_factorials(16, modx(1, P));
    differentiate(table.begin(), table.end());
    EXPECT_EQ(expected, table);
}
