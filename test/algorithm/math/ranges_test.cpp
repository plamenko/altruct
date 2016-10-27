#include "algorithm/math/ranges.h"
#include "structure/math/modulo.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

const int P = 1000000007;
typedef moduloX<int> modx;

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
	EXPECT_EQ(expected, factorials(16, modx(1, P)));
}

TEST(ranges_test, inv_factorials) {
	vector<modx> expected{ 1, 1, 500000004, 166666668, 41666667, 808333339, 301388891, 900198419, 487524805, 831947206, 283194722, 571199524, 380933296, 490841026, 320774361, 821384963 };
	vector<modx> table(16); inv_factorials(table.begin(), table.end(), modx(1, P));
	EXPECT_EQ(expected, table);
	EXPECT_EQ(expected, inv_factorials(16, modx(1, P)));
}

TEST(ranges_test, inverses) {
	vector<modx> expected{ 0, 1, 500000004, 333333336, 250000002, 400000003, 166666668, 142857144, 125000001, 111111112, 700000005, 818181824, 83333334, 153846155, 71428572, 466666670 };
	vector<modx> table(16); inverses(table.begin(), table.end(), modx(1, P));
	EXPECT_EQ(expected, table);
	EXPECT_EQ(expected, inverses(16, modx(1, P)));
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
	auto table = factorials(16, modx(1, P));
	differentiate(table.begin(), table.end());
	EXPECT_EQ(expected, table);
}
