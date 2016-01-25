#include "algorithm/math/recurrence.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

typedef long long ll;
typedef modulo<int, 1000000007> mod;

TEST(recurrence_test, linear_recurrence) {
	std::vector<int> f;
	for (int n = 0; n < 20; n++) {
		f.push_back(linear_recurrence<int>({1, -2, 3, 4, -5}, {2, 3, 5, 7, 11}, n));
	}
	EXPECT_EQ((vector<int> {2, 3, 5, 7, 11, 14, 18, 26, 41, 44, 42, 91, 173, 88, -37, 460, 1035, -509, -1787, 4361}), f);
}

TEST(recurrence_test, fibonacci) {
	std::vector<int> f;
	for (int n = 0; n < 20; n++) {
		f.push_back(fibonacci<int>(n));
	}
	EXPECT_EQ((vector<int> { 0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181}), f);
	EXPECT_EQ(mod(517691607), fibonacci<mod>(1000));
}

TEST(recurrence_test, lucas_l) {
	std::vector<int> f;
	for (int n = 0; n < 20; n++) {
		f.push_back(lucas_l<int>(n));
	}
	EXPECT_EQ((vector<int> {2, 1, 3, 4, 7, 11, 18, 29, 47, 76, 123, 199, 322, 521, 843, 1364, 2207, 3571, 5778, 9349}), f);
}

TEST(recurrence_test, lucas_u_3_2) {
	std::vector<int> f;
	for (int n = 0; n < 20; n++) {
		f.push_back(lucas_u<int>(3, 2, n));
	}
	EXPECT_EQ((vector<int> {0, 1, 3, 7, 15, 31, 63, 127, 255, 511, 1023, 2047, 4095, 8191, 16383, 32767, 65535, 131071, 262143, 524287}), f);
}

TEST(recurrence_test, lucas_v_3_2) {
	std::vector<int> f;
	for (int n = 0; n < 20; n++) {
		f.push_back(lucas_v<int>(3, 2, n));
	}
	EXPECT_EQ((vector<int> {2, 3, 5, 9, 17, 33, 65, 129, 257, 513, 1025, 2049, 4097, 8193, 16385, 32769, 65537, 131073, 262145, 524289}), f);
}

TEST(recurrence_test, lucas_u_11_10) {
	std::vector<int> f;
	for (int n = 0; n < 11; n++) {
		f.push_back(lucas_u<int>(11, 10, n));
	}
	EXPECT_EQ((vector<int> {0, 1, 11, 111, 1111, 11111, 111111, 1111111, 11111111, 111111111, 1111111111}), f);
}

TEST(recurrence_test, bernoulli_b) {
	std::vector<mod> b = bernoulli_b<mod>(10);
	EXPECT_EQ((vector<mod> {mod(1)/1, mod(1)/2, mod(1)/6, 0, -mod(1)/30, 0, mod(1)/42, 0, -mod(1)/30, 0, mod(5)/66}), b);
}
