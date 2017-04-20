#include "algorithm/math/recurrence.h"
#include "structure/math/matrix.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

typedef long long ll;
typedef modulo<int, 1000000007> mod;
typedef matrix<int> mat;

TEST(recurrence_test, linear_recurrence) {
	std::vector<int> f;
	for (int n = 0; n < 20; n++) {
		f.push_back(linear_recurrence<int, int>({1, -2, 3, 4, -5}, {2, 3, 5, 7, 11}, n));
	}
	EXPECT_EQ((vector<int> {2, 3, 5, 7, 11, 14, 18, 26, 41, 44, 42, 91, 173, 88, -37, 460, 1035, -509, -1787, 4361}), f);
}

TEST(recurrence_test, linear_recurrence_next) {
	std::vector<int> f{ 2, 3, 5, 7, 11 };
	while (f.size() < 20) {
		f.push_back(linear_recurrence_next<int>({ 1, -2, 3, 4, -5 }, f));
	}
	EXPECT_EQ((vector<int> {2, 3, 5, 7, 11, 14, 18, 26, 41, 44, 42, 91, 173, 88, -37, 460, 1035, -509, -1787, 4361}), f);
}

TEST(recurrence_test, linear_recurrence_on_matrix) {
	std::vector<mat> a;
	mat a0 = { { 1, 0 }, { 0, 1 } };
	mat a1 = { { 1, 1 }, { 1, 0 } };
	mat a2 = { { 3, 1 }, { 1, 2 } };
	mat a3 = { { 5, 3 }, { 3, 2 } };
	mat a4 = { { 11, 5 }, { 5, 6 } };
	for (int n = 0; n < 5; n++) {
		a.push_back(linear_recurrence<int, mat>({ 1, 2 }, { a0, a1 }, n));
	}
	EXPECT_EQ((vector<mat> {a0, a1, a2, a3, a4}), a);
}

TEST(recurrence_test, linear_recurrence_on_matrix_matrix) {
    std::vector<mat> a;
    mat c0 = { { 2, 0 }, { 0, 4 } };
    mat c1 = { { -3, 5 }, { 7, 0 } };
    mat a0 = { { 1, 0 }, { 0, 1 } };
    mat a1 = { { 1, 1 }, { 1, 0 } };
    mat a2 = { { -1, 7 }, { 11, 0 } };
    mat a3 = { { 0, 11 }, { 51, 7 } };
    mat a4 = { { 58, 1 }, { 197, 77 } };
    for (int n = 0; n < 5; n++) {
        a.push_back(linear_recurrence<mat, mat>({ c0, c1 }, { a0, a1 }, n));
    }
    EXPECT_EQ((vector<mat> {a0, a1, a2, a3, a4}), a);
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

TEST(recurrence_test, berlekamp_massey) {
	// a[n+1] = 17 a[n-0] - 23 a[n-1] + 13 a[n-2] + 45 a[n-3] - 58 a[n-4]
	// x^(n+1) = 17 x^(n-0) - 23 x^(n-1) + 13 x^(n-2) + 45 x^(n-3) - 58 x^(n-4)   / x^(n-4)
	// x^5 = 17 x^4 - 23 x^3 + 13 x^2 + 45 x^1 - 58 x^0
	// 58 x^0 - 45 x^1 - 13 x^2 + 23 x^3 - 17 x^4 + 1 x^5 = 0
	std::vector<mod> a;
	for (int n = 0; n <= 100; n++) {
		a.push_back(linear_recurrence<mod, mod>({ 17, -23, 13, 45, -58 }, { 2, 3, 5, 7, 11 }, n));
	}
	auto p = berlekamp_massey<mod>(a);
	EXPECT_EQ((polynom<mod> { +58, -45, -13, +23, -17, 1 }), p);
	
	// use the characteristic polynomial to calculate the n-th term of the sequence
	typedef moduloX<polynom<mod>> polymod;
	auto x_n = powT(polymod({ 0, 1 }, p), 100); // x^n % p
	mod r = 0;
	for (int i = 0; i <= x_n.v.deg(); i++) {
		r += x_n.v[i] * a[i];
	}
	EXPECT_EQ(a[100], r);
}
