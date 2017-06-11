#include "altruct/algorithm/math/factorization.h"
#include "altruct/algorithm/math/primes.h"
#include "altruct/structure/math/polynom.h"

#include <algorithm>
#include <vector>
#include <map>

#include "gtest/gtest.h"

#include "altruct/algorithm/collections/collections.h"

using namespace std;
using namespace altruct::math;
using namespace altruct::collections;

TEST(primes_test, miller_rabin) {
	int n = 100000;
	vector<char> vq(n);
	int m = primes(nullptr, &vq[0], n);
	vector<char> vr(n);
	for (int i = 0; i < n; i++) {
		vr[i] = miller_rabin(i);
	}
	EXPECT_EQ(vq, vr);
}

TEST(primes_test, pollard_rho) {
	EXPECT_EQ(1, pollard_rho_repeated(1));
	EXPECT_EQ(2, pollard_rho_repeated(2 * 2 * 2 * 2 * 2 * 2));
	EXPECT_EQ(27, pollard_rho_repeated(3 * 3 * 3 * 3 * 3));
	EXPECT_EQ(5, pollard_rho_repeated(5 * 5 * 5 * 7 * 7));
	EXPECT_EQ(5, pollard_rho_repeated(5 * 5 * 5 * 7 * 13 * 13));
	EXPECT_EQ(int64_t(7027), pollard_rho_repeated(int64_t(1657) * 7027));
	EXPECT_EQ(int64_t(21859), pollard_rho_repeated(int64_t(21859) * 45751));
	EXPECT_EQ(int64_t(113903), pollard_rho_repeated(int64_t(87803) * 113903));
	EXPECT_EQ(int64_t(36947), pollard_rho_repeated(int64_t(27259) * 36947));
}

TEST(primes_test, factor_integer_general_purpose) {
	typedef vector<pair<int64_t, int>> fact;
	// smooth
	EXPECT_EQ((fact{ { 2, 3 }, { 3, 5 }, { 5, 2 }, { 7, 4 }, { 13, 1 }, { 17, 2 } }), sorted(factor_integer(int64_t(438399070200LL))));
	// non-square-free
	EXPECT_EQ((fact{ { 2, 2 }, { 79, 1 }, { 24137441, 1 }, { 32046803, 1 } }), sorted(factor_integer(int64_t(244434790061754868LL))));
	// square-free
	EXPECT_EQ((fact{ { 2, 1 }, { 13, 1 }, { 11329, 1 }, { 39367, 1 }, { 11293829, 1 } }), sorted(factor_integer(int64_t(130959935583540622LL))));
	// small-big prime
	EXPECT_EQ((fact{ { 3, 1 }, { 10402882839016853, 1 } }), sorted(factor_integer(int64_t(31208648517050559LL))));
	// big semi-prime
	EXPECT_EQ((fact{ { 181153303, 1 }, { 558255521, 1 } }), sorted(factor_integer(int64_t(101129831547135863LL))));
	// big square
	EXPECT_EQ((fact{ { 549843233, 2 } }), sorted(factor_integer(int64_t(302327580875892289LL))));
	// big power
	EXPECT_EQ((fact{ { 337013, 3 } }), sorted(factor_integer(int64_t(38277182361861197LL))));
	// big prime
	EXPECT_EQ((fact{ { 988359650216386457, 1 } }), sorted(factor_integer(int64_t(988359650216386457LL))));
}

TEST(primes_test, factor_integer_general_purpose_first_1000) {
	typedef vector<pair<int, int>> fact;

	int n = 1000;
	vector<int> vp(n);
	int m = primes(&vp[0], nullptr, n);
	vector<int> vpf(n);
	factor(&vpf[0], n, &vp[0], m);

	for (int i = 1; i < n; i++) {
		fact vf; factor_integer(vf, i, &vpf[0]);
		EXPECT_EQ(sorted(vf), sorted(factor_integer(i)));
	}
}

TEST(primes_test, factor_integer_trial_division_first_1000) {
	typedef vector<pair<int, int>> fact;

	int n = 1000;
	vector<int> vp(n);
	int m = primes(&vp[0], nullptr, n);
	vector<int> vpf(n);
	factor(&vpf[0], n, &vp[0], m);

	for (int i = 1; i < n; i++) {
		fact vf; factor_integer(vf, i, &vpf[0]);
		EXPECT_EQ(sorted(vf), sorted(factor_integer_slow(i)));
	}
}

TEST(primes_test, factor_out) {
    int e1 = 1000;
    EXPECT_EQ(17, factor_out(17, 3, e1));
    EXPECT_EQ(1000, e1);
    int e2 = 1000;
    EXPECT_EQ(1, factor_out(243, 3, e2));
    EXPECT_EQ(1005, e2);
    int e3 = 1000;
    EXPECT_EQ(17, factor_out(243 * 17, 3, e3));
    EXPECT_EQ(1005, e3);
    int e4 = 1000;
    EXPECT_EQ(powT(2LL, 15), factor_out(powT(10LL, 15), 5, e4));
    EXPECT_EQ(1015, e4);
}

TEST(primes_test, from_factorization) {
    EXPECT_EQ(1, from_factorization<int>({}));
    EXPECT_EQ(7593750000000000, (from_factorization<int, int64_t>({ { 2, 10 }, { 3, 5 }, { 5, 15 } })));
}

TEST(primes_test, fraction_reduce) {
    auto gcd_f = [](int x, int y){return gcd(x, y); };
    vector<int> num0{ 2 * 6, 5, 35, 22 };
    vector<int> den0{ 5, 13, 6 * 17 };
    fraction_reduce(num0, den0, gcd_f);
    EXPECT_EQ((vector<int>{2, 1, 35, 22}), num0);
    EXPECT_EQ((vector<int>{1, 13, 17, }), den0);
    vector<int> num1{ 5, 13, 6 * 17 };
    vector<int> den1{ 2 * 6, 5, 35, 22 };
    fraction_reduce(num1, den1, gcd_f);
    EXPECT_EQ((vector<int>{1, 13, 17, }), num1);
    EXPECT_EQ((vector<int>{2, 1, 35, 22}), den1);
}
