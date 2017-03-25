#include "algorithm/math/primes.h"
#include "structure/math/polynom.h"

#include <algorithm>
#include <vector>
#include <map>

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

TEST(primes_test, primes_pq) {
	int n = 30;
	vector<int> vp(n);
	vector<char> vq(n);
	int m = primes(&vp[0], &vq[0], n);
	EXPECT_EQ(10, m);
	EXPECT_EQ((vector<int> { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }), vp);
	EXPECT_EQ((vector<char> { 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1 }), vq);
}

TEST(primes_test, primes_p) {
	int n = 30;
	vector<int> vp(n);
	int m = primes(&vp[0], nullptr, n); vp.resize(m);
	EXPECT_EQ(10, m);
	EXPECT_EQ((vector<int> { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 }), vp);
}

TEST(primes_test, primes_q) {
	int n = 30;
	vector<char> vq(n);
	int m = primes(nullptr, &vq[0], n);
	EXPECT_EQ(10, m);
	EXPECT_EQ((vector<char> { 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1 }), vq);
}

TEST(primes_test, prime_pi) {
	int n = 30;
	vector<int> vp(n);
	int m = primes(&vp[0], nullptr, n);

	vector<int> vpi(n);
	prime_pi(&vpi[0], n, &vp[0], m);
	EXPECT_EQ((vector<int> { 0, 0, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 10 }), vpi);
}

TEST(primes_test, euler_phi) {
	int n = 30;
	vector<int> vp(n);
	int m = primes(&vp[0], nullptr, n);

	vector<int> vphi(n);
	euler_phi(&vphi[0], n, &vp[0], m);
	EXPECT_EQ((vector<int> { 0, 1, 1, 2, 2, 4, 2, 6, 4, 6, 4, 10, 4, 12, 6, 8, 8, 16, 6, 18, 8, 12, 10, 22, 8, 20, 12, 18, 12, 28 }), vphi);
}

TEST(primes_test, moebius_mu) {
	int n = 30;
	vector<int> vmu(n);
	moebius_mu(&vmu[0], n);
	EXPECT_EQ((vector<int> { 0, +1, -1, -1, 0, -1, +1, -1, 0, 0, +1, -1, 0, -1, +1, +1, 0, -1, 0, -1, 0, +1, +1, -1, 0, 0, +1, 0, 0, -1 }), vmu);
}

TEST(primes_test, segmented_q) {
	int b = 20, e = 30;
	int q = isqrt(e) + 1;
	vector<int> vp(q);
	int m = primes(&vp[0], nullptr, q);

	vector<char> vq(e - b);
	segmented_q(&vq[0], b, e, &vp[0], m);
	EXPECT_EQ((vector<char> {0, 0, 0, 1, 0, 0, 0, 0, 0, 1}), vq);

	vq.resize(e);
	segmented_q(&vq[0], 0, e, &vp[0], m);
	EXPECT_EQ((vector<char> {0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1}), vq);
}

TEST(primes_test, segmented_phi) {
	int b = 20, e = 30;
	int q = isqrt(e) + 1;
	vector<int> vp(q);
	int m = primes(&vp[0], nullptr, q);

	vector<int64_t> vphi(e - b);
	vector<int64_t> vtmp(e - b);
	segmented_phi(&vphi[0], &vtmp[0], b, e, &vp[0], m);
	EXPECT_EQ((vector<int64_t> { 8, 12, 10, 22, 8, 20, 12, 18, 12, 28 }), vphi);

	vphi.resize(e);
	vtmp.resize(e);
	segmented_phi(&vphi[0], &vtmp[0], 0, e, &vp[0], m);
	EXPECT_EQ((vector<int64_t> {0, 1, 1, 2, 2, 4, 2, 6, 4, 6, 4, 10, 4, 12, 6, 8, 8, 16, 6, 18, 8, 12, 10, 22, 8, 20, 12, 18, 12, 28}), vphi);
}

TEST(primes_test, segmented_mu) {
	int b = 20, e = 30;
	int q = isqrt(e) + 1;
	vector<int> vp(q);
	int m = primes(&vp[0], nullptr, q);

	vector<int64_t> vmu(e - b);
	segmented_mu(&vmu[0], b, e, &vp[0], m);
	EXPECT_EQ((vector<int64_t> { 0, +1, +1, -1, 0, 0, +1, 0, 0, -1 }), vmu);

	vmu.resize(e);
	segmented_mu(&vmu[0], 0, e, &vp[0], m);
	EXPECT_EQ((vector<int64_t> {0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1, 0, -1, 0, -1, 0, 1, 1, -1, 0, 0, 1, 0, 0, -1}), vmu);
}

TEST(primes_test, divisor_sigma_0) {
	int n = 30;
	vector<int> vds0(n);
	divisor_sigma0(&vds0[0], n);
	EXPECT_EQ((vector<int> {0, 1, 2, 2, 3, 2, 4, 2, 4, 3, 4, 2, 6, 2, 4, 4, 5, 2, 6, 2, 6, 4, 4, 2, 8, 3, 4, 4, 6, 2}), vds0);
}

TEST(primes_test, divisor_sigma_1) {
	int n = 30;
	vector<int64_t> vds1(n);
	divisor_sigma1(&vds1[0], n);
	EXPECT_EQ((vector<int64_t> {0, 1, 3, 4, 7, 6, 12, 8, 15, 13, 18, 12, 28, 14, 24, 24, 31, 18, 39, 20, 42, 32, 36, 24, 60, 31, 42, 40, 56, 30}), vds1);
}

TEST(primes_test, factor) {
	int n = 30;
	vector<int> vp(n);
	int m = primes(&vp[0], nullptr, n);

	vector<int> vpf(n);
	factor(&vpf[0], n, &vp[0], m);
	EXPECT_EQ((vector<int> { 0, 1, 2, 3, 2, 5, 3, 7, 2, 3, 5, 11, 3, 13, 7, 5, 2, 17, 3, 19, 5, 7, 11, 23, 3, 5, 13, 3, 7, 29 }), vpf);
}

TEST(primes_test, factor_integer) {
    int n = 30;
    vector<int> vp(n);
    int m = primes(&vp[0], nullptr, n);
    vector<int> vpf(n);
    factor(&vpf[0], n, &vp[0], m);

    vector<int> vm0(1); factor_integer_to_map(vm0, 0, &vpf[0]);
    EXPECT_EQ((vector<int> {0}), vm0);
    vector<int> vm1(2); factor_integer_to_map(vm1, 1, &vpf[0]);
    EXPECT_EQ((vector<int> {0, 0}), vm1);
    vector<int> vm2(7); factor_integer_to_map(vm2, 6, &vpf[0]);
    EXPECT_EQ((vector<int> {0, 0, 1, 1, 0, 0, 0}), vm2);
    vector<int> vm3(18); factor_integer_to_map(vm3, 17, &vpf[0]);
    EXPECT_EQ((vector<int> {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}), vm3);
    map<int, int> vm4; factor_integer_to_map(vm4, 20, &vpf[0]);
    EXPECT_EQ((map<int, int> {{ 2, 2 }, { 5, 1 }}), vm4);

    vector<pair<int, int>> vf0; factor_integer(vf0, 0, &vpf[0]);
    EXPECT_EQ((vector<pair<int, int>> {}), vf0);
    vector<pair<int, int>> vf1; factor_integer(vf1, 1, &vpf[0]);
    EXPECT_EQ((vector<pair<int, int>> {}), vf1);
    vector<pair<int, int>> vf2; factor_integer(vf2, 2, &vpf[0]);
    EXPECT_EQ((vector<pair<int, int>> {{ 2, 1 } }), vf2);
    vector<pair<int, int>> vf17; factor_integer(vf17, 17, &vpf[0]);
    EXPECT_EQ((vector<pair<int, int>> {{ 17, 1 } }), vf17);
    vector<pair<int, int>> vf20; factor_integer(vf20, 20, &vpf[0]);
    EXPECT_EQ((vector<pair<int, int>> {{ 5, 1 }, { 2, 2 } }), vf20);

    vector<pair<int, int>> vf9800; factor_integer(vf9800, vector<int>{ 20, 14, 35 }, &vpf[0]);
    EXPECT_EQ((vector<pair<int, int>> {{ 5, 2 }, { 2, 3 }, { 7, 2 } }), vf9800);

    vector<int64_t> vd20; divisors(vd20, vf20); sort(vd20.begin(), vd20.end());
    EXPECT_EQ((vector<int64_t>{ 1, 2, 4, 5, 10, 20 }), vd20);
    vector<int64_t> vd9800; divisors(vd9800, vf9800, int64_t(49)); sort(vd9800.begin(), vd9800.end());
    EXPECT_EQ((vector<int64_t>{ 1, 2, 4, 5, 7, 8, 10, 14, 20, 25, 28, 35, 40, 49 }), vd9800);
    vector<pair<int, int>> vf1e9{ { 1000000007, 1 }, { 1000000009, 1 } };
    vector<int64_t> vd1e9; divisors(vd1e9, vf1e9); sort(vd1e9.begin(), vd1e9.end());
    EXPECT_EQ((vector<int64_t>{ 1, 1000000007, 1000000009, int64_t(1000000016000000063LL) }), vd1e9);

    EXPECT_EQ((vector<int>{5, 2}), prime_factors(vf20));
    EXPECT_EQ((vector<int>{5, 2, 7}), prime_factors(vf9800));
    EXPECT_EQ((vector<int>{1, 2}), prime_exponents(vf20));
    EXPECT_EQ((vector<int>{2, 3, 2}), prime_exponents(vf9800));

    EXPECT_EQ(1, divisor_sigma0(vf0));
    EXPECT_EQ(1, divisor_sigma0(vf1));
    EXPECT_EQ(2, divisor_sigma0(vf2));
    EXPECT_EQ(2, divisor_sigma0(vf17));
    EXPECT_EQ(6, divisor_sigma0(vf20));
    EXPECT_EQ(36, divisor_sigma0(vf9800));
}

TEST(primes_test, carmichael_lambda) {
	EXPECT_EQ(3360, euler_phi(vector<pair<int, int>> {{ 5, 2 }, { 2, 3 }, { 7, 2 } }));
	EXPECT_EQ(420, carmichael_lambda(vector<pair<int, int>> {{ 5, 2 }, { 2, 3 }, { 7, 2 } }));
	EXPECT_EQ(1, carmichael_lambda(vector<pair<int, int>> {}));
	EXPECT_EQ(1, carmichael_lambda(vector<pair<int, int>> {{ 2, 1 } }));
	EXPECT_EQ(2, carmichael_lambda(vector<pair<int, int>> {{ 2, 2 } }));
	EXPECT_EQ(2, carmichael_lambda(vector<pair<int, int>> {{ 2, 3 }}));
	EXPECT_EQ(4, carmichael_lambda(vector<pair<int, int>> {{ 2, 4 }}));
	EXPECT_EQ(256, carmichael_lambda(vector<pair<int, int>> {{ 2, 10 }}));
}

TEST(primes_test, squares_r) {
	vector<int> va, vu;
	for (int i = 1; i <= 30; i++) {
		auto vf = factor_integer_slow(i);
		va.push_back(squares_r(vf, false));
		vu.push_back(squares_r(vf, true));
	}
	EXPECT_EQ((vector<int>{4, 4, 0, 4, 8, 0, 0, 4, 4, 8, 0, 0, 8, 0, 0, 4, 8, 4, 0, 8, 0, 0, 0, 0, 12, 8, 0, 0, 8, 0}), va);
	EXPECT_EQ((vector<int>{1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 2, 1, 0, 0, 1, 0}), vu);
}

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

template<typename C> C sorted(C c) { sort(c.begin(), c.end()); return c; }

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

TEST(primes_test, integer_digits) {
	EXPECT_EQ((vector<int>{}), integer_digits(0, 10));
	EXPECT_EQ((vector<int>{1}), integer_digits(1, 10));
	EXPECT_EQ((vector<int>{3, 2, 1}), integer_digits(123, 10));
    EXPECT_EQ((vector<int>{3, 2, 1, 0, 0}), integer_digits(123, 10, 5));
	EXPECT_EQ((vector<int>{11, 7}), integer_digits(123, 16));
	EXPECT_EQ((vector<int>{1, 0, 1, 1, 0, 0, 1}), integer_digits(77, 2));
	EXPECT_EQ(77, (polynom<int>{ 1, 0, 1, 1, 0, 0, 1 }(2)));
	EXPECT_EQ(123, (polynom<int>{ 11, 7 }(16)));
	EXPECT_EQ(123, (polynom<int>{ 3, 2, 1 }(10)));
    EXPECT_EQ(123, (polynom<int>{ 3, 2, 1, 0, 0 }(10)));
}

TEST(primes_test, integer_string) {
    EXPECT_EQ("", integer_string(0, 10));
    EXPECT_EQ("0", integer_string(0, 10, 1));
    EXPECT_EQ("1", integer_string(1, 10));
    EXPECT_EQ("123", integer_string(123, 10));
    EXPECT_EQ("00123", integer_string(123, 10, 5));
    EXPECT_EQ("7b", integer_string(123, 16));
    EXPECT_EQ("1001101", integer_string(77, 2));
}
