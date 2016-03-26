#include "algorithm/math/modulos.h"
#include "algorithm/collections/collections.h"

#include "gtest/gtest.h"

#include <vector>

using namespace std;
using namespace altruct::math;

void crt_test_impl(int a1, int n1, int a2, int n2) {
	int a, n; chinese_remainder(a, n, a1, n1, a2, n2);
	EXPECT_EQ(lcm(n1, n2), n);
	EXPECT_EQ(a1, a % n1);
	EXPECT_EQ(a2, a % n2);
	EXPECT_GE(a, 0);
	EXPECT_LT(a, n);
}

TEST(modulos_test, chinese_remainder) {
	crt_test_impl(0, 10, 5, 13);
	crt_test_impl(5, 10, 3, 13);
	crt_test_impl(5, 10, 3, 14);
	crt_test_impl(4, 10, 6, 14);
	crt_test_impl(6, 14, 6, 14);
	crt_test_impl(102, 65535, 12345, 48888);
}

TEST(modulos_test, garner) {
	typedef moduloX<int> modx;
	typedef std::vector<modx> modv;
	
	modv a = { { 1000, 1009 }, { 1000, 1013 }, { 1000, 1019 } };
	modv a3 = altruct::collections::transform(a, [](const modx& e){ return powT(e, 3); });

	moduloX<int> r0{ 0, 1 };
	for (int i = 0; i < a3.size(); i++) {
		chinese_remainder<int>(r0.v, r0.M, a3[i].v, a3[i].M);
	}
	EXPECT_EQ(1000000000, r0.v);
	EXPECT_EQ(1009*1013*1019, r0.M);

	modv x3 = garner(a3);
	moduloX<int> r1{ 0, 1 };
	for (int i = 0; i < x3.size(); i++) {
		r1.v += r1.M * x3[i].v;
		r1.M *= x3[i].M;
	}
	EXPECT_EQ(1000000000, r1.v);
	EXPECT_EQ(1009 * 1013 * 1019, r1.M);
}

TEST(modulos_test, jacobi) {
	vector<vector<int>> v(1);
	for (int m = 1; m <= 50; m++) {
		v.push_back(vector<int>());
		for (int n = 0; n <= 20; n++) {
			v[m].push_back(jacobi(n, m));
		}
	}
	EXPECT_EQ((vector<int>{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}), v[1]);
	EXPECT_EQ((vector<int>{0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1}), v[3]);
	EXPECT_EQ((vector<int>{0, 1, -1, -1, 1, 0, 1, -1, -1, 1, 0, 1, -1, -1, 1, 0, 1, -1, -1, 1, 0}), v[5]);
	EXPECT_EQ((vector<int>{0, 1, 1, -1, 1, -1, -1, 0, 1, 1, -1, 1, -1, -1, 0, 1, 1, -1, 1, -1, -1}), v[7]);
	EXPECT_EQ((vector<int>{0, 1, -1, 0, 1, 0, 0, -1, -1, 0, 0, 1, 0, -1, 1, 0, 1, -1, 0, 1, 0}), v[45]);
}

TEST(modulos_test, sqrt_cipolla) {
	EXPECT_EQ(0, sqrt_cipolla(0, 17));
	EXPECT_EQ(1, sqrt_cipolla(1, 17));
	EXPECT_EQ(6, sqrt_cipolla(2, 17));
	EXPECT_EQ(15, sqrt_cipolla(4, 17));
	EXPECT_EQ(12, sqrt_cipolla(8, 17));
	EXPECT_EQ(14, sqrt_cipolla(9, 17));
	EXPECT_EQ(8, sqrt_cipolla(13, 17));
	EXPECT_EQ(7, sqrt_cipolla(15, 17));
	EXPECT_EQ(13, sqrt_cipolla(16, 17));
	typedef modulo<int, 17> mod;
	EXPECT_EQ(mod(0), sqT(sqrt_cipolla(mod(0))));
	EXPECT_EQ(mod(1), sqT(sqrt_cipolla(mod(1))));
	EXPECT_EQ(mod(2), sqT(sqrt_cipolla(mod(2))));
	EXPECT_EQ(mod(4), sqT(sqrt_cipolla(mod(4))));
	EXPECT_EQ(mod(8), sqT(sqrt_cipolla(mod(8))));
	EXPECT_EQ(mod(9), sqT(sqrt_cipolla(mod(9))));
	EXPECT_EQ(mod(13), sqT(sqrt_cipolla(mod(13))));
	EXPECT_EQ(mod(15), sqT(sqrt_cipolla(mod(15))));
	EXPECT_EQ(mod(16), sqT(sqrt_cipolla(mod(16))));
}

TEST(modulos_test, sqrt_hensel_lift) {
	EXPECT_EQ(0, sqrt_hensel_lift(0, 17, 5));
	EXPECT_EQ(1, sqrt_hensel_lift(1, 17, 5));
	EXPECT_EQ(461199, sqrt_hensel_lift(2, 17, 5));
	EXPECT_EQ(1419855, sqrt_hensel_lift(4, 17, 5));
	EXPECT_EQ(922398, sqrt_hensel_lift(8, 17, 5));
	EXPECT_EQ(1419854, sqrt_hensel_lift(9, 17, 5));
	EXPECT_EQ(499740, sqrt_hensel_lift(13, 17, 5));
	EXPECT_EQ(1318629, sqrt_hensel_lift(15, 17, 5));
	EXPECT_EQ(1419853, sqrt_hensel_lift(16, 17, 5));
	EXPECT_EQ(883131, sqrt_hensel_lift(12346, 17, 5));
}

TEST(modulos_test, primitive_root) {
	prime_holder prim(100);
	EXPECT_EQ(1, primitive_root(2, 1, vector<int>{ }));
	EXPECT_EQ(2, primitive_root(3, 2, vector<int>{ 2 }));
	EXPECT_EQ(3, primitive_root(4, 2, vector<int>{ 2 }));
	EXPECT_EQ(2, primitive_root(5, 4, vector<int>{ 2 }));
	EXPECT_EQ(5, primitive_root(6, 2, vector<int>{ 2 }));
	EXPECT_EQ(3, primitive_root(7, 6, vector<int>{ 2, 3 }));
	EXPECT_EQ(0, primitive_root(8, 4, vector<int>{ 2 }));
	EXPECT_EQ(2, primitive_root(9, 6, vector<int>{ 2, 3 }));
	EXPECT_EQ(3, primitive_root(10, 4, vector<int>{ 2 }));
	EXPECT_EQ(2, primitive_root(11, 10, vector<int>{ 2, 5 }));
	EXPECT_EQ(5, primitive_root(18, 6, vector<int>{ 2, 3 }));
	vector<int> vg;
	for (int m = 2; m <= 20; m++) {
		vg.push_back(primitive_root(m, prim));
	}
	EXPECT_EQ((vector<int>{1, 2, 3, 2, 5, 3, 0, 2, 3, 2, 0, 2, 3, 0, 0, 3, 5, 2, 0}), vg);
}

TEST(modulos_test, kth_roots) {
	prime_holder prim(100);
	EXPECT_EQ((vector<int>{ 1, 4, 13, 16 }), kth_roots(17, 4, 16, 16, vector<int>{2}));
	EXPECT_EQ((vector<int>{1, 17}), kth_roots(18, 4, 6, 6, vector<int>{2, 3}));
	EXPECT_EQ((vector<int>{1, 7, 13}), kth_roots(18, 3, 6, 6, vector<int>{2, 3}));
	EXPECT_EQ((vector<int>{1, 4, 13, 16}), kth_roots(17, 4, prim));
	EXPECT_EQ((vector<int>{1, 7, 13}), kth_roots(18, 3, prim));
	EXPECT_EQ((vector<int>{1, 17}), kth_roots(18, 4, prim));
	EXPECT_EQ((vector<int>{1}), kth_roots(18, 5, prim));
	EXPECT_EQ((vector<int>{1, 5, 7, 11, 13, 17}), kth_roots(18, 6, prim));
}

TEST(modulos_test, factorial_mod_p) {
	int p = 37;
	vector<int> table(p);
	factorial_mod_p_table(p, &table[0]);
	EXPECT_EQ((vector<int>{1, 1, 2, 6, 24, 9, 17, 8, 27, 21, 25, 16, 7, 17, 16, 18, 29, 12, 31, 34, 14, 35, 30, 24, 21, 7, 34, 30, 26, 14, 13, 33, 20, 31, 18, 1, 36}), table);
	
	long long e = 0;
	EXPECT_EQ(25, factorial_mod_p(10LL, p, &table[0], &e));
	EXPECT_EQ(0, e);
	EXPECT_EQ(31, factorial_mod_p(100LL, p, &table[0], &e));
	EXPECT_EQ(2, e);
	EXPECT_EQ(7, factorial_mod_p(1000LL, p, &table[0], &e));
	EXPECT_EQ(27, e);
	EXPECT_EQ(19, factorial_mod_p(10000LL, p, &table[0], &e));
	EXPECT_EQ(277, e);
	EXPECT_EQ(3, factorial_mod_p(100000LL, p, &table[0], &e));
	EXPECT_EQ(2776, e);
	EXPECT_EQ(30, factorial_mod_p(1000000LL, p, &table[0], &e));
	EXPECT_EQ(27776, e);
}
