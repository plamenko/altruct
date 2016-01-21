#include "structure/math/prime_holder.h"

#include "gtest/gtest.h"

#include <vector>

using namespace std;
using namespace altruct::math;

typedef long long ll;
typedef std::pair<int, int> fact_pair;

TEST(prime_holder_test, primes) {
	prime_holder prim(114);
	
	EXPECT_EQ(114, prim.size());
	EXPECT_EQ(30, prim.primes());
	
	EXPECT_EQ((vector<int>{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113}), prim.p());
	EXPECT_EQ(2, prim.p(0));
	EXPECT_EQ(113, prim.p(29));
	
	vector<char> vq{
		0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0,
		1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0,
		0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0,
		0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0,
		0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1};
	EXPECT_EQ(vq, prim.q());
	EXPECT_EQ(0, prim.q(0));
	EXPECT_EQ(0, prim.q(1));
	EXPECT_EQ(1, prim.q(2));
	EXPECT_EQ(0, prim.q(112));
	EXPECT_EQ(1, prim.q(113));
}

TEST(prime_holder, factor_integer) {
	prime_holder prim(100);

	EXPECT_EQ((vector<fact_pair> {}), prim.factor_integer(0));
	EXPECT_EQ((vector<fact_pair> {}), prim.factor_integer(1));
	EXPECT_EQ((vector<fact_pair> {{ 2, 1 } }), prim.factor_integer(2));
	EXPECT_EQ((vector<fact_pair> {{ 17, 1 } }), prim.factor_integer(17));
	EXPECT_EQ((vector<fact_pair> {{ 2, 2 }, { 5, 1 } }), prim.factor_integer(20));

	EXPECT_EQ((vector<fact_pair> {{ 2, 3 }, { 5, 2 }, { 7, 2 } }), prim.factor_integer(vector<int>{ 20, 14, 35 }));

	EXPECT_EQ((vector<ll>{ 1, 2, 4, 5, 10, 20 }), prim.divisors(20));
	EXPECT_EQ((vector<ll>{ 1, 2, 4, 5 }), prim.divisors(20, 8));
	EXPECT_EQ((vector<ll>{ 1, 2, 4, 5, 8, 10, 20, 40 }), prim.divisors(vector<int>{ 10, 4 }));
	EXPECT_EQ((vector<ll>{ 1, 2, 4, 5, 7, 8, 10, 14, 20, 25, 28, 35, 40, 49 }), prim.divisors(vector<int>{ 20, 14, 35 }, 49));
	EXPECT_EQ((vector<ll>{ 1, 1000000007, 1000000009, 1000000016000000063LL }), prim.divisors(vector<fact_pair>{{ 1000000007, 1 }, { 1000000009, 1 }}));
	EXPECT_EQ((vector<ll>{ 1, 1000000007, 1000000009 }), prim.divisors(vector<fact_pair>{{ 1000000007, 1 }, { 1000000009, 1 }}, 1000000009));
}

// TODO: pf, phi, pi, mu
