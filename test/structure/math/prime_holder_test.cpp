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

TEST(prime_holder_test, other) {
	prime_holder prim(30);
	EXPECT_EQ((vector<int>{0, 1, 2, 3, 2, 5, 3, 7, 2, 3, 5, 11, 3, 13, 7, 5, 2, 17, 3, 19, 5, 7, 11, 23, 3, 5, 13, 3, 7, 29}), prim.pf());
	EXPECT_EQ(0, prim.pf(0));
	EXPECT_EQ(7, prim.pf(28));
	EXPECT_EQ((vector<int>{0, 1, 1, 2, 2, 4, 2, 6, 4, 6, 4, 10, 4, 12, 6, 8, 8, 16, 6, 18, 8, 12, 10, 22, 8, 20, 12, 18, 12, 28}), prim.phi());
	EXPECT_EQ(0, prim.phi(0));
	EXPECT_EQ(12, prim.phi(28));
	EXPECT_EQ((vector<int>{0, 0, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 10}), prim.pi());
	EXPECT_EQ(0, prim.pi(0));
	EXPECT_EQ(10, prim.pi(29));
	EXPECT_EQ((vector<int>{0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1, 0, -1, 0, -1, 0, 1, 1, -1, 0, 0, 1, 0, 0, -1}), prim.mu());
	EXPECT_EQ(0, prim.mu(0));
	EXPECT_EQ(-1, prim.mu(29));
}
