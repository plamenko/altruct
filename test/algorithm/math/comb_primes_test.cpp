#include "algorithm/math/comb_primes.h"

#include <algorithm>
#include <vector>
#include <string>

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

typedef long long ll;

TEST(comb_primes_test, factorial_prime_exponent) {
	EXPECT_EQ(0, factorial_prime_exponent(2, 1));
	
	EXPECT_EQ(1, factorial_prime_exponent(2, 2));
	EXPECT_EQ(0, factorial_prime_exponent(3, 2));
	
	EXPECT_EQ(1, factorial_prime_exponent(2, 3));
	EXPECT_EQ(1, factorial_prime_exponent(3, 3));
	EXPECT_EQ(0, factorial_prime_exponent(5, 3));

	EXPECT_EQ(8, factorial_prime_exponent(2, 10));
	EXPECT_EQ(4, factorial_prime_exponent(3, 10));
	EXPECT_EQ(2, factorial_prime_exponent(5, 10));
	EXPECT_EQ(1, factorial_prime_exponent(7, 10));
	EXPECT_EQ(0, factorial_prime_exponent(11, 10));
}

TEST(comb_primes_test, binomial_prime_exponent) {
	vector<int> vp = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101 };
	vector<int> e = { 1, 4, 1, 2, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0 };
	vector<int> a; for (int p : vp) a.push_back(binomial_prime_exponent(p, 100, 20));
	EXPECT_EQ(e, a);
}

TEST(comb_primes_test, multinomial_prime_exponent) {
	vector<int> vp = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47 };
	vector<int> e = { 4, 5, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 0, 0, 0 };
	vector<int> a; for (int p : vp) a.push_back(multinomial_prime_exponent(p, vector<int>{ 8, 17, 11, 2 }));
	EXPECT_EQ(e, a);
}

TEST(comb_primes_test, elements_multinomial) {
	EXPECT_EQ(1.0, elements_multinomial(vector<string>{}, 1.0));
	EXPECT_EQ(1.0, elements_multinomial(vector<string>{ "aa" }, 1.0));
	EXPECT_EQ(1.0, elements_multinomial(vector<string>{ "aa", "aa", "aa", "aa" }, 1.0));
	EXPECT_EQ(24.0, elements_multinomial(vector<string>{ "aa", "bbb", "c", "dddd" }, 1.0));
	EXPECT_EQ(60.0, elements_multinomial(vector<string>{ "aa", "aa", "b", "b", "b", "ccc" }, 1.0));
	EXPECT_EQ(60, elements_multinomial(vector<string>{ "aa", "aa", "b", "b", "b", "ccc" }, 1));
}
