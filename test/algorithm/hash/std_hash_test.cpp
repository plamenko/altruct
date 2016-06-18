#include "algorithm/hash/std_hash_combine.h"
#include "algorithm/hash/std_tuple_hash.h"

#include "gtest/gtest.h"

#include <string>
#include <unordered_map>

using namespace std;

TEST(std_hash_test, hash_combine) {
	// std::hash is implementation-dependent
	// so the best we can do is to check that
	// the hash changes
	size_t seed = 0x1234;
	size_t seed0 = seed;
	hash_combine(seed, 0x5678);
	EXPECT_NE(seed0, seed);
	EXPECT_NE(0, seed);
	EXPECT_NE(0x5678, seed);
	size_t seed1 = seed;
	hash_combine(seed, 0xaaaaaaaaaaaaaaaa);
	EXPECT_NE(seed1, seed);
	EXPECT_NE(0, seed);
	EXPECT_NE(0xaaaaaaaaaaaaaaaa, seed);
}

TEST(std_hash_test, tuple_hash) {
	unordered_map<tuple<double, string>, int> m;
	m[make_tuple(2.5, "abc")] = 5;
	m[make_tuple(-4, "xxx")] = 11;
	EXPECT_EQ(5, m[make_tuple(2.5, "abc")]);
	EXPECT_EQ(11, m[make_tuple(-4, "xxx")]);
}

TEST(std_hash_test, pair_hash) {
	unordered_map<pair<double, string>, int> m;
	m[make_pair(2.5, "abc")] = 5;
	m[make_pair(-4, "xxx")] = 11;
	EXPECT_EQ(5, m[make_pair(2.5, "abc")]);
	EXPECT_EQ(11, m[make_pair(-4, "xxx")]);
}
