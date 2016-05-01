#include "algorithm/hash/std_hash_combine.h"
#include "algorithm/hash/std_tuple_hash.h"

#include "gtest/gtest.h"

#include <string>
#include <unordered_map>

using namespace std;

TEST(std_hash_test, hash_combine) {
	size_t seed = 0x1234;
	hash_combine(seed, 0x5678);
	EXPECT_EQ(0xccf8a653c6a270baULL, seed);
	hash_combine(seed, 0xaaaaaaaaaaaaaaaa);
	EXPECT_EQ(0x8e390c1b98f4ba31ULL, seed);
}

TEST(std_hash_test, tuple_hash) {
	unordered_map<tuple<double, string>, int> m;
	m[make_tuple(2.5, "abc")] = 5;
	m[make_tuple(-4, "xxx")] = 11;
	EXPECT_EQ(5, m[make_tuple(2.5, "abc")]);
	EXPECT_EQ(11, m[make_tuple(-4, "xxx")]);
}
