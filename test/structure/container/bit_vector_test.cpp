#include "algorithm/math/base.h"
#include "structure/container/bit_vector.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;
using namespace altruct::container;

TEST(bit_vector_test, constructor) {
}
TEST(bit_vector_test, bit_proxy) {
	bit_vector<> bv(10);

	bv[4] = 0; // 0 => 0
	EXPECT_EQ(0, (int)bv[4]);
	bv[4] = 1; // 0 => 1
	EXPECT_EQ(1, (int)bv[4]);
	bv[4] = 1; // 1 => 1
	EXPECT_EQ(1, (int)bv[4]);
	bv[4] = 0; // 1 => 0
	EXPECT_EQ(0, (int)bv[4]);
	
	bv[4] ^= 0; // 0 ^ 0 => 0
	EXPECT_EQ(0, (int)bv[4]);
	bv[4] ^= 1; // 0 ^ 1 => 1
	EXPECT_EQ(1, (int)bv[4]);
	bv[4] ^= 0; // 1 ^ 0 => 1
	EXPECT_EQ(1, (int)bv[4]);
	bv[4] ^= 1; // 1 ^ 1 => 0
	EXPECT_EQ(0, (int)bv[4]);
	
	bv[4] |= 0; // 0 | 0 => 0
	EXPECT_EQ(0, (int)bv[4]);
	bv[4] |= 1; // 0 | 1 => 1
	EXPECT_EQ(1, (int)bv[4]);
	bv[4] |= 1; // 1 | 1 => 1
	EXPECT_EQ(1, (int)bv[4]);
	bv[4] |= 0; // 1 | 0 => 1
	EXPECT_EQ(1, (int)bv[4]);

	bv[4] &= 1; // 1 & 1 => 1
	EXPECT_EQ(1, (int)bv[4]);
	bv[4] &= 0; // 1 & 0 => 0
	EXPECT_EQ(0, (int)bv[4]);
	bv[4] &= 1; // 0 & 1 => 0
	EXPECT_EQ(0, (int)bv[4]);
	bv[4] &= 0; // 0 & 0 => 0
	EXPECT_EQ(0, (int)bv[4]);
}
