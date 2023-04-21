#include "altruct/algorithm/math/intrinsic.h"
#include <limits>

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

template<typename I>
void test_impl(const char* type_str) {
    I r = 0;
    I M = std::numeric_limits<I>::max();
    
    EXPECT_EQ(false, add_overflow(I(30), I(0), &r)) << type_str;
    EXPECT_EQ(30, r) << type_str;
    EXPECT_EQ(false, add_overflow(I(0), I(40), &r)) << type_str;
    EXPECT_EQ(40, r) << type_str;
    EXPECT_EQ(false, add_overflow(I(30), I(40), &r)) << type_str;
    EXPECT_EQ(70, r) << type_str;
    EXPECT_EQ(false, add_overflow(M, I(0), &r)) << type_str;
    EXPECT_EQ(M, r) << type_str;
    EXPECT_EQ(false, add_overflow(I(0), I(0), &r)) << type_str;
    EXPECT_EQ(0, r) << type_str;
    EXPECT_EQ(false, add_overflow(I(0), M, &r)) << type_str;
    EXPECT_EQ(M, r) << type_str;

    EXPECT_EQ(true, add_overflow(M, I(10), &r)) << type_str;
    EXPECT_EQ(9, r) << type_str;
    EXPECT_EQ(true, add_overflow(M, M - I(20), &r)) << type_str;
    EXPECT_EQ(I(-22), r) << type_str;
    EXPECT_EQ(true, add_overflow(I(10), M, &r)) << type_str;
    EXPECT_EQ(9, r) << type_str;
    EXPECT_EQ(true, add_overflow(M - I(20), M, &r)) << type_str;
    EXPECT_EQ(I(-22), r) << type_str;
    EXPECT_EQ(true, add_overflow(M - I(30), M - I(40), &r)) << type_str;
    EXPECT_EQ(I(-72), r) << type_str;
}

TEST(intrinsic_test, add_overflow_uint64) {
    test_impl<uint8_t>("uint8_t");
    test_impl<uint16_t>("uint16_t");
    test_impl<uint32_t>("uint32_t");
    test_impl<uint64_t>("uint64_t");
}
