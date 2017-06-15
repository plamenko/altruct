#include "altruct/algorithm/random/random.h"

#include <map>

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::random;

TEST(random_test, integer_to_double_0_1) {
    double eps = 1e-14;
    // 0 inclusive, 1 inclusive, 32bit
    EXPECT_TRUE(integer_to_double_0_1<uint32_t>(0x00000000U) == 0.0);
    EXPECT_TRUE(integer_to_double_0_1<uint32_t>(0xFFFFFFFFU) == 1.0);
    // 0 inclusive, 1 inclusive, 64bit
    EXPECT_TRUE(integer_to_double_0_1<uint64_t>(0x0000000000000000ULL) == 0.0);
    EXPECT_TRUE(integer_to_double_0_1<uint64_t>(0xFFFFFFFFFFFFFFFFULL) == 1.0);
    // in between values
    EXPECT_NEAR(0.0, integer_to_double_0_1<uint64_t>(0x0000000000000000ULL), eps);
    EXPECT_NEAR(1.0, integer_to_double_0_1<uint64_t>(0xFFFFFFFFFFFFFFFFULL), eps);
    EXPECT_NEAR(0.5, integer_to_double_0_1<uint64_t>(0x7FFFFFFFFFFFFFFFULL), eps);
    EXPECT_NEAR(0.5, integer_to_double_0_1<uint64_t>(0x8000000000000000ULL), eps);
    EXPECT_NEAR(1.0 / 5, integer_to_double_0_1<uint64_t>(0x3333333333333333ULL), eps);
    EXPECT_NEAR(1.0 / 3, integer_to_double_0_1<uint64_t>(0x5555555555555555ULL), eps);
}

TEST(random_test, integer_to_range) {
    EXPECT_EQ(100ULL, integer_to_range<uint64_t>(0, 100, 1100 - 1));
    EXPECT_EQ(123ULL, integer_to_range<uint64_t>(23, 100, 1100 - 1));
    EXPECT_EQ(100ULL, integer_to_range<uint64_t>(1000, 100, 1100 - 1));
    EXPECT_EQ(100ULL, integer_to_range<uint64_t>(5000, 100, 1100 - 1));
    EXPECT_EQ(123ULL, integer_to_range<uint64_t>(5023, 100, 1100 - 1));
}

TEST(random_test, biggest_multiple) {
    // 32bit
    EXPECT_EQ(4200000000U, biggest_multiple<uint32_t>(100000000U));
    EXPECT_EQ(2200000000U, biggest_multiple<uint32_t>(2200000000U));
    // 64bit
    EXPECT_EQ(18400000000000000000ULL, biggest_multiple<uint64_t>(100000000000000000ULL));
    EXPECT_EQ(9300000000000000000ULL, biggest_multiple<uint64_t>(9300000000000000000ULL));
    // powers of two
    EXPECT_EQ(0ULL, biggest_multiple<uint64_t>(0ULL));
    EXPECT_EQ(0ULL, biggest_multiple<uint64_t>(1ULL));
    EXPECT_EQ(0ULL, biggest_multiple<uint64_t>(2ULL));
    EXPECT_EQ(0ULL, biggest_multiple<uint64_t>(0x4000000000000000ULL));
    EXPECT_EQ(0ULL, biggest_multiple<uint64_t>(0x8000000000000000ULL));
    // small numbers
    EXPECT_EQ(0xfffffffffffffff0ULL, biggest_multiple<uint64_t>(30ULL));
    EXPECT_EQ(0xfffffffffffffffaULL, biggest_multiple<uint64_t>(10ULL));
    EXPECT_EQ(0xfffffffffffffff0ULL, biggest_multiple<uint64_t>(100ULL));
    EXPECT_EQ(0xfffffffffffffd98ULL, biggest_multiple<uint64_t>(1000ULL));
    EXPECT_EQ(0xfffffffffff79540ULL, biggest_multiple<uint64_t>(1000000ULL));
    // big numbers
    EXPECT_EQ(0xFFFFFFFFFFFFFFFFULL, biggest_multiple<uint64_t>(0xFFFFFFFFFFFFFFFFULL));
    EXPECT_EQ(0x8000000000000001ULL, biggest_multiple<uint64_t>(0x8000000000000001ULL));
    EXPECT_EQ(0xFFFFFFFFFFFFFFFEULL, biggest_multiple<uint64_t>(0x7FFFFFFFFFFFFFFFULL));
    EXPECT_EQ(0xAAAAAAAAAAAAAAACULL, biggest_multiple<uint64_t>(0x5555555555555556ULL));
    EXPECT_EQ(0xFFFFFFFFFFFFFFFFULL, biggest_multiple<uint64_t>(0x5555555555555555ULL));
}

TEST(random_test, uniform_next) {
    vector<int> values(256);
    for (int i = 0; i < values.size(); i++) {
        values[i] = i;
    }
    random_shuffle(values.begin(), values.end());

    auto next = [&](){
        static int index = 0;
        index %= values.size();
        return values[index++];
    };

    map<uint8_t, int> hist;
    for (int i = 0; i < 100 * 10; i++) {
        //auto r = integer_to_range<uint8_t>(next(), 10, 10 + 100 - 1); // fails
        auto r = uniform_next<uint8_t>(next, 10, 10 + 100 - 1); // passes
        hist[r]++;
    }
    for (const auto& entry : hist) {
        EXPECT_EQ(10, entry.second);
    }
}
