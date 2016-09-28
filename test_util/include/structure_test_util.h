#pragma once

#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>

#include "common_test_util.h"

#include "gtest/gtest.h"

namespace altruct {
namespace test_util {

template<typename TL, typename TR>
void assert_comparison_operators(int expected, const TL& lhs, const TR& rhs, const char* message) {
	ASSERT_EQ(expected == 0, lhs == rhs) << message;
	ASSERT_EQ(expected != 0, lhs != rhs) << message;
	ASSERT_EQ(expected < 0, lhs < rhs) << message;
	ASSERT_EQ(expected >= 0, lhs >= rhs) << message;
	ASSERT_EQ(expected <= 0, lhs <= rhs) << message;
	ASSERT_EQ(expected > 0, lhs > rhs) << message;
}

#define ASSERT_COMPARISON_OPERATORS(expected, lhs, rhs) \
	assert_comparison_operators(expected, lhs, rhs, ALTRUCT_AT)

} // structure
} // altruct
