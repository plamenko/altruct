#pragma once

#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>
#include <string>

#include "common_test_util.h"

#include "gtest/gtest.h"

namespace altruct {
namespace test_util {

namespace {
// googletest doesn't handle ASSERT_EQ(bool, bool) well...
std::string to_str(bool val) { return val ? "true" : "false"; }
}

template<typename TL, typename TR>
void assert_comparison_operators(int expected, const TL& lhs, const TR& rhs, const char* message) {
	ASSERT_EQ(to_str(expected == 0), to_str(lhs == rhs)) << message;
	ASSERT_EQ(to_str(expected != 0), to_str(lhs != rhs)) << message;
	ASSERT_EQ(to_str(expected < 0), to_str(lhs < rhs)) << message;
	ASSERT_EQ(to_str(expected >= 0), to_str(lhs >= rhs)) << message;
	ASSERT_EQ(to_str(expected <= 0), to_str(lhs <= rhs)) << message;
	ASSERT_EQ(to_str(expected > 0), to_str(lhs > rhs)) << message;
}

#define ASSERT_COMPARISON_OPERATORS(expected, lhs, rhs) \
	assert_comparison_operators(expected, lhs, rhs, ALTRUCT_AT)

} // structure
} // altruct
