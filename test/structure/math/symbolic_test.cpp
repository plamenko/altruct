#include "altruct/structure/math/symbolic.h"

#include "structure_test_util.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;
using namespace altruct::test_util;

TEST(symbolic_test, constructor) {
	symbolic s("s");
	EXPECT_EQ("s", s.v);
}

TEST(symbolic_test, operators_arithmetic) {
    symbolic x("x");
    symbolic y("y");
    EXPECT_EQ("(+x)", (+x).v);
    EXPECT_EQ("(-x)", (-x).v);
    EXPECT_EQ("(x+y)", (x + y).v);
    EXPECT_EQ("(x-y)", (x - y).v);
    EXPECT_EQ("(x*y)", (x * y).v);
    EXPECT_EQ("(x/y)", (x / y).v);
    EXPECT_EQ("(x%y)", (x % y).v);
    EXPECT_EQ("(~x)", (~x).v);
    EXPECT_EQ("(x&y)", (x & y).v);
    EXPECT_EQ("(x|y)", (x | y).v);
    EXPECT_EQ("(x^y)", (x ^ y).v);
    EXPECT_EQ("(x<<y)", (x << y).v);
    EXPECT_EQ("(x>>y)", (x >> y).v);
    EXPECT_EQ("(!x)", (!x).v);
    EXPECT_EQ("(x&&y)", (x && y).v);
    EXPECT_EQ("(x||y)", (x || y).v);
}

TEST(symbolic_test, operators_inplace) {
    symbolic x("x");
    symbolic y("y");
    symbolic r("r");
    r = x; r += y;
    EXPECT_EQ("(x+y)", r.v);
    r = x; r -= y;
    EXPECT_EQ("(x-y)", r.v);
    r = x; r *= y;
    EXPECT_EQ("(x*y)", r.v);
    r = x; r /= y;
    EXPECT_EQ("(x/y)", r.v);
    r = x; r %= y;
    EXPECT_EQ("(x%y)", r.v);
    r = x; r &= y;
    EXPECT_EQ("(x&y)", r.v);
    r = x; r |= y;
    EXPECT_EQ("(x|y)", r.v);
    r = x; r ^= y;
    EXPECT_EQ("(x^y)", r.v);
    r = x; r <<= y;
    EXPECT_EQ("(x<<y)", r.v);
    r = x; r >>= y;
    EXPECT_EQ("(x>>y)", r.v);
}

TEST(symbolic_test, operators_inplace_self) {
    symbolic x("x");
    symbolic r("r");
    r = x; r += r;
    EXPECT_EQ("(x+x)", r.v);
    r = x; r -= r;
    EXPECT_EQ("(x-x)", r.v);
    r = x; r *= r;
    EXPECT_EQ("(x*x)", r.v);
    r = x; r /= r;
    EXPECT_EQ("(x/x)", r.v);
    r = x; r %= r;
    EXPECT_EQ("(x%x)", r.v);
    r = x; r &= r;
    EXPECT_EQ("(x&x)", r.v);
    r = x; r |= r;
    EXPECT_EQ("(x|x)", r.v);
    r = x; r ^= r;
    EXPECT_EQ("(x^x)", r.v);
    r = x; r <<= r;
    EXPECT_EQ("(x<<x)", r.v);
    r = x; r >>= r;
    EXPECT_EQ("(x>>x)", r.v);
}

TEST(symbolic_test, casts) {
    symbolic s("s");
    symbolic e0 = zeroOf(s);
    symbolic e1 = identityOf(s);
    EXPECT_EQ("0", e0.v);
    EXPECT_EQ("1", e1.v);
    symbolic s5 = castOf<symbolic>(5);
    EXPECT_EQ("5", s5.v);
}
