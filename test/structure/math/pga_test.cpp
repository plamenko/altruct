#include "altruct/structure/math/pga.h"
#include "altruct/structure/math/symbolic.h"
#include "altruct/io/iostream_overloads.h"

#include "structure_test_util.h"

#include "gtest/gtest.h"

#include <sstream>

using namespace std;
using namespace altruct::math;
using namespace altruct::test_util;

namespace {
template<typename T>
std::string to_string(const T& v) { std::stringstream ss; ss << v; return ss.str(); }
} // namespace

TEST(pga_test, constructor_blade1) {
	auto d1 = pga::blade1<symbolic>();
	EXPECT_EQ("?", d1.e0.v);
	EXPECT_EQ("?", d1.v.x.v);
	EXPECT_EQ("?", d1.v.y.v);
	EXPECT_EQ("?", d1.v.z.v);
	auto s1 = pga::blade1<symbolic>({ "se0" });
	EXPECT_EQ("se0", s1.e0.v);
	EXPECT_EQ("0", s1.v.x.v);
	EXPECT_EQ("0", s1.v.y.v);
	EXPECT_EQ("0", s1.v.z.v);
	auto a1 = pga::blade1<symbolic>({ "ae0" }, { {"avx"}, {"avy"}, {"avz"} });
	EXPECT_EQ("ae0", a1.e0.v);
	EXPECT_EQ("avx", a1.v.x.v);
	EXPECT_EQ("avy", a1.v.y.v);
	EXPECT_EQ("avz", a1.v.z.v);
}

TEST(pga_test, operators_arithmetic_blade1) {
	auto a1 = pga::blade1<symbolic>({ "ae0" }, { {"avx"}, {"avy"}, {"avz"} });
	auto b1 = pga::blade1<symbolic>({ "be0" }, { {"bvx"}, {"bvy"}, {"bvz"} });
	EXPECT_EQ("(+ae0) e0 + (+avx) e1 + (+avy) e2 + (+avz) e3", to_string(+a1));
	EXPECT_EQ("(-ae0) e0 + (-avx) e1 + (-avy) e2 + (-avz) e3", to_string(-a1));
	EXPECT_EQ("(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3", to_string(a1 + b1));
	EXPECT_EQ("(ae0-be0) e0 + (avx-bvx) e1 + (avy-bvy) e2 + (avz-bvz) e3", to_string(a1 - b1));
	EXPECT_EQ("(ae0*bs) e0 + (avx*bs) e1 + (avy*bs) e2 + (avz*bs) e3", to_string(a1 * symbolic("bs")));
	EXPECT_EQ("(ae0/bs) e0 + (avx/bs) e1 + (avy/bs) e2 + (avz/bs) e3", to_string(a1 / symbolic("bs")));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3", to_string(a1.rev()));
}

TEST(pga_test, operators_inplace_blade1) {
	auto a1 = pga::blade1<symbolic>({ "ae0" }, { {"avx"}, {"avy"}, {"avz"} });
	auto b1 = pga::blade1<symbolic>({ "be0" }, { {"bvx"}, {"bvy"}, {"bvz"} });
	auto r = a1; r += b1;
	EXPECT_EQ("(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3", to_string(r));
	r = a1; r -= b1;
	EXPECT_EQ("(ae0-be0) e0 + (avx-bvx) e1 + (avy-bvy) e2 + (avz-bvz) e3", to_string(r));
	r = a1; r *= symbolic("bs");
	EXPECT_EQ("(ae0*bs) e0 + (avx*bs) e1 + (avy*bs) e2 + (avz*bs) e3", to_string(r));
	r = a1; r /= symbolic("bs");
	EXPECT_EQ("(ae0/bs) e0 + (avx/bs) e1 + (avy/bs) e2 + (avz/bs) e3", to_string(r));
	r = a1; r += a1;
	EXPECT_EQ("(ae0+ae0) e0 + (avx+avx) e1 + (avy+avy) e2 + (avz+avz) e3", to_string(r));
	r = a1; r -= a1;
	EXPECT_EQ("(ae0-ae0) e0 + (avx-avx) e1 + (avy-avy) e2 + (avz-avz) e3", to_string(r));
}

TEST(pga_test, operators_multiply) {
    auto a1 = pga::blade1<symbolic>({ "ae0" }, { {"avx"}, {"avy"}, {"avz"} });
    auto b1 = pga::blade1<symbolic>({ "be0" }, { {"bvx"}, {"bvy"}, {"bvz"} });
    EXPECT_EQ(
        "(((avx*bvx)+(avy*bvy))+(avz*bvz)) id + "
        "((avy*bvz)-(avz*bvy)) e23 + "
        "((avz*bvx)-(avx*bvz)) e31 + "
        "((avx*bvy)-(avy*bvx)) e12 + "
        "((ae0*bvx)-(be0*avx)) e01 + "
        "((ae0*bvy)-(be0*avy)) e02 + "
        "((ae0*bvz)-(be0*avz)) e03 + "
        "0 e0123",
        to_string(a1 * b1));
}

//TEST(symbolic_test, casts) {
//    symbolic s("s");
//    symbolic e0 = zeroOf(s);
//    symbolic e1 = identityOf(s);
//    EXPECT_EQ("0", e0.v);
//    EXPECT_EQ("1", e1.v);
//    symbolic s5 = castOf<symbolic>(5);
//    EXPECT_EQ("5", s5.v);
//}
