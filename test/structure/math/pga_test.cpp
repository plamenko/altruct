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
	EXPECT_EQ("0", d1.v.x.v);
	EXPECT_EQ("0", d1.v.y.v);
	EXPECT_EQ("0", d1.v.z.v);
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
	auto a02 = pga::blade02<symbolic>({ "as" }, { {"abiEx"}, {"abiEy"}, {"abiEz"} });
	auto a24 = pga::blade24<symbolic>({ {"abiex"}, {"abiey"}, {"abiez"} }, { "ae0123" });
	auto a3 = pga::blade3<symbolic>({ "ae123" }, { {"atriPx"}, {"atriPy"}, {"atriPz"} });
	auto b1 = pga::blade1<symbolic>({ "be0" }, { {"bvx"}, {"bvy"}, {"bvz"} });
	auto b02 = pga::blade02<symbolic>({ "bs" }, { {"bbiEx"}, {"bbiEy"}, {"bbiEz"} });
	auto b24 = pga::blade24<symbolic>({ {"bbiex"}, {"bbiey"}, {"bbiez"} }, { "be0123" });
	auto b3 = pga::blade3<symbolic>({ "be123" }, { {"btriPx"}, {"btriPy"}, {"btriPz"} });
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
	EXPECT_EQ(
		"(ae0*bs) e0 + "
		"((avx*bs)-((avy*bbiEz)-(avz*bbiEy))) e1 + "
		"((avy*bs)-((avz*bbiEx)-(avx*bbiEz))) e2 + "
		"((avz*bs)-((avx*bbiEy)-(avy*bbiEx))) e3 + "
		"(((avx*bbiEx)+(avy*bbiEy))+(avz*bbiEz)) e123 + "
		"((-ae0)*bbiEx) e032 + "
		"((-ae0)*bbiEy) e013 + "
		"((-ae0)*bbiEz) e021",
		to_string(a1 * b02));
	EXPECT_EQ(
		"((((-avx)*bbiex)+((-avy)*bbiey))+((-avz)*bbiez)) e0 + "
		"0 e1 + "
		"0 e2 + "
		"0 e3 + "
		"0 e123 + "
		"((avx*be0123)+((avy*bbiez)-(avz*bbiey))) e032 + "
		"((avy*be0123)+((avz*bbiex)-(avx*bbiez))) e013 + "
		"((avz*be0123)+((avx*bbiey)-(avy*bbiex))) e021",
		to_string(a1 * b24));
	EXPECT_EQ(
		"0 id + "
		"(avx*be123) e23 + "
		"(avy*be123) e31 + "
		"(avz*be123) e12 + "
		"(((-avy)*btriPz)-((-avz)*btriPy)) e01 + "
		"(((-avz)*btriPx)-((-avx)*btriPz)) e02 + "
		"(((-avx)*btriPy)-((-avy)*btriPx)) e03 + "
		"((ae0*be123)+(((avx*btriPx)+(avy*btriPy))+(avz*btriPz))) e0123",
		to_string(a1 * b3));
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
