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
	auto v1 = pga::blade1<symbolic>({ {"avx"}, {"avy"}, {"avz"} });
	EXPECT_EQ("0", v1.e0.v);
	EXPECT_EQ("avx", v1.v.x.v);
	EXPECT_EQ("avy", v1.v.y.v);
	EXPECT_EQ("avz", v1.v.z.v);
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
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3", to_string(~a1));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3", to_string(a1.rev()));
	EXPECT_EQ("ae0 e123 + avx e032 + avy e013 + avz e021", to_string(!a1));
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

TEST(pga_test, constructor_blade02) {
	auto d02 = pga::blade02<symbolic>();
	EXPECT_EQ("?", d02.s.v);
	EXPECT_EQ("0", d02.biE.x.v);
	EXPECT_EQ("0", d02.biE.y.v);
	EXPECT_EQ("0", d02.biE.z.v);
	auto s02 = pga::blade02<symbolic>({ "ss" });
	EXPECT_EQ("ss", s02.s.v);
	EXPECT_EQ("0", s02.biE.x.v);
	EXPECT_EQ("0", s02.biE.y.v);
	EXPECT_EQ("0", s02.biE.z.v);
	auto v02 = pga::blade02<symbolic>({ {"abiEx"}, {"abiEy"}, {"abiEz"} });
	EXPECT_EQ("0", v02.s.v);
	EXPECT_EQ("abiEx", v02.biE.x.v);
	EXPECT_EQ("abiEy", v02.biE.y.v);
	EXPECT_EQ("abiEz", v02.biE.z.v);
	auto a02 = pga::blade02<symbolic>({ "as" }, { {"abiEx"}, {"abiEy"}, {"abiEz"} });
	EXPECT_EQ("as", a02.s.v);
	EXPECT_EQ("abiEx", a02.biE.x.v);
	EXPECT_EQ("abiEy", a02.biE.y.v);
	EXPECT_EQ("abiEz", a02.biE.z.v);
}

TEST(pga_test, operators_arithmetic_blade02) {
	auto a02 = pga::blade02<symbolic>({ "as" }, { {"abiEx"}, {"abiEy"}, {"abiEz"} });
	auto b02 = pga::blade02<symbolic>({ "bs" }, { {"bbiEx"}, {"bbiEy"}, {"bbiEz"} });
	EXPECT_EQ("(+as) id + (+abiEx) e23 + (+abiEy) e31 + (+abiEz) e12", to_string(+a02));
	EXPECT_EQ("(-as) id + (-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12", to_string(-a02));
	EXPECT_EQ("(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12", to_string(a02 + b02));
	EXPECT_EQ("(as-bs) id + (abiEx-bbiEx) e23 + (abiEy-bbiEy) e31 + (abiEz-bbiEz) e12", to_string(a02 - b02));
	EXPECT_EQ("(as*bs) id + (abiEx*bs) e23 + (abiEy*bs) e31 + (abiEz*bs) e12", to_string(a02 * symbolic("bs")));
	EXPECT_EQ("(as/bs) id + (abiEx/bs) e23 + (abiEy/bs) e31 + (abiEz/bs) e12", to_string(a02 / symbolic("bs")));
	EXPECT_EQ("as id + (-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12", to_string(~a02));
	EXPECT_EQ("as id + (-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12", to_string(a02.rev()));
	EXPECT_EQ("abiEx e01 + abiEy e02 + abiEz e03 + as e0123", to_string(!a02));
}

TEST(pga_test, operators_inplace_blade02) {
	auto a02 = pga::blade02<symbolic>({ "as" }, { {"abiEx"}, {"abiEy"}, {"abiEz"} });
	auto b02 = pga::blade02<symbolic>({ "bs" }, { {"bbiEx"}, {"bbiEy"}, {"bbiEz"} });
	auto r = a02; r += b02;
	EXPECT_EQ("(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12", to_string(r));
	r = a02; r -= b02;
	EXPECT_EQ("(as-bs) id + (abiEx-bbiEx) e23 + (abiEy-bbiEy) e31 + (abiEz-bbiEz) e12", to_string(r));
	r = a02; r *= symbolic("bs");
	EXPECT_EQ("(as*bs) id + (abiEx*bs) e23 + (abiEy*bs) e31 + (abiEz*bs) e12", to_string(r));
	r = a02; r /= symbolic("bs");
	EXPECT_EQ("(as/bs) id + (abiEx/bs) e23 + (abiEy/bs) e31 + (abiEz/bs) e12", to_string(r));
	r = a02; r += a02;
	EXPECT_EQ("(as+as) id + (abiEx+abiEx) e23 + (abiEy+abiEy) e31 + (abiEz+abiEz) e12", to_string(r));
	r = a02; r -= a02;
	EXPECT_EQ("(as-as) id + (abiEx-abiEx) e23 + (abiEy-abiEy) e31 + (abiEz-abiEz) e12", to_string(r));
}

TEST(pga_test, constructor_blade24) {
	auto d24 = pga::blade24<symbolic>();
	EXPECT_EQ("0", d24.bie.x.v);
	EXPECT_EQ("0", d24.bie.y.v);
	EXPECT_EQ("0", d24.bie.z.v);
	EXPECT_EQ("?", d24.e0123.v);
	auto s24 = pga::blade24<symbolic>({ "se0123" });
	EXPECT_EQ("0", s24.bie.x.v);
	EXPECT_EQ("0", s24.bie.y.v);
	EXPECT_EQ("0", s24.bie.z.v);
	EXPECT_EQ("se0123", s24.e0123.v);
	auto v24 = pga::blade24<symbolic>({ {"abiex"}, {"abiey"}, {"abiez"} });
	EXPECT_EQ("abiex", v24.bie.x.v);
	EXPECT_EQ("abiey", v24.bie.y.v);
	EXPECT_EQ("abiez", v24.bie.z.v);
	EXPECT_EQ("0", v24.e0123.v);
	auto a24 = pga::blade24<symbolic>({ {"abiex"}, {"abiey"}, {"abiez"} }, { "ae0123" });
	EXPECT_EQ("abiex", a24.bie.x.v);
	EXPECT_EQ("abiey", a24.bie.y.v);
	EXPECT_EQ("abiez", a24.bie.z.v);
	EXPECT_EQ("ae0123", a24.e0123.v);
}

TEST(pga_test, operators_arithmetic_blade24) {
	auto a24 = pga::blade24<symbolic>({ {"abiex"}, {"abiey"}, {"abiez"} }, { "ae0123" });
	auto b24 = pga::blade24<symbolic>({ {"bbiex"}, {"bbiey"}, {"bbiez"} }, { "be0123" });
	EXPECT_EQ("(+abiex) e01 + (+abiey) e02 + (+abiez) e03 + (+ae0123) e0123", to_string(+a24));
	EXPECT_EQ("(-abiex) e01 + (-abiey) e02 + (-abiez) e03 + (-ae0123) e0123", to_string(-a24));
	EXPECT_EQ("(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123", to_string(a24 + b24));
	EXPECT_EQ("(abiex-bbiex) e01 + (abiey-bbiey) e02 + (abiez-bbiez) e03 + (ae0123-be0123) e0123", to_string(a24 - b24));
	EXPECT_EQ("(abiex*bs) e01 + (abiey*bs) e02 + (abiez*bs) e03 + (ae0123*bs) e0123", to_string(a24 * symbolic("bs")));
	EXPECT_EQ("(abiex/bs) e01 + (abiey/bs) e02 + (abiez/bs) e03 + (ae0123/bs) e0123", to_string(a24 / symbolic("bs")));
	EXPECT_EQ("(-abiex) e01 + (-abiey) e02 + (-abiez) e03 + ae0123 e0123", to_string(~a24));
	EXPECT_EQ("(-abiex) e01 + (-abiey) e02 + (-abiez) e03 + ae0123 e0123", to_string(a24.rev()));
	EXPECT_EQ("ae0123 id + abiex e23 + abiey e31 + abiez e12", to_string(!a24));
}

TEST(pga_test, operators_inplace_blade24) {
	auto a24 = pga::blade24<symbolic>({ {"abiex"}, {"abiey"}, {"abiez"} }, { "ae0123" });
	auto b24 = pga::blade24<symbolic>({ {"bbiex"}, {"bbiey"}, {"bbiez"} }, { "be0123" });
	auto r = a24; r += b24;
	EXPECT_EQ("(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123", to_string(r));
	r = a24; r -= b24;
	EXPECT_EQ("(abiex-bbiex) e01 + (abiey-bbiey) e02 + (abiez-bbiez) e03 + (ae0123-be0123) e0123", to_string(r));
	r = a24; r *= symbolic("bs");
	EXPECT_EQ("(abiex*bs) e01 + (abiey*bs) e02 + (abiez*bs) e03 + (ae0123*bs) e0123", to_string(r));
	r = a24; r /= symbolic("bs");
	EXPECT_EQ("(abiex/bs) e01 + (abiey/bs) e02 + (abiez/bs) e03 + (ae0123/bs) e0123", to_string(r));
	r = a24; r += a24;
	EXPECT_EQ("(abiex+abiex) e01 + (abiey+abiey) e02 + (abiez+abiez) e03 + (ae0123+ae0123) e0123", to_string(r));
	r = a24; r -= a24;
	EXPECT_EQ("(abiex-abiex) e01 + (abiey-abiey) e02 + (abiez-abiez) e03 + (ae0123-ae0123) e0123", to_string(r));
}

TEST(pga_test, constructor_blade3) {
	auto d3 = pga::blade3<symbolic>();
	EXPECT_EQ("?", d3.e123.v);
	EXPECT_EQ("0", d3.triP.x.v);
	EXPECT_EQ("0", d3.triP.y.v);
	EXPECT_EQ("0", d3.triP.z.v);
	auto s3 = pga::blade3<symbolic>({ "se123" });
	EXPECT_EQ("se123", s3.e123.v);
	EXPECT_EQ("0", s3.triP.x.v);
	EXPECT_EQ("0", s3.triP.y.v);
	EXPECT_EQ("0", s3.triP.z.v);
	auto v3 = pga::blade3<symbolic>({ {"atriPx"}, {"atriPy"}, {"atriPz"} });
	EXPECT_EQ("0", v3.e123.v);
	EXPECT_EQ("atriPx", v3.triP.x.v);
	EXPECT_EQ("atriPy", v3.triP.y.v);
	EXPECT_EQ("atriPz", v3.triP.z.v);
	auto a3 = pga::blade3<symbolic>({ "ae123" }, { {"atriPx"}, {"atriPy"}, {"atriPz"} });
	EXPECT_EQ("ae123", a3.e123.v);
	EXPECT_EQ("atriPx", a3.triP.x.v);
	EXPECT_EQ("atriPy", a3.triP.y.v);
	EXPECT_EQ("atriPz", a3.triP.z.v);
}

TEST(pga_test, operators_arithmetic_blade3) {
	auto a3 = pga::blade3<symbolic>({ "ae123" }, { {"atriPx"}, {"atriPy"}, {"atriPz"} });
	auto b3 = pga::blade3<symbolic>({ "be123" }, { {"btriPx"}, {"btriPy"}, {"btriPz"} });
	EXPECT_EQ("(+ae123) e123 + (+atriPx) e032 + (+atriPy) e013 + (+atriPz) e021", to_string(+a3));
	EXPECT_EQ("(-ae123) e123 + (-atriPx) e032 + (-atriPy) e013 + (-atriPz) e021", to_string(-a3));
	EXPECT_EQ("(ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021", to_string(a3 + b3));
	EXPECT_EQ("(ae123-be123) e123 + (atriPx-btriPx) e032 + (atriPy-btriPy) e013 + (atriPz-btriPz) e021", to_string(a3 - b3));
	EXPECT_EQ("(ae123*bs) e123 + (atriPx*bs) e032 + (atriPy*bs) e013 + (atriPz*bs) e021", to_string(a3 * symbolic("bs")));
	EXPECT_EQ("(ae123/bs) e123 + (atriPx/bs) e032 + (atriPy/bs) e013 + (atriPz/bs) e021", to_string(a3 / symbolic("bs")));
	EXPECT_EQ("(-ae123) e123 + (-atriPx) e032 + (-atriPy) e013 + (-atriPz) e021", to_string(~a3));
	EXPECT_EQ("(-ae123) e123 + (-atriPx) e032 + (-atriPy) e013 + (-atriPz) e021", to_string(a3.rev()));
	EXPECT_EQ("ae123 e0 + atriPx e1 + atriPy e2 + atriPz e3", to_string(!a3));
}

TEST(pga_test, operators_inplace_blade3) {
	auto a3 = pga::blade3<symbolic>({ "ae123" }, { {"atriPx"}, {"atriPy"}, {"atriPz"} });
	auto b3 = pga::blade3<symbolic>({ "be123" }, { {"btriPx"}, {"btriPy"}, {"btriPz"} });
	auto r = a3; r += b3;
	EXPECT_EQ("(ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021", to_string(r));
	r = a3; r -= b3;
	EXPECT_EQ("(ae123-be123) e123 + (atriPx-btriPx) e032 + (atriPy-btriPy) e013 + (atriPz-btriPz) e021", to_string(r));
	r = a3; r *= symbolic("bs");
	EXPECT_EQ("(ae123*bs) e123 + (atriPx*bs) e032 + (atriPy*bs) e013 + (atriPz*bs) e021", to_string(r));
	r = a3; r /= symbolic("bs");
	EXPECT_EQ("(ae123/bs) e123 + (atriPx/bs) e032 + (atriPy/bs) e013 + (atriPz/bs) e021", to_string(r));
	r = a3; r += a3;
	EXPECT_EQ("(ae123+ae123) e123 + (atriPx+atriPx) e032 + (atriPy+atriPy) e013 + (atriPz+atriPz) e021", to_string(r));
	r = a3; r -= a3;
	EXPECT_EQ("(ae123-ae123) e123 + (atriPx-atriPx) e032 + (atriPy-atriPy) e013 + (atriPz-atriPz) e021", to_string(r));
}

TEST(pga_test, constructor_multivector) {
	auto dm = pga::multivector<symbolic>();
	EXPECT_EQ("? e0 + 0 e1 + 0 e2 + 0 e3 + ? e123 + 0 e032 + 0 e013 + 0 e021 + ? id + 0 e23 + 0 e31 + 0 e12 + 0 e01 + 0 e02 + 0 e03 + ? e0123", to_string(dm));
	auto a1 = pga::blade1<symbolic>({ "ae0" }, { {"avx"}, {"avy"}, {"avz"} });
	auto a02 = pga::blade02<symbolic>({ "as" }, { {"abiEx"}, {"abiEy"}, {"abiEz"} });
	auto a24 = pga::blade24<symbolic>({ {"abiex"}, {"abiey"}, {"abiez"} }, { "ae0123" });
	auto a3 = pga::blade3<symbolic>({ "ae123" }, { {"atriPx"}, {"atriPy"}, {"atriPz"} });
	auto am = pga::multivector<symbolic>({ a1, a3 }, { a02, a24 });
	EXPECT_EQ(
		"ae0 e0 + avx e1 + avy e2 + avz e3 + "
		"ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021 + "
		"as id + abiEx e23 + abiEy e31 + abiEz e12 + "
		"abiex e01 + abiey e02 + abiez e03 + ae0123 e0123",
		to_string(am));
	auto an = pga::multivector<symbolic>(a1, a02, a24, a3);
	EXPECT_EQ(
		"ae0 e0 + avx e1 + avy e2 + avz e3 + "
		"ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021 + "
		"as id + abiEx e23 + abiEy e31 + abiEz e12 + "
		"abiex e01 + abiey e02 + abiez e03 + ae0123 e0123",
		to_string(an));
}

TEST(pga_test, operators_arithmetic_multivector) {
	auto a1 = pga::blade1<symbolic>({ "ae0" }, { {"avx"}, {"avy"}, {"avz"} });
	auto a02 = pga::blade02<symbolic>({ "as" }, { {"abiEx"}, {"abiEy"}, {"abiEz"} });
	auto a24 = pga::blade24<symbolic>({ {"abiex"}, {"abiey"}, {"abiez"} }, { "ae0123" });
	auto a3 = pga::blade3<symbolic>({ "ae123" }, { {"atriPx"}, {"atriPy"}, {"atriPz"} });
	auto am = pga::multivector<symbolic>(a1, a02, a24, a3);
	auto b1 = pga::blade1<symbolic>({ "be0" }, { {"bvx"}, {"bvy"}, {"bvz"} });
	auto b02 = pga::blade02<symbolic>({ "bs" }, { {"bbiEx"}, {"bbiEy"}, {"bbiEz"} });
	auto b24 = pga::blade24<symbolic>({ {"bbiex"}, {"bbiey"}, {"bbiez"} }, { "be0123" });
	auto b3 = pga::blade3<symbolic>({ "be123" }, { {"btriPx"}, {"btriPy"}, {"btriPz"} });
	auto bm = pga::multivector<symbolic>(b1, b02, b24, b3);
	EXPECT_EQ(
		"(+ae0) e0 + (+avx) e1 + (+avy) e2 + (+avz) e3 + "
		"(+ae123) e123 + (+atriPx) e032 + (+atriPy) e013 + (+atriPz) e021 + "
		"(+as) id + (+abiEx) e23 + (+abiEy) e31 + (+abiEz) e12 + "
		"(+abiex) e01 + (+abiey) e02 + (+abiez) e03 + (+ae0123) e0123",
		to_string(+am));
	EXPECT_EQ(
		"(-ae0) e0 + (-avx) e1 + (-avy) e2 + (-avz) e3 + "
		"(-ae123) e123 + (-atriPx) e032 + (-atriPy) e013 + (-atriPz) e021 + "
		"(-as) id + (-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12 + "
		"(-abiex) e01 + (-abiey) e02 + (-abiez) e03 + (-ae0123) e0123",
		to_string(-am));
	EXPECT_EQ(
		"(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3 + "
		"(ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021 + "
		"(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + "
		"(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123",
		to_string(am + bm));
	EXPECT_EQ(
		"(ae0-be0) e0 + (avx-bvx) e1 + (avy-bvy) e2 + (avz-bvz) e3 + "
		"(ae123-be123) e123 + (atriPx-btriPx) e032 + (atriPy-btriPy) e013 + (atriPz-btriPz) e021 + "
		"(as-bs) id + (abiEx-bbiEx) e23 + (abiEy-bbiEy) e31 + (abiEz-bbiEz) e12 + "
		"(abiex-bbiex) e01 + (abiey-bbiey) e02 + (abiez-bbiez) e03 + (ae0123-be0123) e0123",
		to_string(am - bm));
	EXPECT_EQ(
		"(ae0*bs) e0 + (avx*bs) e1 + (avy*bs) e2 + (avz*bs) e3 + "
		"(ae123*bs) e123 + (atriPx*bs) e032 + (atriPy*bs) e013 + (atriPz*bs) e021 + "
		"(as*bs) id + (abiEx*bs) e23 + (abiEy*bs) e31 + (abiEz*bs) e12 + "
		"(abiex*bs) e01 + (abiey*bs) e02 + (abiez*bs) e03 + (ae0123*bs) e0123",
		to_string(am * symbolic("bs")));
	EXPECT_EQ(
		"(ae0/bs) e0 + (avx/bs) e1 + (avy/bs) e2 + (avz/bs) e3 + "
		"(ae123/bs) e123 + (atriPx/bs) e032 + (atriPy/bs) e013 + (atriPz/bs) e021 + "
		"(as/bs) id + (abiEx/bs) e23 + (abiEy/bs) e31 + (abiEz/bs) e12 + "
		"(abiex/bs) e01 + (abiey/bs) e02 + (abiez/bs) e03 + (ae0123/bs) e0123",
		to_string(am / symbolic("bs")));
	EXPECT_EQ(
		"ae0 e0 + avx e1 + avy e2 + avz e3 + "
		"(-ae123) e123 + (-atriPx) e032 + (-atriPy) e013 + (-atriPz) e021 + "
		"as id + (-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12 + "
		"(-abiex) e01 + (-abiey) e02 + (-abiez) e03 + ae0123 e0123",
		to_string(~am));
	EXPECT_EQ(
		"ae0 e0 + avx e1 + avy e2 + avz e3 + "
		"(-ae123) e123 + (-atriPx) e032 + (-atriPy) e013 + (-atriPz) e021 + "
		"as id + (-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12 + "
		"(-abiex) e01 + (-abiey) e02 + (-abiez) e03 + ae0123 e0123",
		to_string(am.rev()));
	EXPECT_EQ(
		"ae123 e0 + atriPx e1 + atriPy e2 + atriPz e3 + "
		"ae0 e123 + avx e032 + avy e013 + avz e021 + "
		"ae0123 id + abiex e23 + abiey e31 + abiez e12 + "
		"abiEx e01 + abiEy e02 + abiEz e03 + as e0123",
		to_string(!am));
}

TEST(pga_test, operators_inplace_multivector) {
	auto a1 = pga::blade1<symbolic>({ "ae0" }, { {"avx"}, {"avy"}, {"avz"} });
	auto a02 = pga::blade02<symbolic>({ "as" }, { {"abiEx"}, {"abiEy"}, {"abiEz"} });
	auto a24 = pga::blade24<symbolic>({ {"abiex"}, {"abiey"}, {"abiez"} }, { "ae0123" });
	auto a3 = pga::blade3<symbolic>({ "ae123" }, { {"atriPx"}, {"atriPy"}, {"atriPz"} });
	auto am = pga::multivector<symbolic>(a1, a02, a24, a3);
	auto b1 = pga::blade1<symbolic>({ "be0" }, { {"bvx"}, {"bvy"}, {"bvz"} });
	auto b02 = pga::blade02<symbolic>({ "bs" }, { {"bbiEx"}, {"bbiEy"}, {"bbiEz"} });
	auto b24 = pga::blade24<symbolic>({ {"bbiex"}, {"bbiey"}, {"bbiez"} }, { "be0123" });
	auto b3 = pga::blade3<symbolic>({ "be123" }, { {"btriPx"}, {"btriPy"}, {"btriPz"} });
	auto bm = pga::multivector<symbolic>(b1, b02, b24, b3);
	auto r = am; r += bm;
	EXPECT_EQ(
		"(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3 + "
		"(ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021 + "
		"(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + "
		"(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123",
		to_string(r));
	r = am; r -= bm;
	EXPECT_EQ(
		"(ae0-be0) e0 + (avx-bvx) e1 + (avy-bvy) e2 + (avz-bvz) e3 + "
		"(ae123-be123) e123 + (atriPx-btriPx) e032 + (atriPy-btriPy) e013 + (atriPz-btriPz) e021 + "
		"(as-bs) id + (abiEx-bbiEx) e23 + (abiEy-bbiEy) e31 + (abiEz-bbiEz) e12 + "
		"(abiex-bbiex) e01 + (abiey-bbiey) e02 + (abiez-bbiez) e03 + (ae0123-be0123) e0123",
		to_string(r));
	r = am; r *= symbolic("bs");
	EXPECT_EQ(
		"(ae0*bs) e0 + (avx*bs) e1 + (avy*bs) e2 + (avz*bs) e3 + "
		"(ae123*bs) e123 + (atriPx*bs) e032 + (atriPy*bs) e013 + (atriPz*bs) e021 + "
		"(as*bs) id + (abiEx*bs) e23 + (abiEy*bs) e31 + (abiEz*bs) e12 + "
		"(abiex*bs) e01 + (abiey*bs) e02 + (abiez*bs) e03 + (ae0123*bs) e0123",
		to_string(r));
	r = am; r /= symbolic("bs");
	EXPECT_EQ(
		"(ae0/bs) e0 + (avx/bs) e1 + (avy/bs) e2 + (avz/bs) e3 + "
		"(ae123/bs) e123 + (atriPx/bs) e032 + (atriPy/bs) e013 + (atriPz/bs) e021 + "
		"(as/bs) id + (abiEx/bs) e23 + (abiEy/bs) e31 + (abiEz/bs) e12 + "
		"(abiex/bs) e01 + (abiey/bs) e02 + (abiez/bs) e03 + (ae0123/bs) e0123",
		to_string(r));
	r = am; r += am;
	EXPECT_EQ(
		"(ae0+ae0) e0 + (avx+avx) e1 + (avy+avy) e2 + (avz+avz) e3 + "
		"(ae123+ae123) e123 + (atriPx+atriPx) e032 + (atriPy+atriPy) e013 + (atriPz+atriPz) e021 + "
		"(as+as) id + (abiEx+abiEx) e23 + (abiEy+abiEy) e31 + (abiEz+abiEz) e12 + "
		"(abiex+abiex) e01 + (abiey+abiey) e02 + (abiez+abiez) e03 + (ae0123+ae0123) e0123",
		to_string(r));
	r = am; r -= am;
	EXPECT_EQ(
		"(ae0-ae0) e0 + (avx-avx) e1 + (avy-avy) e2 + (avz-avz) e3 + "
		"(ae123-ae123) e123 + (atriPx-atriPx) e032 + (atriPy-atriPy) e013 + (atriPz-atriPz) e021 + "
		"(as-as) id + (abiEx-abiEx) e23 + (abiEy-abiEy) e31 + (abiEz-abiEz) e12 + "
		"(abiex-abiex) e01 + (abiey-abiey) e02 + (abiez-abiez) e03 + (ae0123-ae0123) e0123",
		to_string(r));
}

TEST(pga_test, operators_dual) {
	auto a1 = pga::blade1<symbolic>({ "ae0" }, { {"avx"}, {"avy"}, {"avz"} });
	auto a02 = pga::blade02<symbolic>({ "as" }, { {"abiEx"}, {"abiEy"}, {"abiEz"} });
	auto a24 = pga::blade24<symbolic>({ {"abiex"}, {"abiey"}, {"abiez"} }, { "ae0123" });
	auto a3 = pga::blade3<symbolic>({ "ae123" }, { {"atriPx"}, {"atriPy"}, {"atriPz"} });
	auto am = pga::multivector<symbolic>(a1, a02, a24, a3);
	EXPECT_EQ("ae0 e123 + avx e032 + avy e013 + avz e021", to_string(!a1));
	EXPECT_EQ("abiEx e01 + abiEy e02 + abiEz e03 + as e0123", to_string(!a02));
	EXPECT_EQ("ae0123 id + abiex e23 + abiey e31 + abiez e12", to_string(!a24));
	EXPECT_EQ("ae123 e0 + atriPx e1 + atriPy e2 + atriPz e3", to_string(!a3));
	EXPECT_EQ(
		"ae123 e0 + atriPx e1 + atriPy e2 + atriPz e3 + "
		"ae0 e123 + avx e032 + avy e013 + avz e021 + "
		"ae0123 id + abiex e23 + abiey e31 + abiez e12 + "
		"abiEx e01 + abiEy e02 + abiEz e03 + as e0123",
		to_string(!am));
}

TEST(pga_test, get) {
	auto z = pga::zero<symbolic>();
	auto a1 = pga::blade1<symbolic>({ "ae0" }, { {"avx"}, {"avy"}, {"avz"} });
	auto a02 = pga::blade02<symbolic>({ "as" }, { {"abiEx"}, {"abiEy"}, {"abiEz"} });
	auto a24 = pga::blade24<symbolic>({ {"abiex"}, {"abiey"}, {"abiez"} }, { "ae0123" });
	auto a3 = pga::blade3<symbolic>({ "ae123" }, { {"atriPx"}, {"atriPy"}, {"atriPz"} });
	auto a13 = pga::blade13<symbolic>(a1, a3);
	auto a024 = pga::blade024<symbolic>(a02, a24);
	auto am = pga::multivector<symbolic>(a1, a02, a24, a3);

	EXPECT_EQ("0", to_string(pga::get<decltype(z)>::b1(z)));
	EXPECT_EQ("0", to_string(pga::get<decltype(z)>::b02(z)));
	EXPECT_EQ("0", to_string(pga::get<decltype(z)>::b24(z)));
	EXPECT_EQ("0", to_string(pga::get<decltype(z)>::b3(z)));

	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3", to_string(pga::get<decltype(a1)>::b1(a1)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a1)>::b02(a1)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a1)>::b24(a1)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a1)>::b3(a1)));

	EXPECT_EQ("0", to_string(pga::get<decltype(a02)>::b1(a02)));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12", to_string(pga::get<decltype(a02)>::b02(a02)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a02)>::b24(a02)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a02)>::b3(a02)));

	EXPECT_EQ("0", to_string(pga::get<decltype(a24)>::b1(a24)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a24)>::b02(a24)));
	EXPECT_EQ("abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(pga::get<decltype(a24)>::b24(a24)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a24)>::b3(a24)));

	EXPECT_EQ("0", to_string(pga::get<decltype(a3)>::b1(a3)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a3)>::b02(a3)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a3)>::b24(a3)));
	EXPECT_EQ("ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(pga::get<decltype(a3)>::b3(a3)));

	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3", to_string(pga::get<decltype(a13)>::b1(a13)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a13)>::b02(a13)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a13)>::b24(a13)));
	EXPECT_EQ("ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(pga::get<decltype(a13)>::b3(a13)));

	EXPECT_EQ("0", to_string(pga::get<decltype(a024)>::b1(a024)));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12", to_string(pga::get<decltype(a024)>::b02(a024)));
	EXPECT_EQ("abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(pga::get<decltype(a024)>::b24(a024)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a024)>::b3(a024)));

	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3", to_string(pga::get<decltype(am)>::b1(am)));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12", to_string(pga::get<decltype(am)>::b02(am)));
	EXPECT_EQ("abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(pga::get<decltype(am)>::b24(am)));
	EXPECT_EQ("ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(pga::get<decltype(am)>::b3(am)));
}

TEST(pga_test, combine) {
	auto z = pga::zero<symbolic>();
	auto a1 = pga::blade1<symbolic>({ "ae0" }, { {"avx"}, {"avy"}, {"avz"} });
	auto a02 = pga::blade02<symbolic>({ "as" }, { {"abiEx"}, {"abiEy"}, {"abiEz"} });
	auto a24 = pga::blade24<symbolic>({ {"abiex"}, {"abiey"}, {"abiez"} }, { "ae0123" });
	auto a3 = pga::blade3<symbolic>({ "ae123" }, { {"atriPx"}, {"atriPy"}, {"atriPz"} });
	auto b1 = pga::blade1<symbolic>({ "be0" }, { {"bvx"}, {"bvy"}, {"bvz"} });
	auto b02 = pga::blade02<symbolic>({ "bs" }, { {"bbiEx"}, {"bbiEy"}, {"bbiEz"} });
	auto b24 = pga::blade24<symbolic>({ {"bbiex"}, {"bbiey"}, {"bbiez"} }, { "be0123" });
	auto b3 = pga::blade3<symbolic>({ "be123" }, { {"btriPx"}, {"btriPy"}, {"btriPz"} });
	auto a13 = pga::blade13<symbolic>(a1, a3);
	auto a024 = pga::blade024<symbolic>(a02, a24);
	auto b13 = pga::blade13<symbolic>(b1, b3);
	auto b024 = pga::blade024<symbolic>(b02, b24);
	auto am = pga::multivector<symbolic>(a1, a02, a24, a3);
	auto bm = pga::multivector<symbolic>(b1, b02, b24, b3);

	std::string sa1 = "ae0 e0 + avx e1 + avy e2 + avz e3";
	std::string sa3 = "ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021";
	std::string sb02 = "bs id + bbiEx e23 + bbiEy e31 + bbiEz e12";
	std::string sb24 = "bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123";
	std::string sz1 = "0 e0 + 0 e1 + 0 e2 + 0 e3";
	std::string sz3 = "0 e123 + 0 e032 + 0 e013 + 0 e021";
	std::string sz02 = "0 id + 0 e23 + 0 e31 + 0 e12";
	std::string sz24 = "0 e01 + 0 e02 + 0 e03 + 0 e0123";

	EXPECT_EQ("0", to_string(pga::combine_primitive(z, z)));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3", to_string(pga::combine_primitive(a1, z)));
	EXPECT_EQ("be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(pga::combine_primitive(z, b3)));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12", to_string(pga::combine_primitive(a02, z)));
	EXPECT_EQ("bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(pga::combine_primitive(z, b24)));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(pga::combine_primitive(a1, b3)));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(pga::combine_primitive(a02, b24)));

	EXPECT_EQ("0", to_string(pga::combine_multivector(z, z)));
	EXPECT_EQ(sb02, to_string(pga::combine_multivector(z, b02)));
	EXPECT_EQ(sb24, to_string(pga::combine_multivector(z, b24)));
	EXPECT_EQ(sb02 + " + " + sb24, to_string(pga::combine_multivector(z, b024)));
	EXPECT_EQ(sa1, to_string(pga::combine_multivector(a1, z)));
	EXPECT_EQ(sa1 + " + " + sz3 + " + " + sb02 + " + " + sz24, to_string(pga::combine_multivector(a1, b02)));
	EXPECT_EQ(sa1 + " + " + sz3 + " + " + sz02 + " + " + sb24, to_string(pga::combine_multivector(a1, b24)));
	EXPECT_EQ(sa1 + " + " + sz3 + " + " + sb02 + " + " + sb24, to_string(pga::combine_multivector(a1, b024)));
	EXPECT_EQ(sa3, to_string(pga::combine_multivector(a3, z)));
	EXPECT_EQ(sz1 + " + " + sa3 + " + " + sb02 + " + " + sz24, to_string(pga::combine_multivector(a3, b02)));
	EXPECT_EQ(sz1 + " + " + sa3 + " + " + sz02 + " + " + sb24, to_string(pga::combine_multivector(a3, b24)));
	EXPECT_EQ(sz1 + " + " + sa3 + " + " + sb02 + " + " + sb24, to_string(pga::combine_multivector(a3, b024)));
	EXPECT_EQ(sa1 + " + " + sa3 + "", to_string(pga::combine_multivector(a13, z)));
	EXPECT_EQ(sa1 + " + " + sa3 + " + " + sb02 + " + " + sz24, to_string(pga::combine_multivector(a13, b02)));
	EXPECT_EQ(sa1 + " + " + sa3 + " + " + sz02 + " + " + sb24, to_string(pga::combine_multivector(a13, b24)));
	EXPECT_EQ(sa1 + " + " + sa3 + " + " + sb02 + " + " + sb24, to_string(pga::combine_multivector(a13, b024)));
}

TEST(pga_test, operators_add) {
	auto z = pga::zero<symbolic>();
	auto a1 = pga::blade1<symbolic>({ "ae0" }, { {"avx"}, {"avy"}, {"avz"} });
	auto a02 = pga::blade02<symbolic>({ "as" }, { {"abiEx"}, {"abiEy"}, {"abiEz"} });
	auto a24 = pga::blade24<symbolic>({ {"abiex"}, {"abiey"}, {"abiez"} }, { "ae0123" });
	auto a3 = pga::blade3<symbolic>({ "ae123" }, { {"atriPx"}, {"atriPy"}, {"atriPz"} });
	auto a13 = pga::blade13<symbolic>(a1, a3);
	auto a024 = pga::blade024<symbolic>(a02, a24);
	auto am = pga::multivector<symbolic>(a1, a02, a24, a3);
	auto b1 = pga::blade1<symbolic>({ "be0" }, { {"bvx"}, {"bvy"}, {"bvz"} });
	auto b02 = pga::blade02<symbolic>({ "bs" }, { {"bbiEx"}, {"bbiEy"}, {"bbiEz"} });
	auto b24 = pga::blade24<symbolic>({ {"bbiex"}, {"bbiey"}, {"bbiez"} }, { "be0123" });
	auto b3 = pga::blade3<symbolic>({ "be123" }, { {"btriPx"}, {"btriPy"}, {"btriPz"} });
	auto b13 = pga::blade13<symbolic>(b1, b3);
	auto b024 = pga::blade024<symbolic>(b02, b24);
	auto bm = pga::multivector<symbolic>(b1, b02, b24, b3);

	EXPECT_EQ("0", to_string(z + z));
	EXPECT_EQ("be0 e0 + bvx e1 + bvy e2 + bvz e3", to_string(z + b1));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12", to_string(z + b02));
	EXPECT_EQ("bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(z + b24));
	EXPECT_EQ("be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(z + b3));
	EXPECT_EQ("be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(z + b13));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(z + b024));
	EXPECT_EQ("be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021 + "
		"bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(z + bm));

	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3", to_string(a1 + z));
	EXPECT_EQ("(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3", to_string(a1 + b1));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3 + 0 e123 + 0 e032 + 0 e013 + 0 e021 + "
		"bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123", to_string(a1 + b02));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3 + 0 e123 + 0 e032 + 0 e013 + 0 e021 + 0 "
		"id + 0 e23 + 0 e31 + 0 e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(a1 + b24));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a1 + b3));
	EXPECT_EQ("(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a1 + b13));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3 + 0 e123 + 0 e032 + 0 e013 + 0 e021 + "
		"bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(a1 + b024));
	EXPECT_EQ("(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021 + "
		"bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(a1 + bm));

	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12", to_string(a02 + z));
	EXPECT_EQ("be0 e0 + bvx e1 + bvy e2 + bvz e3 + 0 e123 + 0 e032 + 0 e013 + 0 e021 + "
		"as id + abiEx e23 + abiEy e31 + abiEz e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123", to_string(a02 + b1));
	EXPECT_EQ("(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12", to_string(a02 + b02));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(a02 + b24));
	EXPECT_EQ("0 e0 + 0 e1 + 0 e2 + 0 e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021 + "
		"as id + abiEx e23 + abiEy e31 + abiEz e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123", to_string(a02 + b3));
	EXPECT_EQ("be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021 + "
		"as id + abiEx e23 + abiEy e31 + abiEz e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123", to_string(a02 + b13));
	EXPECT_EQ("(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(a02 + b024));
	EXPECT_EQ("be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021 + "
		"(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(a02 + bm));

	EXPECT_EQ("abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(a24 + z));
	EXPECT_EQ("be0 e0 + bvx e1 + bvy e2 + bvz e3 + 0 e123 + 0 e032 + 0 e013 + 0 e021 + "
		"0 id + 0 e23 + 0 e31 + 0 e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(a24 + b1));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(a24 + b02));
	EXPECT_EQ("(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123", to_string(a24 + b24));
	EXPECT_EQ("0 e0 + 0 e1 + 0 e2 + 0 e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021 + "
		"0 id + 0 e23 + 0 e31 + 0 e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(a24 + b3));
	EXPECT_EQ("be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021 + "
		"0 id + 0 e23 + 0 e31 + 0 e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(a24 + b13));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123", to_string(a24 + b024));
	EXPECT_EQ("be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021 + "
		"bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123", to_string(a24 + bm));

	EXPECT_EQ("ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a3 + z));
	EXPECT_EQ("be0 e0 + bvx e1 + bvy e2 + bvz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a3 + b1));
	EXPECT_EQ("0 e0 + 0 e1 + 0 e2 + 0 e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021 + "
		"bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123", to_string(a3 + b02));
	EXPECT_EQ("0 e0 + 0 e1 + 0 e2 + 0 e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021 + "
		"0 id + 0 e23 + 0 e31 + 0 e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(a3 + b24));
	EXPECT_EQ("(ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021", to_string(a3 + b3));
	EXPECT_EQ("be0 e0 + bvx e1 + bvy e2 + bvz e3 + (ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021", to_string(a3 + b13));
	EXPECT_EQ("0 e0 + 0 e1 + 0 e2 + 0 e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021 + "
		"bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(a3 + b024));
	EXPECT_EQ("be0 e0 + bvx e1 + bvy e2 + bvz e3 + (ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021 + "
		"bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(a3 + bm));

	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a13 + z));
	EXPECT_EQ("(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a13 + b1));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021 + "
		"bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123", to_string(a13 + b02));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021 + "
		"0 id + 0 e23 + 0 e31 + 0 e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(a13 + b24));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3 + (ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021", to_string(a13 + b3));
	EXPECT_EQ("(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3 + (ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021", to_string(a13 + b13));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021 + "
		"bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(a13 + b024));
	EXPECT_EQ("(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3 + (ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021 + "
		"bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(a13 + bm));

	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(a024 + z));
	EXPECT_EQ("be0 e0 + bvx e1 + bvy e2 + bvz e3 + 0 e123 + 0 e032 + 0 e013 + 0 e021 + "
		"as id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(a024 + b1));
	EXPECT_EQ("(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(a024 + b02));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123", to_string(a024 + b24));
	EXPECT_EQ("0 e0 + 0 e1 + 0 e2 + 0 e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021 + "
		"as id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(a024 + b3));
	EXPECT_EQ("be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021 + "
		"as id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(a024 + b13));
	EXPECT_EQ("(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + "
		"(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123", to_string(a024 + b024));
	EXPECT_EQ("be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021 + "
		"(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + "
		"(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123", to_string(a024 + bm));

	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021 + "
		"as id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(am + z));
	EXPECT_EQ("(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021 + "
		"as id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(am + b1));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021 + "
		"(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(am + b02));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021 + "
		"as id + abiEx e23 + abiEy e31 + abiEz e12 + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123", to_string(am + b24));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3 + (ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021 + "
		"as id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(am + b3));
	EXPECT_EQ("(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3 + "
		"(ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021 + "
		"as id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(am + b13));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021 + "
		"(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + "
		"(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123", to_string(am + b024));
	EXPECT_EQ("(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3 + "
		"(ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021 + "
		"(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + "
		"(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123", to_string(am + bm));
}

TEST(pga_test, operators_multiply) {
	auto as = symbolic("as");
    auto a1 = pga::blade1<symbolic>({ "ae0" }, { {"avx"}, {"avy"}, {"avz"} });
	auto a02 = pga::blade02<symbolic>({ "as" }, { {"abiEx"}, {"abiEy"}, {"abiEz"} });
	auto a24 = pga::blade24<symbolic>({ {"abiex"}, {"abiey"}, {"abiez"} }, { "ae0123" });
	auto a3 = pga::blade3<symbolic>({ "ae123" }, { {"atriPx"}, {"atriPy"}, {"atriPz"} });
	auto b1 = pga::blade1<symbolic>({ "be0" }, { {"bvx"}, {"bvy"}, {"bvz"} });
	auto b02 = pga::blade02<symbolic>({ "bs" }, { {"bbiEx"}, {"bbiEy"}, {"bbiEz"} });
	auto b24 = pga::blade24<symbolic>({ {"bbiex"}, {"bbiey"}, {"bbiez"} }, { "be0123" });
	auto b3 = pga::blade3<symbolic>({ "be123" }, { {"btriPx"}, {"btriPy"}, {"btriPz"} });
	auto am = pga::multivector<symbolic>(a1, a02, a24, a3);
	auto bm = pga::multivector<symbolic>(b1, b02, b24, b3);
	
	EXPECT_EQ("(as*be0) e0 + (as*bvx) e1 + (as*bvy) e2 + (as*bvz) e3", to_string(as * b1));
	EXPECT_EQ("(as*bs) id + (as*bbiEx) e23 + (as*bbiEy) e31 + (as*bbiEz) e12", to_string(as * b02));
	EXPECT_EQ("(as*bbiex) e01 + (as*bbiey) e02 + (as*bbiez) e03 + (as*be0123) e0123", to_string(as * b24));
	EXPECT_EQ("(as*be123) e123 + (as*btriPx) e032 + (as*btriPy) e013 + (as*btriPz) e021", to_string(as * b3));

	EXPECT_EQ(
		"(((avx*bvx)+(avy*bvy))+(avz*bvz)) id + "
		"((avy*bvz)-(avz*bvy)) e23 + "
		"((avz*bvx)-(avx*bvz)) e31 + "
		"((avx*bvy)-(avy*bvx)) e12 + "
		"((ae0*bvx)-(avx*be0)) e01 + "
		"((ae0*bvy)-(avy*be0)) e02 + "
		"((ae0*bvz)-(avz*be0)) e03 + "
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
		"(-(((avx*bbiex)+(avy*bbiey))+(avz*bbiez))) e0 + "
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

	EXPECT_EQ(
		"(as*be0) e0 + "
		"((as*bvx)-((abiEy*bvz)-(abiEz*bvy))) e1 + "
		"((as*bvy)-((abiEz*bvx)-(abiEx*bvz))) e2 + "
		"((as*bvz)-((abiEx*bvy)-(abiEy*bvx))) e3 + "
		"(((abiEx*bvx)+(abiEy*bvy))+(abiEz*bvz)) e123 + "
		"((-abiEx)*be0) e032 + "
		"((-abiEy)*be0) e013 + "
		"((-abiEz)*be0) e021",
		to_string(a02 * b1));
	EXPECT_EQ(
		"((as*bs)-(((abiEx*bbiEx)+(abiEy*bbiEy))+(abiEz*bbiEz))) id + "
		"(((abiEx*bs)+(as*bbiEx))-((abiEy*bbiEz)-(abiEz*bbiEy))) e23 + "
		"(((abiEy*bs)+(as*bbiEy))-((abiEz*bbiEx)-(abiEx*bbiEz))) e31 + "
		"(((abiEz*bs)+(as*bbiEz))-((abiEx*bbiEy)-(abiEy*bbiEx))) e12",
		to_string(a02 * b02));
	EXPECT_EQ(
		"(((as*bbiex)-(abiEx*be0123))-((abiEy*bbiez)-(abiEz*bbiey))) e01 + "
		"(((as*bbiey)-(abiEy*be0123))-((abiEz*bbiex)-(abiEx*bbiez))) e02 + "
		"(((as*bbiez)-(abiEz*be0123))-((abiEx*bbiey)-(abiEy*bbiex))) e03 + "
		"((as*be0123)+(((abiEx*bbiex)+(abiEy*bbiey))+(abiEz*bbiez))) e0123",
		to_string(a02 * b24));
	EXPECT_EQ(
		"(((abiEx*btriPx)+(abiEy*btriPy))+(abiEz*btriPz)) e0 + "
		"((-abiEx)*be123) e1 + "
		"((-abiEy)*be123) e2 + "
		"((-abiEz)*be123) e3 + "
		"(as*be123) e123 + "
		"((as*btriPx)-((abiEy*btriPz)-(abiEz*btriPy))) e032 + "
		"((as*btriPy)-((abiEz*btriPx)-(abiEx*btriPz))) e013 + "
		"((as*btriPz)-((abiEx*btriPy)-(abiEy*btriPx))) e021",
		to_string(a02 * b3));

	EXPECT_EQ(
		"(((abiex*bvx)+(abiey*bvy))+(abiez*bvz)) e0 + "
		"0 e1 + "
		"0 e2 + "
		"0 e3 + "
		"0 e123 + "
		"(((-ae0123)*bvx)-((abiey*bvz)-(abiez*bvy))) e032 + "
		"(((-ae0123)*bvy)-((abiez*bvx)-(abiex*bvz))) e013 + "
		"(((-ae0123)*bvz)-((abiex*bvy)-(abiey*bvx))) e021",
		to_string(a24 * b1));
	EXPECT_EQ(
		"(((abiex*bs)-(ae0123*bbiEx))-((abiey*bbiEz)-(abiez*bbiEy))) e01 + "
		"(((abiey*bs)-(ae0123*bbiEy))-((abiez*bbiEx)-(abiex*bbiEz))) e02 + "
		"(((abiez*bs)-(ae0123*bbiEz))-((abiex*bbiEy)-(abiey*bbiEx))) e03 + "
		"((ae0123*bs)+(((abiex*bbiEx)+(abiey*bbiEy))+(abiez*bbiEz))) e0123",
		to_string(a24 * b02));
	EXPECT_EQ(
		"0",
		to_string(a24 * b24));
	EXPECT_EQ(
		"((-ae0123)*be123) e0 + "
		"0 e1 + "
		"0 e2 + "
		"0 e3 + "
		"0 e123 + "
		"((-abiex)*be123) e032 + "
		"((-abiey)*be123) e013 + "
		"((-abiez)*be123) e021",
		to_string(a24 * b3));

	EXPECT_EQ(
		"0 id + "
		"(ae123*bvx) e23 + "
		"(ae123*bvy) e31 + "
		"(ae123*bvz) e12 + "
		"((atriPy*bvz)-(atriPz*bvy)) e01 + "
		"((atriPz*bvx)-(atriPx*bvz)) e02 + "
		"((atriPx*bvy)-(atriPy*bvx)) e03 + "
		"(((-ae123)*be0)-(((atriPx*bvx)+(atriPy*bvy))+(atriPz*bvz))) e0123",
		to_string(a3 * b1));
	EXPECT_EQ(
		"(((atriPx*bbiEx)+(atriPy*bbiEy))+(atriPz*bbiEz)) e0 + "
		"((-ae123)*bbiEx) e1 + "
		"((-ae123)*bbiEy) e2 + "
		"((-ae123)*bbiEz) e3 + "
		"(ae123*bs) e123 + "
		"((atriPx*bs)-((atriPy*bbiEz)-(atriPz*bbiEy))) e032 + "
		"((atriPy*bs)-((atriPz*bbiEx)-(atriPx*bbiEz))) e013 + "
		"((atriPz*bs)-((atriPx*bbiEy)-(atriPy*bbiEx))) e021",
		to_string(a3 * b02));
	EXPECT_EQ(
		"(ae123*be0123) e0 + "
		"0 e1 + "
		"0 e2 + "
		"0 e3 + "
		"0 e123 + "
		"(ae123*bbiex) e032 + "
		"(ae123*bbiey) e013 + "
		"(ae123*bbiez) e021",
		to_string(a3 * b24));
	EXPECT_EQ(
		"((-ae123)*be123) id + "
		"0 e23 + "
		"0 e31 + "
		"0 e12 + "
		"((atriPx*be123)-(ae123*btriPx)) e01 + "
		"((atriPy*be123)-(ae123*btriPy)) e02 + "
		"((atriPz*be123)-(ae123*btriPz)) e03 + "
		"0 e0123",
		to_string(a3 * b3));

	EXPECT_EQ(
		"((((as*be0)+(((abiex*bvx)+(abiey*bvy))+(abiez*bvz)))+((((abiEx*btriPx)+(abiEy*btriPy))+(abiEz*btriPz))+((-ae0123)*be123)))+(((ae0*bs)+(((atriPx*bbiEx)+(atriPy*bbiEy))+(atriPz*bbiEz)))+((-(((avx*bbiex)+(avy*bbiey))+(avz*bbiez)))+(ae123*be0123)))) e0 + "
		"(((((as*bvx)-((abiEy*bvz)-(abiEz*bvy)))+0)+(((-abiEx)*be123)+0))+((((avx*bs)-((avy*bbiEz)-(avz*bbiEy)))+((-ae123)*bbiEx))+(0+0))) e1 + "
		"(((((as*bvy)-((abiEz*bvx)-(abiEx*bvz)))+0)+(((-abiEy)*be123)+0))+((((avy*bs)-((avz*bbiEx)-(avx*bbiEz)))+((-ae123)*bbiEy))+(0+0))) e2 + "
		"(((((as*bvz)-((abiEx*bvy)-(abiEy*bvx)))+0)+(((-abiEz)*be123)+0))+((((avz*bs)-((avx*bbiEy)-(avy*bbiEx)))+((-ae123)*bbiEz))+(0+0))) e3 + "
		"((((((abiEx*bvx)+(abiEy*bvy))+(abiEz*bvz))+0)+((as*be123)+0))+(((((avx*bbiEx)+(avy*bbiEy))+(avz*bbiEz))+(ae123*bs))+(0+0))) e123 + "
		"(((((-abiEx)*be0)+(((-ae0123)*bvx)-((abiey*bvz)-(abiez*bvy))))+(((as*btriPx)-((abiEy*btriPz)-(abiEz*btriPy)))+((-abiex)*be123)))+((((-ae0)*bbiEx)+((atriPx*bs)-((atriPy*bbiEz)-(atriPz*bbiEy))))+(((avx*be0123)+((avy*bbiez)-(avz*bbiey)))+(ae123*bbiex)))) e032 + "
		"(((((-abiEy)*be0)+(((-ae0123)*bvy)-((abiez*bvx)-(abiex*bvz))))+(((as*btriPy)-((abiEz*btriPx)-(abiEx*btriPz)))+((-abiey)*be123)))+((((-ae0)*bbiEy)+((atriPy*bs)-((atriPz*bbiEx)-(atriPx*bbiEz))))+(((avy*be0123)+((avz*bbiex)-(avx*bbiez)))+(ae123*bbiey)))) e013 + "
		"(((((-abiEz)*be0)+(((-ae0123)*bvz)-((abiex*bvy)-(abiey*bvx))))+(((as*btriPz)-((abiEx*btriPy)-(abiEy*btriPx)))+((-abiez)*be123)))+((((-ae0)*bbiEz)+((atriPz*bs)-((atriPx*bbiEy)-(atriPy*bbiEx))))+(((avz*be0123)+((avx*bbiey)-(avy*bbiex)))+(ae123*bbiez)))) e021 + "
		"((((((avx*bvx)+(avy*bvy))+(avz*bvz))+0)+(0+((-ae123)*be123)))+(((as*bs)-(((abiEx*bbiEx)+(abiEy*bbiEy))+(abiEz*bbiEz)))+0)) id + "
		"(((((avy*bvz)-(avz*bvy))+(ae123*bvx))+((avx*be123)+0))+((((abiEx*bs)+(as*bbiEx))-((abiEy*bbiEz)-(abiEz*bbiEy)))+0)) e23 + "
		"(((((avz*bvx)-(avx*bvz))+(ae123*bvy))+((avy*be123)+0))+((((abiEy*bs)+(as*bbiEy))-((abiEz*bbiEx)-(abiEx*bbiEz)))+0)) e31 + "
		"(((((avx*bvy)-(avy*bvx))+(ae123*bvz))+((avz*be123)+0))+((((abiEz*bs)+(as*bbiEz))-((abiEx*bbiEy)-(abiEy*bbiEx)))+0)) e12 + "
		"(((((ae0*bvx)-(avx*be0))+((atriPy*bvz)-(atriPz*bvy)))+((((-avy)*btriPz)-((-avz)*btriPy))+((atriPx*be123)-(ae123*btriPx))))+((((abiex*bs)-(ae0123*bbiEx))-((abiey*bbiEz)-(abiez*bbiEy)))+(((as*bbiex)-(abiEx*be0123))-((abiEy*bbiez)-(abiEz*bbiey))))) e01 + "
		"(((((ae0*bvy)-(avy*be0))+((atriPz*bvx)-(atriPx*bvz)))+((((-avz)*btriPx)-((-avx)*btriPz))+((atriPy*be123)-(ae123*btriPy))))+((((abiey*bs)-(ae0123*bbiEy))-((abiez*bbiEx)-(abiex*bbiEz)))+(((as*bbiey)-(abiEy*be0123))-((abiEz*bbiex)-(abiEx*bbiez))))) e02 + "
		"(((((ae0*bvz)-(avz*be0))+((atriPx*bvy)-(atriPy*bvx)))+((((-avx)*btriPy)-((-avy)*btriPx))+((atriPz*be123)-(ae123*btriPz))))+((((abiez*bs)-(ae0123*bbiEz))-((abiex*bbiEy)-(abiey*bbiEx)))+(((as*bbiez)-(abiEz*be0123))-((abiEx*bbiey)-(abiEy*bbiex))))) e03 + "
		"(((0+(((-ae123)*be0)-(((atriPx*bvx)+(atriPy*bvy))+(atriPz*bvz))))+(((ae0*be123)+(((avx*btriPx)+(avy*btriPy))+(avz*btriPz)))+0))+(((ae0123*bs)+(((abiex*bbiEx)+(abiey*bbiEy))+(abiez*bbiEz)))+((as*be0123)+(((abiEx*bbiex)+(abiEy*bbiey))+(abiEz*bbiez))))) e0123",
		to_string(am * bm));
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
