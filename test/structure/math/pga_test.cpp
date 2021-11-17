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
	EXPECT_EQ("as id + (-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12", to_string(a02.rev()));
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
	EXPECT_EQ("(-abiex) e01 + (-abiey) e02 + (-abiez) e03 + ae0123 e0123", to_string(a24.rev()));
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
	EXPECT_EQ("(-ae123) e123 + (-atriPx) e032 + (-atriPy) e013 + (-atriPz) e021", to_string(a3.rev()));
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
		to_string(am.rev()));
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
	
	EXPECT_EQ("(as*be0) e0 + (as*bvx) e1 + (as*bvy) e2 + (as*bvz) e3", to_string(as * b1));
	EXPECT_EQ("(as*bs) id + (as*bbiEx) e23 + (as*bbiEy) e31 + (as*bbiEz) e12", to_string(as * b02));
	EXPECT_EQ("(as*bbiex) e01 + (as*bbiey) e02 + (as*bbiez) e03 + (as*be0123) e0123", to_string(as * b24));
	EXPECT_EQ("(as*be123) e123 + (as*btriPx) e032 + (as*btriPy) e013 + (as*btriPz) e021", to_string(as * b3));

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
