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

std::string vec_to_string(const std::vector<std::string>& v, std::string sep = " + ") {
	std::stringstream ss;
	bool first = true;
	for (const auto& e : v) {
		if (!first) ss << sep;
		ss << e;
		first = false;
	}
	return ss.str();
}

auto make_z() { return pga::zero<symbolic>(); }

auto make_as() { return symbolic("as"); }
auto make_a0() { return pga::blade0<symbolic>({ "as" }); }
auto make_a1() { return pga::blade1<symbolic>({ "ae0" }, { {"avx"}, {"avy"}, {"avz"} }); }
auto make_a2E() { return pga::blade2E<symbolic>({ {"abiEx"}, {"abiEy"}, {"abiEz"} }); }
auto make_a2e() { return pga::blade2e<symbolic>({ {"abiex"}, {"abiey"}, {"abiez"} }); }
auto make_a3() { return pga::blade3<symbolic>({ "ae123" }, { {"atriPx"}, {"atriPy"}, {"atriPz"} }); }
auto make_a4() { return pga::blade4<symbolic>({ "ae0123" }); }
auto make_a02E() { return pga::blade02E<symbolic>(make_a0(), make_a2E()); }
auto make_a02e() { return pga::blade02e<symbolic>(make_a0(), make_a2e()); }
auto make_a22() { return pga::blade22<symbolic>(make_a2E(), make_a2e()); }
auto make_a2E4() { return pga::blade2E4<symbolic>(make_a2E(), make_a4()); }
auto make_a2e4() { return pga::blade2e4<symbolic>(make_a2e(), make_a4()); }
auto make_a024() { return pga::blade024<symbolic>(make_a02E(), make_a2e4()); }
auto make_a13() { return pga::blade13<symbolic>(make_a1(), make_a3()); }
auto make_am() { return pga::multivector<symbolic>(make_a024(), make_a13()); }

auto make_bs() { return symbolic("bs"); }
auto make_b0() { return pga::blade0<symbolic>({ "bs" }); }
auto make_b1() { return pga::blade1<symbolic>({ "be0" }, { {"bvx"}, {"bvy"}, {"bvz"} }); }
auto make_b2E() { return pga::blade2E<symbolic>({ {"bbiEx"}, {"bbiEy"}, {"bbiEz"} }); }
auto make_b2e() { return pga::blade2e<symbolic>({ {"bbiex"}, {"bbiey"}, {"bbiez"} }); }
auto make_b3() { return pga::blade3<symbolic>({ "be123" }, { {"btriPx"}, {"btriPy"}, {"btriPz"} }); }
auto make_b4() { return pga::blade4<symbolic>({ "be0123" }); }
auto make_b02E() { return pga::blade02E<symbolic>(make_b0(), make_b2E()); }
auto make_b02e() { return pga::blade02e<symbolic>(make_b0(), make_b2e()); }
auto make_b22() { return pga::blade22<symbolic>(make_b2E(), make_b2e()); }
auto make_b2E4() { return pga::blade2E4<symbolic>(make_b2E(), make_b4()); }
auto make_b2e4() { return pga::blade2e4<symbolic>(make_b2e(), make_b4()); }
auto make_b024() { return pga::blade024<symbolic>(make_b02E(), make_b2e4()); }
auto make_b13() { return pga::blade13<symbolic>(make_b1(), make_b3()); }
auto make_bm() { return pga::multivector<symbolic>(make_b024(), make_b13()); }

} // namespace

TEST(pga_test, constructor_blade0) {
	pga::blade0<symbolic> d0;
	EXPECT_EQ("?", d0.s.v);
	pga::blade0<symbolic> a0({ "as" });
	EXPECT_EQ("as", a0.s.v);
}

TEST(pga_test, operators_arithmetic_blade0) {
	auto a0 = make_a0();
	auto b0 = make_b0();
	EXPECT_EQ("(+as) id", to_string(+a0));
	EXPECT_EQ("(-as) id", to_string(-a0));
	EXPECT_EQ("(as+bs) id", to_string(a0 + b0));
	EXPECT_EQ("(as-bs) id", to_string(a0 - b0));
	EXPECT_EQ("(as*bs) id", to_string(a0 * make_bs()));
	EXPECT_EQ("(as/bs) id", to_string(a0 / make_bs()));
	EXPECT_EQ("as id", to_string(~a0));
	EXPECT_EQ("as id", to_string(a0.rev()));
	EXPECT_EQ("as e0123", to_string(!a0));
}

TEST(pga_test, operators_inplace_blade0) {
	auto a0 = make_a0();
	auto b0 = make_b0();
	auto r = a0; r += b0;
	EXPECT_EQ("(as+bs) id", to_string(r));
	r = a0; r -= b0;
	EXPECT_EQ("(as-bs) id", to_string(r));
	r = a0; r *= make_bs();
	EXPECT_EQ("(as*bs) id", to_string(r));
	r = a0; r /= make_bs();
	EXPECT_EQ("(as/bs) id", to_string(r));
	r = a0; r += a0;
	EXPECT_EQ("(as+as) id", to_string(r));
	r = a0; r -= a0;
	EXPECT_EQ("(as-as) id", to_string(r));
}

TEST(pga_test, constructor_blade1) {
	pga::blade1<symbolic> d1;
	EXPECT_EQ("?", d1.e0.v);
	EXPECT_EQ("0", d1.v.x.v);
	EXPECT_EQ("0", d1.v.y.v);
	EXPECT_EQ("0", d1.v.z.v);
	pga::blade1<symbolic> s1({ "ae0" });
	EXPECT_EQ("ae0", s1.e0.v);
	EXPECT_EQ("0", s1.v.x.v);
	EXPECT_EQ("0", s1.v.y.v);
	EXPECT_EQ("0", s1.v.z.v);
	pga::blade1<symbolic> v1({ {"avx"}, {"avy"}, {"avz"} });
	EXPECT_EQ("0", v1.e0.v);
	EXPECT_EQ("avx", v1.v.x.v);
	EXPECT_EQ("avy", v1.v.y.v);
	EXPECT_EQ("avz", v1.v.z.v);
	pga::blade1<symbolic> a1({ "ae0" }, { {"avx"}, {"avy"}, {"avz"} });
	EXPECT_EQ("ae0", a1.e0.v);
	EXPECT_EQ("avx", a1.v.x.v);
	EXPECT_EQ("avy", a1.v.y.v);
	EXPECT_EQ("avz", a1.v.z.v);
}

TEST(pga_test, operators_arithmetic_blade1) {
	auto a1 = make_a1();
	auto b1 = make_b1();
	EXPECT_EQ("(+ae0) e0 + (+avx) e1 + (+avy) e2 + (+avz) e3", to_string(+a1));
	EXPECT_EQ("(-ae0) e0 + (-avx) e1 + (-avy) e2 + (-avz) e3", to_string(-a1));
	EXPECT_EQ("(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3", to_string(a1 + b1));
	EXPECT_EQ("(ae0-be0) e0 + (avx-bvx) e1 + (avy-bvy) e2 + (avz-bvz) e3", to_string(a1 - b1));
	EXPECT_EQ("(ae0*bs) e0 + (avx*bs) e1 + (avy*bs) e2 + (avz*bs) e3", to_string(a1 * make_bs()));
	EXPECT_EQ("(ae0/bs) e0 + (avx/bs) e1 + (avy/bs) e2 + (avz/bs) e3", to_string(a1 / make_bs()));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3", to_string(~a1));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3", to_string(a1.rev()));
	EXPECT_EQ("ae0 e123 + avx e032 + avy e013 + avz e021", to_string(!a1));
}

TEST(pga_test, operators_inplace_blade1) {
	auto a1 = make_a1();
	auto b1 = make_b1();
	auto r = a1; r += b1;
	EXPECT_EQ("(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3", to_string(r));
	r = a1; r -= b1;
	EXPECT_EQ("(ae0-be0) e0 + (avx-bvx) e1 + (avy-bvy) e2 + (avz-bvz) e3", to_string(r));
	r = a1; r *= make_bs();
	EXPECT_EQ("(ae0*bs) e0 + (avx*bs) e1 + (avy*bs) e2 + (avz*bs) e3", to_string(r));
	r = a1; r /= make_bs();
	EXPECT_EQ("(ae0/bs) e0 + (avx/bs) e1 + (avy/bs) e2 + (avz/bs) e3", to_string(r));
	r = a1; r += a1;
	EXPECT_EQ("(ae0+ae0) e0 + (avx+avx) e1 + (avy+avy) e2 + (avz+avz) e3", to_string(r));
	r = a1; r -= a1;
	EXPECT_EQ("(ae0-ae0) e0 + (avx-avx) e1 + (avy-avy) e2 + (avz-avz) e3", to_string(r));
}

TEST(pga_test, constructor_blade2E) {
	pga::blade2E<symbolic> d2E;
	EXPECT_EQ("0", d2E.biE.x.v);
	EXPECT_EQ("0", d2E.biE.y.v);
	EXPECT_EQ("0", d2E.biE.z.v);
	pga::blade2E<symbolic> a2E({ {"abiEx"}, {"abiEy"}, {"abiEz"} });
	EXPECT_EQ("abiEx", a2E.biE.x.v);
	EXPECT_EQ("abiEy", a2E.biE.y.v);
	EXPECT_EQ("abiEz", a2E.biE.z.v);
}

TEST(pga_test, operators_arithmetic_blade2E) {
	auto a2E = make_a2E();
	auto b2E = make_b2E();
	EXPECT_EQ("(+abiEx) e23 + (+abiEy) e31 + (+abiEz) e12", to_string(+a2E));
	EXPECT_EQ("(-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12", to_string(-a2E));
	EXPECT_EQ("(abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12", to_string(a2E + b2E));
	EXPECT_EQ("(abiEx-bbiEx) e23 + (abiEy-bbiEy) e31 + (abiEz-bbiEz) e12", to_string(a2E - b2E));
	EXPECT_EQ("(abiEx*bs) e23 + (abiEy*bs) e31 + (abiEz*bs) e12", to_string(a2E * make_bs()));
	EXPECT_EQ("(abiEx/bs) e23 + (abiEy/bs) e31 + (abiEz/bs) e12", to_string(a2E / make_bs()));
	EXPECT_EQ("(-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12", to_string(~a2E));
	EXPECT_EQ("(-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12", to_string(a2E.rev()));
	EXPECT_EQ("abiEx e01 + abiEy e02 + abiEz e03", to_string(!a2E));
}

TEST(pga_test, operators_inplace_blade2E) {
	auto a2E = make_a2E();
	auto b2E = make_b2E();
	auto r = a2E; r += b2E;
	EXPECT_EQ("(abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12", to_string(r));
	r = a2E; r -= b2E;
	EXPECT_EQ("(abiEx-bbiEx) e23 + (abiEy-bbiEy) e31 + (abiEz-bbiEz) e12", to_string(r));
	r = a2E; r *= make_bs();
	EXPECT_EQ("(abiEx*bs) e23 + (abiEy*bs) e31 + (abiEz*bs) e12", to_string(r));
	r = a2E; r /= make_bs();
	EXPECT_EQ("(abiEx/bs) e23 + (abiEy/bs) e31 + (abiEz/bs) e12", to_string(r));
	r = a2E; r += a2E;
	EXPECT_EQ("(abiEx+abiEx) e23 + (abiEy+abiEy) e31 + (abiEz+abiEz) e12", to_string(r));
	r = a2E; r -= a2E;
	EXPECT_EQ("(abiEx-abiEx) e23 + (abiEy-abiEy) e31 + (abiEz-abiEz) e12", to_string(r));
}

TEST(pga_test, constructor_blade2e) {
	pga::blade2e<symbolic> d2e;
	EXPECT_EQ("0", d2e.bie.x.v);
	EXPECT_EQ("0", d2e.bie.y.v);
	EXPECT_EQ("0", d2e.bie.z.v);
	pga::blade2e<symbolic> a2e({ {"abiex"}, {"abiey"}, {"abiez"} });
	EXPECT_EQ("abiex", a2e.bie.x.v);
	EXPECT_EQ("abiey", a2e.bie.y.v);
	EXPECT_EQ("abiez", a2e.bie.z.v);
}

TEST(pga_test, operators_arithmetic_blade2e) {
	auto a2e = make_a2e();
	auto b2e = make_b2e();
	EXPECT_EQ("(+abiex) e01 + (+abiey) e02 + (+abiez) e03", to_string(+a2e));
	EXPECT_EQ("(-abiex) e01 + (-abiey) e02 + (-abiez) e03", to_string(-a2e));
	EXPECT_EQ("(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03", to_string(a2e + b2e));
	EXPECT_EQ("(abiex-bbiex) e01 + (abiey-bbiey) e02 + (abiez-bbiez) e03", to_string(a2e - b2e));
	EXPECT_EQ("(abiex*bs) e01 + (abiey*bs) e02 + (abiez*bs) e03", to_string(a2e * make_bs()));
	EXPECT_EQ("(abiex/bs) e01 + (abiey/bs) e02 + (abiez/bs) e03", to_string(a2e / make_bs()));
	EXPECT_EQ("(-abiex) e01 + (-abiey) e02 + (-abiez) e03", to_string(~a2e));
	EXPECT_EQ("(-abiex) e01 + (-abiey) e02 + (-abiez) e03", to_string(a2e.rev()));
	EXPECT_EQ("abiex e23 + abiey e31 + abiez e12", to_string(!a2e));
}

TEST(pga_test, operators_inplace_blade2e) {
	auto a2e = make_a2e();
	auto b2e = make_b2e();
	auto r = a2e; r += b2e;
	EXPECT_EQ("(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03", to_string(r));
	r = a2e; r -= b2e;
	EXPECT_EQ("(abiex-bbiex) e01 + (abiey-bbiey) e02 + (abiez-bbiez) e03", to_string(r));
	r = a2e; r *= make_bs();
	EXPECT_EQ("(abiex*bs) e01 + (abiey*bs) e02 + (abiez*bs) e03", to_string(r));
	r = a2e; r /= make_bs();
	EXPECT_EQ("(abiex/bs) e01 + (abiey/bs) e02 + (abiez/bs) e03", to_string(r));
	r = a2e; r += a2e;
	EXPECT_EQ("(abiex+abiex) e01 + (abiey+abiey) e02 + (abiez+abiez) e03", to_string(r));
	r = a2e; r -= a2e;
	EXPECT_EQ("(abiex-abiex) e01 + (abiey-abiey) e02 + (abiez-abiez) e03", to_string(r));
}

TEST(pga_test, constructor_blade3) {
	pga::blade3<symbolic> d3;
	EXPECT_EQ("?", d3.e123.v);
	EXPECT_EQ("0", d3.triP.x.v);
	EXPECT_EQ("0", d3.triP.y.v);
	EXPECT_EQ("0", d3.triP.z.v);
	pga::blade3<symbolic> s3({ "ae123" });
	EXPECT_EQ("ae123", s3.e123.v);
	EXPECT_EQ("0", s3.triP.x.v);
	EXPECT_EQ("0", s3.triP.y.v);
	EXPECT_EQ("0", s3.triP.z.v);
	pga::blade3<symbolic> v3({ {"atriPx"}, {"atriPy"}, {"atriPz"} });
	EXPECT_EQ("0", v3.e123.v);
	EXPECT_EQ("atriPx", v3.triP.x.v);
	EXPECT_EQ("atriPy", v3.triP.y.v);
	EXPECT_EQ("atriPz", v3.triP.z.v);
	pga::blade3<symbolic> a3({ "ae123" }, { {"atriPx"}, {"atriPy"}, {"atriPz"} });
	EXPECT_EQ("ae123", a3.e123.v);
	EXPECT_EQ("atriPx", a3.triP.x.v);
	EXPECT_EQ("atriPy", a3.triP.y.v);
	EXPECT_EQ("atriPz", a3.triP.z.v);
}

TEST(pga_test, operators_arithmetic_blade3) {
	auto a3 = make_a3();
	auto b3 = make_b3();
	EXPECT_EQ("(+ae123) e123 + (+atriPx) e032 + (+atriPy) e013 + (+atriPz) e021", to_string(+a3));
	EXPECT_EQ("(-ae123) e123 + (-atriPx) e032 + (-atriPy) e013 + (-atriPz) e021", to_string(-a3));
	EXPECT_EQ("(ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021", to_string(a3 + b3));
	EXPECT_EQ("(ae123-be123) e123 + (atriPx-btriPx) e032 + (atriPy-btriPy) e013 + (atriPz-btriPz) e021", to_string(a3 - b3));
	EXPECT_EQ("(ae123*bs) e123 + (atriPx*bs) e032 + (atriPy*bs) e013 + (atriPz*bs) e021", to_string(a3 * make_bs()));
	EXPECT_EQ("(ae123/bs) e123 + (atriPx/bs) e032 + (atriPy/bs) e013 + (atriPz/bs) e021", to_string(a3 / make_bs()));
	EXPECT_EQ("(-ae123) e123 + (-atriPx) e032 + (-atriPy) e013 + (-atriPz) e021", to_string(~a3));
	EXPECT_EQ("(-ae123) e123 + (-atriPx) e032 + (-atriPy) e013 + (-atriPz) e021", to_string(a3.rev()));
	EXPECT_EQ("ae123 e0 + atriPx e1 + atriPy e2 + atriPz e3", to_string(!a3));
}

TEST(pga_test, operators_inplace_blade3) {
	auto a3 = make_a3();
	auto b3 = make_b3();
	auto r = a3; r += b3;
	EXPECT_EQ("(ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021", to_string(r));
	r = a3; r -= b3;
	EXPECT_EQ("(ae123-be123) e123 + (atriPx-btriPx) e032 + (atriPy-btriPy) e013 + (atriPz-btriPz) e021", to_string(r));
	r = a3; r *= make_bs();
	EXPECT_EQ("(ae123*bs) e123 + (atriPx*bs) e032 + (atriPy*bs) e013 + (atriPz*bs) e021", to_string(r));
	r = a3; r /= make_bs();
	EXPECT_EQ("(ae123/bs) e123 + (atriPx/bs) e032 + (atriPy/bs) e013 + (atriPz/bs) e021", to_string(r));
	r = a3; r += a3;
	EXPECT_EQ("(ae123+ae123) e123 + (atriPx+atriPx) e032 + (atriPy+atriPy) e013 + (atriPz+atriPz) e021", to_string(r));
	r = a3; r -= a3;
	EXPECT_EQ("(ae123-ae123) e123 + (atriPx-atriPx) e032 + (atriPy-atriPy) e013 + (atriPz-atriPz) e021", to_string(r));
}

TEST(pga_test, constructor_blade4) {
	pga::blade4<symbolic> d4;
	EXPECT_EQ("?", d4.e0123.v);
	pga::blade4<symbolic> a4({ "ae0123" });
	EXPECT_EQ("ae0123", a4.e0123.v);
}

TEST(pga_test, operators_arithmetic_blade4) {
	auto a4 = make_a4();
	auto b4 = make_b4();
	EXPECT_EQ("(+ae0123) e0123", to_string(+a4));
	EXPECT_EQ("(-ae0123) e0123", to_string(-a4));
	EXPECT_EQ("(ae0123+be0123) e0123", to_string(a4 + b4));
	EXPECT_EQ("(ae0123-be0123) e0123", to_string(a4 - b4));
	EXPECT_EQ("(ae0123*bs) e0123", to_string(a4 * make_bs()));
	EXPECT_EQ("(ae0123/bs) e0123", to_string(a4 / make_bs()));
	EXPECT_EQ("ae0123 e0123", to_string(~a4));
	EXPECT_EQ("ae0123 e0123", to_string(a4.rev()));
	EXPECT_EQ("ae0123 id", to_string(!a4));
}

TEST(pga_test, operators_inplace_blade4) {
	auto a4 = make_a4();
	auto b4 = make_b4();
	auto r = a4; r += b4;
	EXPECT_EQ("(ae0123+be0123) e0123", to_string(r));
	r = a4; r -= b4;
	EXPECT_EQ("(ae0123-be0123) e0123", to_string(r));
	r = a4; r *= make_bs();
	EXPECT_EQ("(ae0123*bs) e0123", to_string(r));
	r = a4; r /= make_bs();
	EXPECT_EQ("(ae0123/bs) e0123", to_string(r));
	r = a4; r += a4;
	EXPECT_EQ("(ae0123+ae0123) e0123", to_string(r));
	r = a4; r -= a4;
	EXPECT_EQ("(ae0123-ae0123) e0123", to_string(r));
}

TEST(pga_test, constructor_blade02E) {
	pga::blade02E<symbolic> d02;
	EXPECT_EQ("?", d02.b0.s.v);
	EXPECT_EQ("0", d02.b2E.biE.x.v);
	EXPECT_EQ("0", d02.b2E.biE.y.v);
	EXPECT_EQ("0", d02.b2E.biE.z.v);
	pga::blade02E<symbolic> s02(make_a0());
	EXPECT_EQ("as", s02.b0.s.v);
	EXPECT_EQ("0", s02.b2E.biE.x.v);
	EXPECT_EQ("0", s02.b2E.biE.y.v);
	EXPECT_EQ("0", s02.b2E.biE.z.v);
	pga::blade02E<symbolic> v02(make_a2E());
	EXPECT_EQ("0", v02.b0.s.v);
	EXPECT_EQ("abiEx", v02.b2E.biE.x.v);
	EXPECT_EQ("abiEy", v02.b2E.biE.y.v);
	EXPECT_EQ("abiEz", v02.b2E.biE.z.v);
	pga::blade02E<symbolic> a02(make_a0(), make_a2E());
	EXPECT_EQ("as", a02.b0.s.v);
	EXPECT_EQ("abiEx", a02.b2E.biE.x.v);
	EXPECT_EQ("abiEy", a02.b2E.biE.y.v);
	EXPECT_EQ("abiEz", a02.b2E.biE.z.v);
}

TEST(pga_test, operators_arithmetic_blade02E) {
	auto a02 = make_a02E();
	auto b02 = make_b02E();
	EXPECT_EQ("(+as) id + (+abiEx) e23 + (+abiEy) e31 + (+abiEz) e12", to_string(+a02));
	EXPECT_EQ("(-as) id + (-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12", to_string(-a02));
	EXPECT_EQ("(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12", to_string(a02 + b02));
	EXPECT_EQ("(as-bs) id + (abiEx-bbiEx) e23 + (abiEy-bbiEy) e31 + (abiEz-bbiEz) e12", to_string(a02 - b02));
	EXPECT_EQ("(as*bs) id + (abiEx*bs) e23 + (abiEy*bs) e31 + (abiEz*bs) e12", to_string(a02 * make_bs()));
	EXPECT_EQ("(as/bs) id + (abiEx/bs) e23 + (abiEy/bs) e31 + (abiEz/bs) e12", to_string(a02 / make_bs()));
	EXPECT_EQ("as id + (-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12", to_string(~a02));
	EXPECT_EQ("as id + (-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12", to_string(a02.rev()));
	EXPECT_EQ("abiEx e01 + abiEy e02 + abiEz e03 + as e0123", to_string(!a02));
}

TEST(pga_test, operators_inplace_blade02E) {
	auto a02 = make_a02E();
	auto b02 = make_b02E();
	auto r = a02; r += b02;
	EXPECT_EQ("(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12", to_string(r));
	r = a02; r -= b02;
	EXPECT_EQ("(as-bs) id + (abiEx-bbiEx) e23 + (abiEy-bbiEy) e31 + (abiEz-bbiEz) e12", to_string(r));
	r = a02; r *= make_bs();
	EXPECT_EQ("(as*bs) id + (abiEx*bs) e23 + (abiEy*bs) e31 + (abiEz*bs) e12", to_string(r));
	r = a02; r /= make_bs();
	EXPECT_EQ("(as/bs) id + (abiEx/bs) e23 + (abiEy/bs) e31 + (abiEz/bs) e12", to_string(r));
	r = a02; r += a02;
	EXPECT_EQ("(as+as) id + (abiEx+abiEx) e23 + (abiEy+abiEy) e31 + (abiEz+abiEz) e12", to_string(r));
	r = a02; r -= a02;
	EXPECT_EQ("(as-as) id + (abiEx-abiEx) e23 + (abiEy-abiEy) e31 + (abiEz-abiEz) e12", to_string(r));
}

TEST(pga_test, constructor_blade02e) {
	pga::blade02e<symbolic> d02;
	EXPECT_EQ("?", d02.b0.s.v);
	EXPECT_EQ("0", d02.b2e.bie.x.v);
	EXPECT_EQ("0", d02.b2e.bie.y.v);
	EXPECT_EQ("0", d02.b2e.bie.z.v);
	pga::blade02e<symbolic> s02(make_a0());
	EXPECT_EQ("as", s02.b0.s.v);
	EXPECT_EQ("0", s02.b2e.bie.x.v);
	EXPECT_EQ("0", s02.b2e.bie.y.v);
	EXPECT_EQ("0", s02.b2e.bie.z.v);
	pga::blade02e<symbolic> v02(make_a2e());
	EXPECT_EQ("0", v02.b0.s.v);
	EXPECT_EQ("abiex", v02.b2e.bie.x.v);
	EXPECT_EQ("abiey", v02.b2e.bie.y.v);
	EXPECT_EQ("abiez", v02.b2e.bie.z.v);
	pga::blade02e<symbolic> a02(make_a0(), make_a2e());
	EXPECT_EQ("as", a02.b0.s.v);
	EXPECT_EQ("abiex", a02.b2e.bie.x.v);
	EXPECT_EQ("abiey", a02.b2e.bie.y.v);
	EXPECT_EQ("abiez", a02.b2e.bie.z.v);
}

TEST(pga_test, operators_arithmetic_blade02e) {
	auto a02 = make_a02e();
	auto b02 = make_b02e();
	EXPECT_EQ("(+as) id + (+abiex) e01 + (+abiey) e02 + (+abiez) e03", to_string(+a02));
	EXPECT_EQ("(-as) id + (-abiex) e01 + (-abiey) e02 + (-abiez) e03", to_string(-a02));
	EXPECT_EQ("(as+bs) id + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03", to_string(a02 + b02));
	EXPECT_EQ("(as-bs) id + (abiex-bbiex) e01 + (abiey-bbiey) e02 + (abiez-bbiez) e03", to_string(a02 - b02));
	EXPECT_EQ("(as*bs) id + (abiex*bs) e01 + (abiey*bs) e02 + (abiez*bs) e03", to_string(a02 * make_bs()));
	EXPECT_EQ("(as/bs) id + (abiex/bs) e01 + (abiey/bs) e02 + (abiez/bs) e03", to_string(a02 / make_bs()));
	EXPECT_EQ("as id + (-abiex) e01 + (-abiey) e02 + (-abiez) e03", to_string(~a02));
	EXPECT_EQ("as id + (-abiex) e01 + (-abiey) e02 + (-abiez) e03", to_string(a02.rev()));
	EXPECT_EQ("abiex e23 + abiey e31 + abiez e12 + as e0123", to_string(!a02));
}

TEST(pga_test, operators_inplace_blade02e) {
	auto a02 = make_a02e();
	auto b02 = make_b02e();
	auto r = a02; r += b02;
	EXPECT_EQ("(as+bs) id + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03", to_string(r));
	r = a02; r -= b02;
	EXPECT_EQ("(as-bs) id + (abiex-bbiex) e01 + (abiey-bbiey) e02 + (abiez-bbiez) e03", to_string(r));
	r = a02; r *= make_bs();
	EXPECT_EQ("(as*bs) id + (abiex*bs) e01 + (abiey*bs) e02 + (abiez*bs) e03", to_string(r));
	r = a02; r /= make_bs();
	EXPECT_EQ("(as/bs) id + (abiex/bs) e01 + (abiey/bs) e02 + (abiez/bs) e03", to_string(r));
	r = a02; r += a02;
	EXPECT_EQ("(as+as) id + (abiex+abiex) e01 + (abiey+abiey) e02 + (abiez+abiez) e03", to_string(r));
	r = a02; r -= a02;
	EXPECT_EQ("(as-as) id + (abiex-abiex) e01 + (abiey-abiey) e02 + (abiez-abiez) e03", to_string(r));
}

TEST(pga_test, constructor_blade22) {
	pga::blade22<symbolic> d22;
	EXPECT_EQ("0", d22.b2E.biE.x.v);
	EXPECT_EQ("0", d22.b2E.biE.y.v);
	EXPECT_EQ("0", d22.b2E.biE.z.v);
	EXPECT_EQ("0", d22.b2e.bie.x.v);
	EXPECT_EQ("0", d22.b2e.bie.y.v);
	EXPECT_EQ("0", d22.b2e.bie.z.v);
	pga::blade22<symbolic> E22(make_a2E());
	EXPECT_EQ("abiEx", E22.b2E.biE.x.v);
	EXPECT_EQ("abiEy", E22.b2E.biE.y.v);
	EXPECT_EQ("abiEz", E22.b2E.biE.z.v);
	EXPECT_EQ("0", E22.b2e.bie.x.v);
	EXPECT_EQ("0", E22.b2e.bie.y.v);
	EXPECT_EQ("0", E22.b2e.bie.z.v);
	pga::blade22<symbolic> e22(make_a2e());
	EXPECT_EQ("0", e22.b2E.biE.x.v);
	EXPECT_EQ("0", e22.b2E.biE.y.v);
	EXPECT_EQ("0", e22.b2E.biE.z.v);
	EXPECT_EQ("abiex", e22.b2e.bie.x.v);
	EXPECT_EQ("abiey", e22.b2e.bie.y.v);
	EXPECT_EQ("abiez", e22.b2e.bie.z.v);
	pga::blade22<symbolic> a22(make_a2E(), make_a2e());
	EXPECT_EQ("abiEx", a22.b2E.biE.x.v);
	EXPECT_EQ("abiEy", a22.b2E.biE.y.v);
	EXPECT_EQ("abiEz", a22.b2E.biE.z.v);
	EXPECT_EQ("abiex", a22.b2e.bie.x.v);
	EXPECT_EQ("abiey", a22.b2e.bie.y.v);
	EXPECT_EQ("abiez", a22.b2e.bie.z.v);
}

TEST(pga_test, operators_arithmetic_blade22) {
	auto a22 = make_a22();
	auto b22 = make_b22();
	EXPECT_EQ("(+abiEx) e23 + (+abiEy) e31 + (+abiEz) e12 + (+abiex) e01 + (+abiey) e02 + (+abiez) e03", to_string(+a22));
	EXPECT_EQ("(-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12 + (-abiex) e01 + (-abiey) e02 + (-abiez) e03", to_string(-a22));
	EXPECT_EQ("(abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03", to_string(a22 + b22));
	EXPECT_EQ("(abiEx-bbiEx) e23 + (abiEy-bbiEy) e31 + (abiEz-bbiEz) e12 + (abiex-bbiex) e01 + (abiey-bbiey) e02 + (abiez-bbiez) e03", to_string(a22 - b22));
	EXPECT_EQ("(abiEx*bs) e23 + (abiEy*bs) e31 + (abiEz*bs) e12 + (abiex*bs) e01 + (abiey*bs) e02 + (abiez*bs) e03", to_string(a22 * make_bs()));
	EXPECT_EQ("(abiEx/bs) e23 + (abiEy/bs) e31 + (abiEz/bs) e12 + (abiex/bs) e01 + (abiey/bs) e02 + (abiez/bs) e03", to_string(a22 / make_bs()));
	EXPECT_EQ("(-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12 + (-abiex) e01 + (-abiey) e02 + (-abiez) e03", to_string(~a22));
	EXPECT_EQ("(-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12 + (-abiex) e01 + (-abiey) e02 + (-abiez) e03", to_string(a22.rev()));
	EXPECT_EQ("abiex e23 + abiey e31 + abiez e12 + abiEx e01 + abiEy e02 + abiEz e03", to_string(!a22));
}

TEST(pga_test, operators_inplace_blade22) {
	auto a22 = make_a22();
	auto b22 = make_b22();
	auto r = a22; r += b22;
	EXPECT_EQ("(abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03", to_string(r));
	r = a22; r -= b22;
	EXPECT_EQ("(abiEx-bbiEx) e23 + (abiEy-bbiEy) e31 + (abiEz-bbiEz) e12 + (abiex-bbiex) e01 + (abiey-bbiey) e02 + (abiez-bbiez) e03", to_string(r));
	r = a22; r *= make_bs();
	EXPECT_EQ("(abiEx*bs) e23 + (abiEy*bs) e31 + (abiEz*bs) e12 + (abiex*bs) e01 + (abiey*bs) e02 + (abiez*bs) e03", to_string(r));
	r = a22; r /= make_bs();
	EXPECT_EQ("(abiEx/bs) e23 + (abiEy/bs) e31 + (abiEz/bs) e12 + (abiex/bs) e01 + (abiey/bs) e02 + (abiez/bs) e03", to_string(r));
	r = a22; r += a22;
	EXPECT_EQ("(abiEx+abiEx) e23 + (abiEy+abiEy) e31 + (abiEz+abiEz) e12 + (abiex+abiex) e01 + (abiey+abiey) e02 + (abiez+abiez) e03", to_string(r));
	r = a22; r -= a22;
	EXPECT_EQ("(abiEx-abiEx) e23 + (abiEy-abiEy) e31 + (abiEz-abiEz) e12 + (abiex-abiex) e01 + (abiey-abiey) e02 + (abiez-abiez) e03", to_string(r));
}

TEST(pga_test, constructor_blade2E4) {
	pga::blade2E4<symbolic> d24;
	EXPECT_EQ("0", d24.b2E.biE.x.v);
	EXPECT_EQ("0", d24.b2E.biE.y.v);
	EXPECT_EQ("0", d24.b2E.biE.z.v);
	EXPECT_EQ("?", d24.b4.e0123.v);
	pga::blade2E4<symbolic> s24(make_a4());
	EXPECT_EQ("0", s24.b2E.biE.x.v);
	EXPECT_EQ("0", s24.b2E.biE.y.v);
	EXPECT_EQ("0", s24.b2E.biE.z.v);
	EXPECT_EQ("ae0123", s24.b4.e0123.v);
	pga::blade2E4<symbolic> v24(make_a2E());
	EXPECT_EQ("abiEx", v24.b2E.biE.x.v);
	EXPECT_EQ("abiEy", v24.b2E.biE.y.v);
	EXPECT_EQ("abiEz", v24.b2E.biE.z.v);
	EXPECT_EQ("0", v24.b4.e0123.v);
	pga::blade2E4<symbolic> a24(make_a2E(), make_a4());
	EXPECT_EQ("abiEx", a24.b2E.biE.x.v);
	EXPECT_EQ("abiEy", a24.b2E.biE.y.v);
	EXPECT_EQ("abiEz", a24.b2E.biE.z.v);
	EXPECT_EQ("ae0123", a24.b4.e0123.v);
}

TEST(pga_test, operators_arithmetic_blade2E4) {
	auto a24 = make_a2E4();
	auto b24 = make_b2E4();
	EXPECT_EQ("(+abiEx) e23 + (+abiEy) e31 + (+abiEz) e12 + (+ae0123) e0123", to_string(+a24));
	EXPECT_EQ("(-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12 + (-ae0123) e0123", to_string(-a24));
	EXPECT_EQ("(abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + (ae0123+be0123) e0123", to_string(a24 + b24));
	EXPECT_EQ("(abiEx-bbiEx) e23 + (abiEy-bbiEy) e31 + (abiEz-bbiEz) e12 + (ae0123-be0123) e0123", to_string(a24 - b24));
	EXPECT_EQ("(abiEx*bs) e23 + (abiEy*bs) e31 + (abiEz*bs) e12 + (ae0123*bs) e0123", to_string(a24 * make_bs()));
	EXPECT_EQ("(abiEx/bs) e23 + (abiEy/bs) e31 + (abiEz/bs) e12 + (ae0123/bs) e0123", to_string(a24 / make_bs()));
	EXPECT_EQ("(-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12 + ae0123 e0123", to_string(~a24));
	EXPECT_EQ("(-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12 + ae0123 e0123", to_string(a24.rev()));
	EXPECT_EQ("ae0123 id + abiEx e01 + abiEy e02 + abiEz e03", to_string(!a24));
}

TEST(pga_test, operators_inplace_blade2E4) {
	auto a24 = make_a2E4();
	auto b24 = make_b2E4();
	auto r = a24; r += b24;
	EXPECT_EQ("(abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + (ae0123+be0123) e0123", to_string(r));
	r = a24; r -= b24;
	EXPECT_EQ("(abiEx-bbiEx) e23 + (abiEy-bbiEy) e31 + (abiEz-bbiEz) e12 + (ae0123-be0123) e0123", to_string(r));
	r = a24; r *= make_bs();
	EXPECT_EQ("(abiEx*bs) e23 + (abiEy*bs) e31 + (abiEz*bs) e12 + (ae0123*bs) e0123", to_string(r));
	r = a24; r /= make_bs();
	EXPECT_EQ("(abiEx/bs) e23 + (abiEy/bs) e31 + (abiEz/bs) e12 + (ae0123/bs) e0123", to_string(r));
	r = a24; r += a24;
	EXPECT_EQ("(abiEx+abiEx) e23 + (abiEy+abiEy) e31 + (abiEz+abiEz) e12 + (ae0123+ae0123) e0123", to_string(r));
	r = a24; r -= a24;
	EXPECT_EQ("(abiEx-abiEx) e23 + (abiEy-abiEy) e31 + (abiEz-abiEz) e12 + (ae0123-ae0123) e0123", to_string(r));
}

TEST(pga_test, constructor_blade2e4) {
	pga::blade2e4<symbolic> d24;
	EXPECT_EQ("0", d24.b2e.bie.x.v);
	EXPECT_EQ("0", d24.b2e.bie.y.v);
	EXPECT_EQ("0", d24.b2e.bie.z.v);
	EXPECT_EQ("?", d24.b4.e0123.v);
	pga::blade2e4<symbolic> s24(make_a4());
	EXPECT_EQ("0", s24.b2e.bie.x.v);
	EXPECT_EQ("0", s24.b2e.bie.y.v);
	EXPECT_EQ("0", s24.b2e.bie.z.v);
	EXPECT_EQ("ae0123", s24.b4.e0123.v);
	pga::blade2e4<symbolic> v24(make_a2e());
	EXPECT_EQ("abiex", v24.b2e.bie.x.v);
	EXPECT_EQ("abiey", v24.b2e.bie.y.v);
	EXPECT_EQ("abiez", v24.b2e.bie.z.v);
	EXPECT_EQ("0", v24.b4.e0123.v);
	pga::blade2e4<symbolic> a24(make_a2e(), make_a4());
	EXPECT_EQ("abiex", a24.b2e.bie.x.v);
	EXPECT_EQ("abiey", a24.b2e.bie.y.v);
	EXPECT_EQ("abiez", a24.b2e.bie.z.v);
	EXPECT_EQ("ae0123", a24.b4.e0123.v);
}

TEST(pga_test, operators_arithmetic_blade2e4) {
	auto a24 = make_a2e4();
	auto b24 = make_b2e4();
	EXPECT_EQ("(+abiex) e01 + (+abiey) e02 + (+abiez) e03 + (+ae0123) e0123", to_string(+a24));
	EXPECT_EQ("(-abiex) e01 + (-abiey) e02 + (-abiez) e03 + (-ae0123) e0123", to_string(-a24));
	EXPECT_EQ("(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123", to_string(a24 + b24));
	EXPECT_EQ("(abiex-bbiex) e01 + (abiey-bbiey) e02 + (abiez-bbiez) e03 + (ae0123-be0123) e0123", to_string(a24 - b24));
	EXPECT_EQ("(abiex*bs) e01 + (abiey*bs) e02 + (abiez*bs) e03 + (ae0123*bs) e0123", to_string(a24 * make_bs()));
	EXPECT_EQ("(abiex/bs) e01 + (abiey/bs) e02 + (abiez/bs) e03 + (ae0123/bs) e0123", to_string(a24 / make_bs()));
	EXPECT_EQ("(-abiex) e01 + (-abiey) e02 + (-abiez) e03 + ae0123 e0123", to_string(~a24));
	EXPECT_EQ("(-abiex) e01 + (-abiey) e02 + (-abiez) e03 + ae0123 e0123", to_string(a24.rev()));
	EXPECT_EQ("ae0123 id + abiex e23 + abiey e31 + abiez e12", to_string(!a24));
}

TEST(pga_test, operators_inplace_blade2e4) {
	auto a24 = make_a2e4();
	auto b24 = make_b2e4();
	auto r = a24; r += b24;
	EXPECT_EQ("(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123", to_string(r));
	r = a24; r -= b24;
	EXPECT_EQ("(abiex-bbiex) e01 + (abiey-bbiey) e02 + (abiez-bbiez) e03 + (ae0123-be0123) e0123", to_string(r));
	r = a24; r *= make_bs();
	EXPECT_EQ("(abiex*bs) e01 + (abiey*bs) e02 + (abiez*bs) e03 + (ae0123*bs) e0123", to_string(r));
	r = a24; r /= make_bs();
	EXPECT_EQ("(abiex/bs) e01 + (abiey/bs) e02 + (abiez/bs) e03 + (ae0123/bs) e0123", to_string(r));
	r = a24; r += a24;
	EXPECT_EQ("(abiex+abiex) e01 + (abiey+abiey) e02 + (abiez+abiez) e03 + (ae0123+ae0123) e0123", to_string(r));
	r = a24; r -= a24;
	EXPECT_EQ("(abiex-abiex) e01 + (abiey-abiey) e02 + (abiez-abiez) e03 + (ae0123-ae0123) e0123", to_string(r));
}

TEST(pga_test, constructor_blade024) {
	pga::blade024<symbolic> d024;
	EXPECT_EQ("?", d024.b02.b0.s.v);
	EXPECT_EQ("0", d024.b02.b2E.biE.x.v);
	EXPECT_EQ("0", d024.b02.b2E.biE.y.v);
	EXPECT_EQ("0", d024.b02.b2E.biE.z.v);
	EXPECT_EQ("0", d024.b24.b2e.bie.x.v);
	EXPECT_EQ("0", d024.b24.b2e.bie.y.v);
	EXPECT_EQ("0", d024.b24.b2e.bie.z.v);
	EXPECT_EQ("?", d024.b24.b4.e0123.v);
	pga::blade024<symbolic> v0(make_a0());
	EXPECT_EQ("as", v0.b02.b0.s.v);
	EXPECT_EQ("0", v0.b02.b2E.biE.x.v);
	EXPECT_EQ("0", v0.b02.b2E.biE.y.v);
	EXPECT_EQ("0", v0.b02.b2E.biE.z.v);
	EXPECT_EQ("0", v0.b24.b2e.bie.x.v);
	EXPECT_EQ("0", v0.b24.b2e.bie.y.v);
	EXPECT_EQ("0", v0.b24.b2e.bie.z.v);
	EXPECT_EQ("0", v0.b24.b4.e0123.v);
	pga::blade024<symbolic> v2E(make_a2E());
	EXPECT_EQ("0", v2E.b02.b0.s.v);
	EXPECT_EQ("abiEx", v2E.b02.b2E.biE.x.v);
	EXPECT_EQ("abiEy", v2E.b02.b2E.biE.y.v);
	EXPECT_EQ("abiEz", v2E.b02.b2E.biE.z.v);
	EXPECT_EQ("0", v2E.b24.b2e.bie.x.v);
	EXPECT_EQ("0", v2E.b24.b2e.bie.y.v);
	EXPECT_EQ("0", v2E.b24.b2e.bie.z.v);
	EXPECT_EQ("0", v2E.b24.b4.e0123.v);
	pga::blade024<symbolic> v2e(make_a2e());
	EXPECT_EQ("0", v2e.b02.b0.s.v);
	EXPECT_EQ("0", v2e.b02.b2E.biE.x.v);
	EXPECT_EQ("0", v2e.b02.b2E.biE.y.v);
	EXPECT_EQ("0", v2e.b02.b2E.biE.z.v);
	EXPECT_EQ("abiex", v2e.b24.b2e.bie.x.v);
	EXPECT_EQ("abiey", v2e.b24.b2e.bie.y.v);
	EXPECT_EQ("abiez", v2e.b24.b2e.bie.z.v);
	EXPECT_EQ("0", v2e.b24.b4.e0123.v);
	pga::blade024<symbolic> v4(make_a4());
	EXPECT_EQ("0", v4.b02.b0.s.v);
	EXPECT_EQ("0", v4.b02.b2E.biE.x.v);
	EXPECT_EQ("0", v4.b02.b2E.biE.y.v);
	EXPECT_EQ("0", v4.b02.b2E.biE.z.v);
	EXPECT_EQ("0", v4.b24.b2e.bie.x.v);
	EXPECT_EQ("0", v4.b24.b2e.bie.y.v);
	EXPECT_EQ("0", v4.b24.b2e.bie.z.v);
	EXPECT_EQ("ae0123", v4.b24.b4.e0123.v);
	pga::blade024<symbolic> v02E(make_a02E());
	EXPECT_EQ("as", v02E.b02.b0.s.v);
	EXPECT_EQ("abiEx", v02E.b02.b2E.biE.x.v);
	EXPECT_EQ("abiEy", v02E.b02.b2E.biE.y.v);
	EXPECT_EQ("abiEz", v02E.b02.b2E.biE.z.v);
	EXPECT_EQ("0", v02E.b24.b2e.bie.x.v);
	EXPECT_EQ("0", v02E.b24.b2e.bie.y.v);
	EXPECT_EQ("0", v02E.b24.b2e.bie.z.v);
	EXPECT_EQ("0", v02E.b24.b4.e0123.v);
	pga::blade024<symbolic> v02e(make_a02e());
	EXPECT_EQ("as", v02e.b02.b0.s.v);
	EXPECT_EQ("0", v02e.b02.b2E.biE.x.v);
	EXPECT_EQ("0", v02e.b02.b2E.biE.y.v);
	EXPECT_EQ("0", v02e.b02.b2E.biE.z.v);
	EXPECT_EQ("abiex", v02e.b24.b2e.bie.x.v);
	EXPECT_EQ("abiey", v02e.b24.b2e.bie.y.v);
	EXPECT_EQ("abiez", v02e.b24.b2e.bie.z.v);
	EXPECT_EQ("0", v02e.b24.b4.e0123.v);
	pga::blade024<symbolic> v22(make_a22());
	EXPECT_EQ("0", v22.b02.b0.s.v);
	EXPECT_EQ("abiEx", v22.b02.b2E.biE.x.v);
	EXPECT_EQ("abiEy", v22.b02.b2E.biE.y.v);
	EXPECT_EQ("abiEz", v22.b02.b2E.biE.z.v);
	EXPECT_EQ("abiex", v22.b24.b2e.bie.x.v);
	EXPECT_EQ("abiey", v22.b24.b2e.bie.y.v);
	EXPECT_EQ("abiez", v22.b24.b2e.bie.z.v);
	EXPECT_EQ("0", v22.b24.b4.e0123.v);
	pga::blade024<symbolic> v2E4(make_a2E4());
	EXPECT_EQ("0", v2E4.b02.b0.s.v);
	EXPECT_EQ("abiEx", v2E4.b02.b2E.biE.x.v);
	EXPECT_EQ("abiEy", v2E4.b02.b2E.biE.y.v);
	EXPECT_EQ("abiEz", v2E4.b02.b2E.biE.z.v);
	EXPECT_EQ("0", v2E4.b24.b2e.bie.x.v);
	EXPECT_EQ("0", v2E4.b24.b2e.bie.y.v);
	EXPECT_EQ("0", v2E4.b24.b2e.bie.z.v);
	EXPECT_EQ("ae0123", v2E4.b24.b4.e0123.v);
	pga::blade024<symbolic> v2e4(make_a2e4());
	EXPECT_EQ("0", v2e4.b02.b0.s.v);
	EXPECT_EQ("0", v2e4.b02.b2E.biE.x.v);
	EXPECT_EQ("0", v2e4.b02.b2E.biE.y.v);
	EXPECT_EQ("0", v2e4.b02.b2E.biE.z.v);
	EXPECT_EQ("abiex", v2e4.b24.b2e.bie.x.v);
	EXPECT_EQ("abiey", v2e4.b24.b2e.bie.y.v);
	EXPECT_EQ("abiez", v2e4.b24.b2e.bie.z.v);
	EXPECT_EQ("ae0123", v2e4.b24.b4.e0123.v);
	pga::blade024<symbolic> a024(make_a02E(), make_a2e4());
	EXPECT_EQ("as", a024.b02.b0.s.v);
	EXPECT_EQ("abiEx", a024.b02.b2E.biE.x.v);
	EXPECT_EQ("abiEy", a024.b02.b2E.biE.y.v);
	EXPECT_EQ("abiEz", a024.b02.b2E.biE.z.v);
	EXPECT_EQ("abiex", a024.b24.b2e.bie.x.v);
	EXPECT_EQ("abiey", a024.b24.b2e.bie.y.v);
	EXPECT_EQ("abiez", a024.b24.b2e.bie.z.v);
	EXPECT_EQ("ae0123", a024.b24.b4.e0123.v);
}

TEST(pga_test, operators_arithmetic_blade024) {
	auto a024 = make_a024();
	auto b024 = make_b024();
	EXPECT_EQ("(+as) id + (+abiEx) e23 + (+abiEy) e31 + (+abiEz) e12 + (+abiex) e01 + (+abiey) e02 + (+abiez) e03 + (+ae0123) e0123", to_string(+a024));
	EXPECT_EQ("(-as) id + (-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12 + (-abiex) e01 + (-abiey) e02 + (-abiez) e03 + (-ae0123) e0123", to_string(-a024));
	EXPECT_EQ("(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123", to_string(a024 + b024));
	EXPECT_EQ("(as-bs) id + (abiEx-bbiEx) e23 + (abiEy-bbiEy) e31 + (abiEz-bbiEz) e12 + (abiex-bbiex) e01 + (abiey-bbiey) e02 + (abiez-bbiez) e03 + (ae0123-be0123) e0123", to_string(a024 - b024));
	EXPECT_EQ("(as*bs) id + (abiEx*bs) e23 + (abiEy*bs) e31 + (abiEz*bs) e12 + (abiex*bs) e01 + (abiey*bs) e02 + (abiez*bs) e03 + (ae0123*bs) e0123", to_string(a024 * make_bs()));
	EXPECT_EQ("(as/bs) id + (abiEx/bs) e23 + (abiEy/bs) e31 + (abiEz/bs) e12 + (abiex/bs) e01 + (abiey/bs) e02 + (abiez/bs) e03 + (ae0123/bs) e0123", to_string(a024 / make_bs()));
	EXPECT_EQ("as id + (-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12 + (-abiex) e01 + (-abiey) e02 + (-abiez) e03 + ae0123 e0123", to_string(~a024));
	EXPECT_EQ("as id + (-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12 + (-abiex) e01 + (-abiey) e02 + (-abiez) e03 + ae0123 e0123", to_string(a024.rev()));
	EXPECT_EQ("ae0123 id + abiex e23 + abiey e31 + abiez e12 + abiEx e01 + abiEy e02 + abiEz e03 + as e0123", to_string(!a024));
}

TEST(pga_test, operators_inplace_blade024) {
	auto a024 = make_a024();
	auto b024 = make_b024();
	auto r = a024; r += b024;
	EXPECT_EQ("(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123", to_string(r));
	r = a024; r -= b024;
	EXPECT_EQ("(as-bs) id + (abiEx-bbiEx) e23 + (abiEy-bbiEy) e31 + (abiEz-bbiEz) e12 + (abiex-bbiex) e01 + (abiey-bbiey) e02 + (abiez-bbiez) e03 + (ae0123-be0123) e0123", to_string(r));
	r = a024; r *= make_bs();
	EXPECT_EQ("(as*bs) id + (abiEx*bs) e23 + (abiEy*bs) e31 + (abiEz*bs) e12 + (abiex*bs) e01 + (abiey*bs) e02 + (abiez*bs) e03 + (ae0123*bs) e0123", to_string(r));
	r = a024; r /= make_bs();
	EXPECT_EQ("(as/bs) id + (abiEx/bs) e23 + (abiEy/bs) e31 + (abiEz/bs) e12 + (abiex/bs) e01 + (abiey/bs) e02 + (abiez/bs) e03 + (ae0123/bs) e0123", to_string(r));
	r = a024; r += a024;
	EXPECT_EQ("(as+as) id + (abiEx+abiEx) e23 + (abiEy+abiEy) e31 + (abiEz+abiEz) e12 + (abiex+abiex) e01 + (abiey+abiey) e02 + (abiez+abiez) e03 + (ae0123+ae0123) e0123", to_string(r));
	r = a024; r -= a024;
	EXPECT_EQ("(as-as) id + (abiEx-abiEx) e23 + (abiEy-abiEy) e31 + (abiEz-abiEz) e12 + (abiex-abiex) e01 + (abiey-abiey) e02 + (abiez-abiez) e03 + (ae0123-ae0123) e0123", to_string(r));
}

TEST(pga_test, constructor_blade13) {
	pga::blade13<symbolic> d13;
	EXPECT_EQ("?", d13.b1.e0.v);
	EXPECT_EQ("0", d13.b1.v.x.v);
	EXPECT_EQ("0", d13.b1.v.y.v);
	EXPECT_EQ("0", d13.b1.v.z.v);
	EXPECT_EQ("?", d13.b3.e123.v);
	EXPECT_EQ("0", d13.b3.triP.x.v);
	EXPECT_EQ("0", d13.b3.triP.y.v);
	EXPECT_EQ("0", d13.b3.triP.z.v);
	pga::blade13<symbolic> s13(make_a1());
	EXPECT_EQ("ae0", s13.b1.e0.v);
	EXPECT_EQ("avx", s13.b1.v.x.v);
	EXPECT_EQ("avy", s13.b1.v.y.v);
	EXPECT_EQ("avz", s13.b1.v.z.v);
	EXPECT_EQ("0", s13.b3.e123.v);
	EXPECT_EQ("0", s13.b3.triP.x.v);
	EXPECT_EQ("0", s13.b3.triP.y.v);
	EXPECT_EQ("0", s13.b3.triP.z.v);
	pga::blade13<symbolic> v13(make_a3());
	EXPECT_EQ("0", v13.b1.e0.v);
	EXPECT_EQ("0", v13.b1.v.x.v);
	EXPECT_EQ("0", v13.b1.v.y.v);
	EXPECT_EQ("0", v13.b1.v.z.v);
	EXPECT_EQ("ae123", v13.b3.e123.v);
	EXPECT_EQ("atriPx", v13.b3.triP.x.v);
	EXPECT_EQ("atriPy", v13.b3.triP.y.v);
	EXPECT_EQ("atriPz", v13.b3.triP.z.v);
	pga::blade13<symbolic> a13(make_a1(), make_a3());
	EXPECT_EQ("ae0", a13.b1.e0.v);
	EXPECT_EQ("avx", a13.b1.v.x.v);
	EXPECT_EQ("avy", a13.b1.v.y.v);
	EXPECT_EQ("avz", a13.b1.v.z.v);
	EXPECT_EQ("ae123", a13.b3.e123.v);
	EXPECT_EQ("atriPx", a13.b3.triP.x.v);
	EXPECT_EQ("atriPy", a13.b3.triP.y.v);
	EXPECT_EQ("atriPz", a13.b3.triP.z.v);
}

TEST(pga_test, operators_arithmetic_blade13) {
	auto a13 = make_a13();
	auto b13 = make_b13();
	EXPECT_EQ("(+ae0) e0 + (+avx) e1 + (+avy) e2 + (+avz) e3 + (+ae123) e123 + (+atriPx) e032 + (+atriPy) e013 + (+atriPz) e021", to_string(+a13));
	EXPECT_EQ("(-ae0) e0 + (-avx) e1 + (-avy) e2 + (-avz) e3 + (-ae123) e123 + (-atriPx) e032 + (-atriPy) e013 + (-atriPz) e021", to_string(-a13));
	EXPECT_EQ("(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3 + (ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021", to_string(a13 + b13));
	EXPECT_EQ("(ae0-be0) e0 + (avx-bvx) e1 + (avy-bvy) e2 + (avz-bvz) e3 + (ae123-be123) e123 + (atriPx-btriPx) e032 + (atriPy-btriPy) e013 + (atriPz-btriPz) e021", to_string(a13 - b13));
	EXPECT_EQ("(ae0*bs) e0 + (avx*bs) e1 + (avy*bs) e2 + (avz*bs) e3 + (ae123*bs) e123 + (atriPx*bs) e032 + (atriPy*bs) e013 + (atriPz*bs) e021", to_string(a13 * make_bs()));
	EXPECT_EQ("(ae0/bs) e0 + (avx/bs) e1 + (avy/bs) e2 + (avz/bs) e3 + (ae123/bs) e123 + (atriPx/bs) e032 + (atriPy/bs) e013 + (atriPz/bs) e021", to_string(a13 / make_bs()));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3 + (-ae123) e123 + (-atriPx) e032 + (-atriPy) e013 + (-atriPz) e021", to_string(~a13));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3 + (-ae123) e123 + (-atriPx) e032 + (-atriPy) e013 + (-atriPz) e021", to_string(a13.rev()));
	EXPECT_EQ("ae123 e0 + atriPx e1 + atriPy e2 + atriPz e3 + ae0 e123 + avx e032 + avy e013 + avz e021", to_string(!a13));
}

TEST(pga_test, operators_inplace_blade13) {
	auto a13 = make_a13();
	auto b13 = make_b13();
	auto r = a13; r += b13;
	EXPECT_EQ("(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3 + (ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021", to_string(r));
	r = a13; r -= b13;
	EXPECT_EQ("(ae0-be0) e0 + (avx-bvx) e1 + (avy-bvy) e2 + (avz-bvz) e3 + (ae123-be123) e123 + (atriPx-btriPx) e032 + (atriPy-btriPy) e013 + (atriPz-btriPz) e021", to_string(r));
	r = a13; r *= make_bs();
	EXPECT_EQ("(ae0*bs) e0 + (avx*bs) e1 + (avy*bs) e2 + (avz*bs) e3 + (ae123*bs) e123 + (atriPx*bs) e032 + (atriPy*bs) e013 + (atriPz*bs) e021", to_string(r));
	r = a13; r /= make_bs();
	EXPECT_EQ("(ae0/bs) e0 + (avx/bs) e1 + (avy/bs) e2 + (avz/bs) e3 + (ae123/bs) e123 + (atriPx/bs) e032 + (atriPy/bs) e013 + (atriPz/bs) e021", to_string(r));
	r = a13; r += a13;
	EXPECT_EQ("(ae0+ae0) e0 + (avx+avx) e1 + (avy+avy) e2 + (avz+avz) e3 + (ae123+ae123) e123 + (atriPx+atriPx) e032 + (atriPy+atriPy) e013 + (atriPz+atriPz) e021", to_string(r));
	r = a13; r -= a13;
	EXPECT_EQ("(ae0-ae0) e0 + (avx-avx) e1 + (avy-avy) e2 + (avz-avz) e3 + (ae123-ae123) e123 + (atriPx-atriPx) e032 + (atriPy-atriPy) e013 + (atriPz-atriPz) e021", to_string(r));
}

TEST(pga_test, constructor_multivector) {
	pga::multivector<symbolic> dm;
	EXPECT_EQ("? id + 0 e23 + 0 e31 + 0 e12 + 0 e01 + 0 e02 + 0 e03 + ? e0123 + ? e0 + 0 e1 + 0 e2 + 0 e3 + ? e123 + 0 e032 + 0 e013 + 0 e021", to_string(dm));
	pga::multivector<symbolic> am13(make_a13());
	EXPECT_EQ(
		"0 id + 0 e23 + 0 e31 + 0 e12 + "
		"0 e01 + 0 e02 + 0 e03 + 0 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + "
		"ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021",
		to_string(am13));
	pga::multivector<symbolic> am024(make_a024());
	EXPECT_EQ(
		"as id + abiEx e23 + abiEy e31 + abiEz e12 + "
		"abiex e01 + abiey e02 + abiez e03 + ae0123 e0123 + "
		"0 e0 + 0 e1 + 0 e2 + 0 e3 + "
		"0 e123 + 0 e032 + 0 e013 + 0 e021",
		to_string(am024));
	pga::multivector<symbolic> am(make_a024(), make_a13());
	EXPECT_EQ(
		"as id + abiEx e23 + abiEy e31 + abiEz e12 + "
		"abiex e01 + abiey e02 + abiez e03 + ae0123 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + "
		"ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021",
		to_string(am));
	pga::multivector<symbolic> an(make_a1(), make_a02E(), make_a2e4(), make_a3());
	EXPECT_EQ(
		"as id + abiEx e23 + abiEy e31 + abiEz e12 + "
		"abiex e01 + abiey e02 + abiez e03 + ae0123 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + "
		"ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021",
		to_string(an));
	pga::multivector<symbolic> ao(make_a1(), make_a0(), make_a2E(), make_a2e(), make_a4(), make_a3());
	EXPECT_EQ(
		"as id + abiEx e23 + abiEy e31 + abiEz e12 + "
		"abiex e01 + abiey e02 + abiez e03 + ae0123 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + "
		"ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021",
		to_string(ao));
	pga::multivector<symbolic> ap(make_a0(), make_a1(), make_a2E(), make_a2e(), make_a3(), make_a4());
	EXPECT_EQ(
		"as id + abiEx e23 + abiEy e31 + abiEz e12 + "
		"abiex e01 + abiey e02 + abiez e03 + ae0123 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + "
		"ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021",
		to_string(ap));
}

TEST(pga_test, operators_arithmetic_multivector) {
	auto am = make_am();
	auto bm = make_bm();
	EXPECT_EQ(
		"(+as) id + (+abiEx) e23 + (+abiEy) e31 + (+abiEz) e12 + "
		"(+abiex) e01 + (+abiey) e02 + (+abiez) e03 + (+ae0123) e0123 + "
		"(+ae0) e0 + (+avx) e1 + (+avy) e2 + (+avz) e3 + "
		"(+ae123) e123 + (+atriPx) e032 + (+atriPy) e013 + (+atriPz) e021",
		to_string(+am));
	EXPECT_EQ(
		"(-as) id + (-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12 + "
		"(-abiex) e01 + (-abiey) e02 + (-abiez) e03 + (-ae0123) e0123 + "
		"(-ae0) e0 + (-avx) e1 + (-avy) e2 + (-avz) e3 + "
		"(-ae123) e123 + (-atriPx) e032 + (-atriPy) e013 + (-atriPz) e021",
		to_string(-am));
	EXPECT_EQ(
		"(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + "
		"(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123 + "
		"(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3 + "
		"(ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021",
		to_string(am + bm));
	EXPECT_EQ(
		"(as-bs) id + (abiEx-bbiEx) e23 + (abiEy-bbiEy) e31 + (abiEz-bbiEz) e12 + "
		"(abiex-bbiex) e01 + (abiey-bbiey) e02 + (abiez-bbiez) e03 + (ae0123-be0123) e0123 + "
		"(ae0-be0) e0 + (avx-bvx) e1 + (avy-bvy) e2 + (avz-bvz) e3 + "
		"(ae123-be123) e123 + (atriPx-btriPx) e032 + (atriPy-btriPy) e013 + (atriPz-btriPz) e021",
		to_string(am - bm));
	EXPECT_EQ(
		"(as*bs) id + (abiEx*bs) e23 + (abiEy*bs) e31 + (abiEz*bs) e12 + "
		"(abiex*bs) e01 + (abiey*bs) e02 + (abiez*bs) e03 + (ae0123*bs) e0123 + "
		"(ae0*bs) e0 + (avx*bs) e1 + (avy*bs) e2 + (avz*bs) e3 + "
		"(ae123*bs) e123 + (atriPx*bs) e032 + (atriPy*bs) e013 + (atriPz*bs) e021",
		to_string(am * make_bs()));
	EXPECT_EQ(
		"(as/bs) id + (abiEx/bs) e23 + (abiEy/bs) e31 + (abiEz/bs) e12 + "
		"(abiex/bs) e01 + (abiey/bs) e02 + (abiez/bs) e03 + (ae0123/bs) e0123 + "
		"(ae0/bs) e0 + (avx/bs) e1 + (avy/bs) e2 + (avz/bs) e3 + "
		"(ae123/bs) e123 + (atriPx/bs) e032 + (atriPy/bs) e013 + (atriPz/bs) e021",
		to_string(am / make_bs()));
	EXPECT_EQ(
		"as id + (-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12 + "
		"(-abiex) e01 + (-abiey) e02 + (-abiez) e03 + ae0123 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + "
		"(-ae123) e123 + (-atriPx) e032 + (-atriPy) e013 + (-atriPz) e021",
		to_string(~am));
	EXPECT_EQ(
		"as id + (-abiEx) e23 + (-abiEy) e31 + (-abiEz) e12 + "
		"(-abiex) e01 + (-abiey) e02 + (-abiez) e03 + ae0123 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + "
		"(-ae123) e123 + (-atriPx) e032 + (-atriPy) e013 + (-atriPz) e021",
		to_string(am.rev()));
	EXPECT_EQ(
		"ae0123 id + abiex e23 + abiey e31 + abiez e12 + "
		"abiEx e01 + abiEy e02 + abiEz e03 + as e0123 + "
		"ae123 e0 + atriPx e1 + atriPy e2 + atriPz e3 + "
		"ae0 e123 + avx e032 + avy e013 + avz e021",
		to_string(!am));
}

TEST(pga_test, operators_inplace_multivector) {
	auto am = make_am();
	auto bm = make_bm();
	auto r = am; r += bm;
	EXPECT_EQ(
		"(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + "
		"(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123 + "
		"(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3 + "
		"(ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021",
		to_string(r));
	r = am; r -= bm;
	EXPECT_EQ(
		"(as-bs) id + (abiEx-bbiEx) e23 + (abiEy-bbiEy) e31 + (abiEz-bbiEz) e12 + "
		"(abiex-bbiex) e01 + (abiey-bbiey) e02 + (abiez-bbiez) e03 + (ae0123-be0123) e0123 + "
		"(ae0-be0) e0 + (avx-bvx) e1 + (avy-bvy) e2 + (avz-bvz) e3 + "
		"(ae123-be123) e123 + (atriPx-btriPx) e032 + (atriPy-btriPy) e013 + (atriPz-btriPz) e021",
		to_string(r));
	r = am; r *= make_bs();
	EXPECT_EQ(
		"(as*bs) id + (abiEx*bs) e23 + (abiEy*bs) e31 + (abiEz*bs) e12 + "
		"(abiex*bs) e01 + (abiey*bs) e02 + (abiez*bs) e03 + (ae0123*bs) e0123 + "
		"(ae0*bs) e0 + (avx*bs) e1 + (avy*bs) e2 + (avz*bs) e3 + "
		"(ae123*bs) e123 + (atriPx*bs) e032 + (atriPy*bs) e013 + (atriPz*bs) e021",
		to_string(r));
	r = am; r /= make_bs();
	EXPECT_EQ(
		"(as/bs) id + (abiEx/bs) e23 + (abiEy/bs) e31 + (abiEz/bs) e12 + "
		"(abiex/bs) e01 + (abiey/bs) e02 + (abiez/bs) e03 + (ae0123/bs) e0123 + "
		"(ae0/bs) e0 + (avx/bs) e1 + (avy/bs) e2 + (avz/bs) e3 + "
		"(ae123/bs) e123 + (atriPx/bs) e032 + (atriPy/bs) e013 + (atriPz/bs) e021",
		to_string(r));
	r = am; r += am;
	EXPECT_EQ(
		"(as+as) id + (abiEx+abiEx) e23 + (abiEy+abiEy) e31 + (abiEz+abiEz) e12 + "
		"(abiex+abiex) e01 + (abiey+abiey) e02 + (abiez+abiez) e03 + (ae0123+ae0123) e0123 + "
		"(ae0+ae0) e0 + (avx+avx) e1 + (avy+avy) e2 + (avz+avz) e3 + "
		"(ae123+ae123) e123 + (atriPx+atriPx) e032 + (atriPy+atriPy) e013 + (atriPz+atriPz) e021",
		to_string(r));
	r = am; r -= am;
	EXPECT_EQ(
		"(as-as) id + (abiEx-abiEx) e23 + (abiEy-abiEy) e31 + (abiEz-abiEz) e12 + "
		"(abiex-abiex) e01 + (abiey-abiey) e02 + (abiez-abiez) e03 + (ae0123-ae0123) e0123 + "
		"(ae0-ae0) e0 + (avx-avx) e1 + (avy-avy) e2 + (avz-avz) e3 + "
		"(ae123-ae123) e123 + (atriPx-atriPx) e032 + (atriPy-atriPy) e013 + (atriPz-atriPz) e021",
		to_string(r));
}

TEST(pga_test, operators_dual) {
	EXPECT_EQ("as e0123", to_string(!make_a0()));
	EXPECT_EQ("ae0 e123 + avx e032 + avy e013 + avz e021", to_string(!make_a1()));
	EXPECT_EQ("abiEx e01 + abiEy e02 + abiEz e03", to_string(!make_a2E()));
	EXPECT_EQ("abiex e23 + abiey e31 + abiez e12", to_string(!make_a2e()));
	EXPECT_EQ("ae123 e0 + atriPx e1 + atriPy e2 + atriPz e3", to_string(!make_a3()));
	EXPECT_EQ("ae0123 id", to_string(!make_a4()));
	EXPECT_EQ("abiEx e01 + abiEy e02 + abiEz e03 + as e0123", to_string(!make_a02E()));
	EXPECT_EQ("abiex e23 + abiey e31 + abiez e12 + as e0123", to_string(!make_a02e()));
	EXPECT_EQ("abiex e23 + abiey e31 + abiez e12 + abiEx e01 + abiEy e02 + abiEz e03", to_string(!make_a22()));
	EXPECT_EQ("ae0123 id + abiEx e01 + abiEy e02 + abiEz e03", to_string(!make_a2E4()));
	EXPECT_EQ("ae0123 id + abiex e23 + abiey e31 + abiez e12", to_string(!make_a2e4()));
	EXPECT_EQ("ae0123 id + abiex e23 + abiey e31 + abiez e12 + abiEx e01 + abiEy e02 + abiEz e03 + as e0123", to_string(!make_a024()));
	EXPECT_EQ("ae123 e0 + atriPx e1 + atriPy e2 + atriPz e3 + ae0 e123 + avx e032 + avy e013 + avz e021", to_string(!make_a13()));
	EXPECT_EQ(
		"ae0123 id + abiex e23 + abiey e31 + abiez e12 + "
		"abiEx e01 + abiEy e02 + abiEz e03 + as e0123 + "
		"ae123 e0 + atriPx e1 + atriPy e2 + atriPz e3 + "
		"ae0 e123 + avx e032 + avy e013 + avz e021",
		to_string(!make_am()));
}

TEST(pga_test, get) {
	auto z = make_z();
	auto a0 = make_a0(); auto a1 = make_a1(); auto a2E = make_a2E();auto a2e = make_a2e(); auto a3 = make_a3(); auto a4 = make_a4();
	auto a02E = make_a02E(); auto a02e = make_a02e(); auto a22 = make_a22(); auto a2E4 = make_a2E4(); auto a2e4 = make_a2e4(); auto a024 = make_a024(); auto a13 = make_a13(); auto am = make_am();

	EXPECT_EQ("0", to_string(pga::get<decltype(z)>::b0(z)));
	EXPECT_EQ("0", to_string(pga::get<decltype(z)>::b1(z)));
	EXPECT_EQ("0", to_string(pga::get<decltype(z)>::b2E(z)));
	EXPECT_EQ("0", to_string(pga::get<decltype(z)>::b2e(z)));
	EXPECT_EQ("0", to_string(pga::get<decltype(z)>::b3(z)));
	EXPECT_EQ("0", to_string(pga::get<decltype(z)>::b4(z)));

	EXPECT_EQ("as id", to_string(pga::get<decltype(a0)>::b0(a0)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a0)>::b1(a0)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a0)>::b2E(a0)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a0)>::b2e(a0)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a0)>::b3(a0)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a0)>::b4(a0)));

	EXPECT_EQ("0", to_string(pga::get<decltype(a1)>::b0(a1)));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3", to_string(pga::get<decltype(a1)>::b1(a1)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a1)>::b2E(a1)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a1)>::b2e(a1)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a1)>::b3(a1)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a1)>::b4(a1)));

	EXPECT_EQ("0", to_string(pga::get<decltype(a2E)>::b0(a2E)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a2E)>::b1(a2E)));
	EXPECT_EQ("abiEx e23 + abiEy e31 + abiEz e12", to_string(pga::get<decltype(a2E)>::b2E(a2E)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a2E)>::b2e(a2E)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a2E)>::b3(a2E)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a2E)>::b4(a2E)));

	EXPECT_EQ("0", to_string(pga::get<decltype(a2e)>::b0(a2e)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a2e)>::b1(a2e)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a2e)>::b2E(a2e)));
	EXPECT_EQ("abiex e01 + abiey e02 + abiez e03", to_string(pga::get<decltype(a2e)>::b2e(a2e)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a2e)>::b3(a2e)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a2e)>::b4(a2e)));

	EXPECT_EQ("0", to_string(pga::get<decltype(a3)>::b0(a3)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a3)>::b1(a3)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a3)>::b2E(a3)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a3)>::b2e(a3)));
	EXPECT_EQ("ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(pga::get<decltype(a3)>::b3(a3)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a3)>::b4(a3)));

	EXPECT_EQ("0", to_string(pga::get<decltype(a4)>::b0(a4)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a4)>::b1(a4)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a4)>::b2E(a4)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a4)>::b2e(a4)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a4)>::b3(a4)));
	EXPECT_EQ("ae0123 e0123", to_string(pga::get<decltype(a4)>::b4(a4)));

	EXPECT_EQ("as id", to_string(pga::get<decltype(a02E)>::b0(a02E)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a02E)>::b1(a02E)));
	EXPECT_EQ("abiEx e23 + abiEy e31 + abiEz e12", to_string(pga::get<decltype(a02E)>::b2E(a02E)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a02E)>::b2e(a02E)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a02E)>::b3(a02E)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a02E)>::b4(a02E)));

	EXPECT_EQ("as id", to_string(pga::get<decltype(a02e)>::b0(a02e)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a02e)>::b1(a02e)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a02e)>::b2E(a02e)));
	EXPECT_EQ("abiex e01 + abiey e02 + abiez e03", to_string(pga::get<decltype(a02e)>::b2e(a02e)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a02e)>::b3(a02e)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a02e)>::b4(a02e)));

	EXPECT_EQ("0", to_string(pga::get<decltype(a22)>::b0(a22)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a22)>::b1(a22)));
	EXPECT_EQ("abiEx e23 + abiEy e31 + abiEz e12", to_string(pga::get<decltype(a22)>::b2E(a22)));
	EXPECT_EQ("abiex e01 + abiey e02 + abiez e03", to_string(pga::get<decltype(a22)>::b2e(a22)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a22)>::b3(a22)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a22)>::b4(a22)));

	EXPECT_EQ("0", to_string(pga::get<decltype(a2E4)>::b0(a2E4)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a2E4)>::b1(a2E4)));
	EXPECT_EQ("abiEx e23 + abiEy e31 + abiEz e12", to_string(pga::get<decltype(a2E4)>::b2E(a2E4)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a2E4)>::b2e(a2E4)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a2E4)>::b3(a2E4)));
	EXPECT_EQ("ae0123 e0123", to_string(pga::get<decltype(a2E4)>::b4(a2E4)));

	EXPECT_EQ("0", to_string(pga::get<decltype(a2e4)>::b0(a2e4)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a2e4)>::b1(a2e4)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a2e4)>::b2E(a2e4)));
	EXPECT_EQ("abiex e01 + abiey e02 + abiez e03", to_string(pga::get<decltype(a2e4)>::b2e(a2e4)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a2e4)>::b3(a2e4)));
	EXPECT_EQ("ae0123 e0123", to_string(pga::get<decltype(a2e4)>::b4(a2e4)));

	EXPECT_EQ("as id", to_string(pga::get<decltype(a024)>::b0(a024)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a024)>::b1(a024)));
	EXPECT_EQ("abiEx e23 + abiEy e31 + abiEz e12", to_string(pga::get<decltype(a024)>::b2E(a024)));
	EXPECT_EQ("abiex e01 + abiey e02 + abiez e03", to_string(pga::get<decltype(a024)>::b2e(a024)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a024)>::b3(a024)));
	EXPECT_EQ("ae0123 e0123", to_string(pga::get<decltype(a024)>::b4(a024)));

	EXPECT_EQ("0", to_string(pga::get<decltype(a13)>::b0(a13)));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3", to_string(pga::get<decltype(a13)>::b1(a13)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a13)>::b2E(a13)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a13)>::b2e(a13)));
	EXPECT_EQ("ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(pga::get<decltype(a13)>::b3(a13)));
	EXPECT_EQ("0", to_string(pga::get<decltype(a13)>::b4(a13)));

	EXPECT_EQ("as id", to_string(pga::get<decltype(am)>::b0(am)));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3", to_string(pga::get<decltype(am)>::b1(am)));
	EXPECT_EQ("abiEx e23 + abiEy e31 + abiEz e12", to_string(pga::get<decltype(am)>::b2E(am)));
	EXPECT_EQ("abiex e01 + abiey e02 + abiez e03", to_string(pga::get<decltype(am)>::b2e(am)));
	EXPECT_EQ("ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(pga::get<decltype(am)>::b3(am)));
	EXPECT_EQ("ae0123 e0123", to_string(pga::get<decltype(am)>::b4(am)));
}

TEST(pga_test, combine) {
	auto z = make_z();
	auto a0 = make_a0(); auto a1 = make_a1(); auto a2E = make_a2E(); auto a2e = make_a2e(); auto a3 = make_a3(); auto a4 = make_a4();
	auto a02 = make_a02E(); auto a22 = make_a22(); auto a24 = make_a2e4(); auto a024 = make_a024(); auto a13 = make_a13(); auto am = make_am();

	std::string sz = "0";
	std::string sz0 = "0 id";
	std::string sz1 = "0 e0 + 0 e1 + 0 e2 + 0 e3";
	std::string sz2E = "0 e23 + 0 e31 + 0 e12";
	std::string sz2e = "0 e01 + 0 e02 + 0 e03";
	std::string sz3 = "0 e123 + 0 e032 + 0 e013 + 0 e021";
	std::string sz4 = "0 e0123";

	std::string sa0 = "as id";
	std::string sa1 = "ae0 e0 + avx e1 + avy e2 + avz e3";
	std::string sa2E = "abiEx e23 + abiEy e31 + abiEz e12";
	std::string sa2e = "abiex e01 + abiey e02 + abiez e03";
	std::string sa3 = "ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021";
	std::string sa4 = "ae0123 e0123";

	EXPECT_EQ(vec_to_string({ sz }), to_string(pga::combine024(z, z, z, z)));
	EXPECT_EQ(vec_to_string({ sa0 }), to_string(pga::combine024(a0, z, z, z)));
	EXPECT_EQ(vec_to_string({ sa2E }), to_string(pga::combine024(z, a2E, z, z)));
	EXPECT_EQ(vec_to_string({ sa2e }), to_string(pga::combine024(z, z, a2e, z)));
	EXPECT_EQ(vec_to_string({ sa4 }), to_string(pga::combine024(z, z, z, a4)));
	EXPECT_EQ(vec_to_string({ sa0, sa2E }), to_string(pga::combine024(a0, a2E, z, z)));
	EXPECT_EQ(vec_to_string({ sa2E, sa2e }), to_string(pga::combine024(z, a2E, a2e, z)));
	EXPECT_EQ(vec_to_string({ sa2e, sa4 }), to_string(pga::combine024(z, z, a2e, a4)));
	EXPECT_EQ(vec_to_string({ sa2E, sa4 }), to_string(pga::combine024(z, a2E, z, a4)));
	EXPECT_EQ(vec_to_string({ sa0, sz2E, sz2e, sa4 }), to_string(pga::combine024(a0, z, z, a4)));
	EXPECT_EQ(vec_to_string({ sa0, sa2e }), to_string(pga::combine024(a0, z, a2e, z)));
	EXPECT_EQ(vec_to_string({ sz0, sa2E, sa2e, sa4 }), to_string(pga::combine024(z, a2E, a2e, a4)));
	EXPECT_EQ(vec_to_string({ sa0, sz2E, sa2e, sa4 }), to_string(pga::combine024(a0, z, a2e, a4)));
	EXPECT_EQ(vec_to_string({ sa0, sa2E, sz2e, sa4 }), to_string(pga::combine024(a0, a2E, z, a4)));
	EXPECT_EQ(vec_to_string({ sa0, sa2E, sa2e, sz4 }), to_string(pga::combine024(a0, a2E, a2e, z)));
	EXPECT_EQ(vec_to_string({ sa0, sa2E, sa2e, sa4 }), to_string(pga::combine024(a0, a2E, a2e, a4)));

	EXPECT_EQ(vec_to_string({ sz }), to_string(pga::combine13(z, z)));
	EXPECT_EQ(vec_to_string({ sa1 }), to_string(pga::combine13(a1, z)));
	EXPECT_EQ(vec_to_string({ sa3 }), to_string(pga::combine13(z, a3)));
	EXPECT_EQ(vec_to_string({ sa1, sa3 }), to_string(pga::combine13(a1, a3)));

	EXPECT_EQ(vec_to_string({ sz }), to_string(pga::combine_multivector(z, z)));
	EXPECT_EQ(vec_to_string({ sa0 }), to_string(pga::combine_multivector(a0, z)));
	EXPECT_EQ(vec_to_string({ sa2E }), to_string(pga::combine_multivector(a2E, z)));
	EXPECT_EQ(vec_to_string({ sa2e }), to_string(pga::combine_multivector(a2e, z)));
	EXPECT_EQ(vec_to_string({ sa4 }), to_string(pga::combine_multivector(a4, z)));
	EXPECT_EQ(vec_to_string({ sa0, sa2E }), to_string(pga::combine_multivector(a02, z)));
	EXPECT_EQ(vec_to_string({ sa2E, sa2e }), to_string(pga::combine_multivector(a22, z)));
	EXPECT_EQ(vec_to_string({ sa2e, sa4 }), to_string(pga::combine_multivector(a24, z)));
	EXPECT_EQ(vec_to_string({ sa0, sa2E, sa2e, sa4 }), to_string(pga::combine_multivector(a024, z)));

	EXPECT_EQ(vec_to_string({ sa1 }), to_string(pga::combine_multivector(z, a1)));
	EXPECT_EQ(vec_to_string({ sa0, sz2E, sz2e, sz4, sa1, sz3 }), to_string(pga::combine_multivector(a0, a1)));
	EXPECT_EQ(vec_to_string({ sz0, sa2E, sz2e, sz4, sa1, sz3 }), to_string(pga::combine_multivector(a2E, a1)));
	EXPECT_EQ(vec_to_string({ sz0, sz2E, sa2e, sz4, sa1, sz3 }), to_string(pga::combine_multivector(a2e, a1)));
	EXPECT_EQ(vec_to_string({ sz0, sz2E, sz2e, sa4, sa1, sz3 }), to_string(pga::combine_multivector(a4, a1)));
	EXPECT_EQ(vec_to_string({ sa0, sa2E, sz2e, sz4, sa1, sz3 }), to_string(pga::combine_multivector(a02, a1)));
	EXPECT_EQ(vec_to_string({ sz0, sa2E, sa2e, sz4, sa1, sz3 }), to_string(pga::combine_multivector(a22, a1)));
	EXPECT_EQ(vec_to_string({ sz0, sz2E, sa2e, sa4, sa1, sz3 }), to_string(pga::combine_multivector(a24, a1)));
	EXPECT_EQ(vec_to_string({ sa0, sa2E, sa2e, sa4, sa1, sz3 }), to_string(pga::combine_multivector(a024, a1)));

	EXPECT_EQ(vec_to_string({ sa3 }), to_string(pga::combine_multivector(z, a3)));
	EXPECT_EQ(vec_to_string({ sa0, sz2E, sz2e, sz4, sz1, sa3 }), to_string(pga::combine_multivector(a0, a3)));
	EXPECT_EQ(vec_to_string({ sz0, sa2E, sz2e, sz4, sz1, sa3 }), to_string(pga::combine_multivector(a2E, a3)));
	EXPECT_EQ(vec_to_string({ sz0, sz2E, sa2e, sz4, sz1, sa3 }), to_string(pga::combine_multivector(a2e, a3)));
	EXPECT_EQ(vec_to_string({ sz0, sz2E, sz2e, sa4, sz1, sa3 }), to_string(pga::combine_multivector(a4, a3)));
	EXPECT_EQ(vec_to_string({ sa0, sa2E, sz2e, sz4, sz1, sa3 }), to_string(pga::combine_multivector(a02, a3)));
	EXPECT_EQ(vec_to_string({ sz0, sa2E, sa2e, sz4, sz1, sa3 }), to_string(pga::combine_multivector(a22, a3)));
	EXPECT_EQ(vec_to_string({ sz0, sz2E, sa2e, sa4, sz1, sa3 }), to_string(pga::combine_multivector(a24, a3)));
	EXPECT_EQ(vec_to_string({ sa0, sa2E, sa2e, sa4, sz1, sa3 }), to_string(pga::combine_multivector(a024, a3)));

	EXPECT_EQ(vec_to_string({ sa1, sa3 }), to_string(pga::combine_multivector(z, a13)));
	EXPECT_EQ(vec_to_string({ sa0, sz2E, sz2e, sz4, sa1, sa3 }), to_string(pga::combine_multivector(a0, a13)));
	EXPECT_EQ(vec_to_string({ sz0, sa2E, sz2e, sz4, sa1, sa3 }), to_string(pga::combine_multivector(a2E, a13)));
	EXPECT_EQ(vec_to_string({ sz0, sz2E, sa2e, sz4, sa1, sa3 }), to_string(pga::combine_multivector(a2e, a13)));
	EXPECT_EQ(vec_to_string({ sz0, sz2E, sz2e, sa4, sa1, sa3 }), to_string(pga::combine_multivector(a4, a13)));
	EXPECT_EQ(vec_to_string({ sa0, sa2E, sz2e, sz4, sa1, sa3 }), to_string(pga::combine_multivector(a02, a13)));
	EXPECT_EQ(vec_to_string({ sz0, sa2E, sa2e, sz4, sa1, sa3 }), to_string(pga::combine_multivector(a22, a13)));
	EXPECT_EQ(vec_to_string({ sz0, sz2E, sa2e, sa4, sa1, sa3 }), to_string(pga::combine_multivector(a24, a13)));
	EXPECT_EQ(vec_to_string({ sa0, sa2E, sa2e, sa4, sa1, sa3 }), to_string(pga::combine_multivector(a024, a13)));
}

TEST(pga_test, operators_add) {
	auto z = make_z();
	auto a0 = make_a0(); auto a1 = make_a1(); auto a2E = make_a2E(); auto a2e = make_a2e(); auto a3 = make_a3(); auto a4 = make_a4();
	auto a02 = make_a02E(); auto a22 = make_a22(); auto a24 = make_a2e4(); auto a024 = make_a024(); auto a13 = make_a13(); auto am = make_am();
	auto b0 = make_b0(); auto b1 = make_b1(); auto b2E = make_b2E(); auto b2e = make_b2e(); auto b3 = make_b3(); auto b4 = make_b4();
	auto b02 = make_b02E(); auto b22 = make_b22(); auto b24 = make_b2e4(); auto b024 = make_b024(); auto b13 = make_b13(); auto bm = make_bm();
	
	EXPECT_EQ("0", to_string(z + z));
	EXPECT_EQ("bs id", to_string(z + b0));
	EXPECT_EQ("be0 e0 + bvx e1 + bvy e2 + bvz e3", to_string(z + b1));
	EXPECT_EQ("bbiEx e23 + bbiEy e31 + bbiEz e12", to_string(z + b2E));
	EXPECT_EQ("bbiex e01 + bbiey e02 + bbiez e03", to_string(z + b2e));
	EXPECT_EQ("be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(z + b3));
	EXPECT_EQ("be0123 e0123", to_string(z + b4));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12", to_string(z + b02));
	EXPECT_EQ("bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03", to_string(z + b22));
	EXPECT_EQ("bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(z + b24));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(z + b024));
	EXPECT_EQ("be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(z + b13));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(z + bm));

	EXPECT_EQ("as id", to_string(a0 + z));
	EXPECT_EQ("(as+bs) id", to_string(a0 + b0));
	EXPECT_EQ("as id + 0 e23 + 0 e31 + 0 e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + 0 e123 + 0 e032 + 0 e013 + 0 e021", to_string(a0 + b1));
	EXPECT_EQ("as id + bbiEx e23 + bbiEy e31 + bbiEz e12", to_string(a0 + b2E));
	EXPECT_EQ("as id + bbiex e01 + bbiey e02 + bbiez e03", to_string(a0 + b2e));
	EXPECT_EQ("as id + 0 e23 + 0 e31 + 0 e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123 + "
		"0 e0 + 0 e1 + 0 e2 + 0 e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a0 + b3));
	EXPECT_EQ("as id + 0 e23 + 0 e31 + 0 e12 + 0 e01 + 0 e02 + 0 e03 + be0123 e0123", to_string(a0 + b4));
	EXPECT_EQ("(as+bs) id + bbiEx e23 + bbiEy e31 + bbiEz e12", to_string(a0 + b02));
	EXPECT_EQ("as id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + 0 e0123", to_string(a0 + b22));
	EXPECT_EQ("as id + 0 e23 + 0 e31 + 0 e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(a0 + b24));
	EXPECT_EQ("(as+bs) id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(a0 + b024));
	EXPECT_EQ("as id + 0 e23 + 0 e31 + 0 e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a0 + b13));
	EXPECT_EQ("(as+bs) id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a0 + bm));

	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3", to_string(a1 + z));
	EXPECT_EQ("bs id + 0 e23 + 0 e31 + 0 e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + 0 e123 + 0 e032 + 0 e013 + 0 e021", to_string(a1 + b0));
	EXPECT_EQ("(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3", to_string(a1 + b1));
	EXPECT_EQ("0 id + bbiEx e23 + bbiEy e31 + bbiEz e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + 0 e123 + 0 e032 + 0 e013 + 0 e021", to_string(a1 + b2E));
	EXPECT_EQ("0 id + 0 e23 + 0 e31 + 0 e12 + bbiex e01 + bbiey e02 + bbiez e03 + 0 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + 0 e123 + 0 e032 + 0 e013 + 0 e021", to_string(a1 + b2e));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a1 + b3));
	EXPECT_EQ("0 id + 0 e23 + 0 e31 + 0 e12 + 0 e01 + 0 e02 + 0 e03 + be0123 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + 0 e123 + 0 e032 + 0 e013 + 0 e021", to_string(a1 + b4));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + 0 e123 + 0 e032 + 0 e013 + 0 e021", to_string(a1 + b02));
	EXPECT_EQ("0 id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + 0 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + 0 e123 + 0 e032 + 0 e013 + 0 e021", to_string(a1 + b22));
	EXPECT_EQ("0 id + 0 e23 + 0 e31 + 0 e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + 0 e123 + 0 e032 + 0 e013 + 0 e021", to_string(a1 + b24));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + 0 e123 + 0 e032 + 0 e013 + 0 e021", to_string(a1 + b024));
	EXPECT_EQ("(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a1 + b13));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123 + "
		"(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a1 + bm));

	EXPECT_EQ("abiEx e23 + abiEy e31 + abiEz e12", to_string(a2E + z));
	EXPECT_EQ("bs id + abiEx e23 + abiEy e31 + abiEz e12", to_string(a2E + b0));
	EXPECT_EQ("0 id + abiEx e23 + abiEy e31 + abiEz e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + 0 e123 + 0 e032 + 0 e013 + 0 e021", to_string(a2E + b1));
	EXPECT_EQ("(abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12", to_string(a2E + b2E));
	EXPECT_EQ("abiEx e23 + abiEy e31 + abiEz e12 + bbiex e01 + bbiey e02 + bbiez e03", to_string(a2E + b2e));
	EXPECT_EQ("0 id + abiEx e23 + abiEy e31 + abiEz e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123 + "
		"0 e0 + 0 e1 + 0 e2 + 0 e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a2E + b3));
	EXPECT_EQ("abiEx e23 + abiEy e31 + abiEz e12 + be0123 e0123", to_string(a2E + b4));
	EXPECT_EQ("bs id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12", to_string(a2E + b02));
	EXPECT_EQ("(abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + bbiex e01 + bbiey e02 + bbiez e03", to_string(a2E + b22));
	EXPECT_EQ("0 id + abiEx e23 + abiEy e31 + abiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(a2E + b24));
	EXPECT_EQ("bs id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(a2E + b024));
	EXPECT_EQ("0 id + abiEx e23 + abiEy e31 + abiEz e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a2E + b13));
	EXPECT_EQ("bs id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a2E + bm));

	EXPECT_EQ("abiex e01 + abiey e02 + abiez e03", to_string(a2e + z));
	EXPECT_EQ("bs id + abiex e01 + abiey e02 + abiez e03", to_string(a2e + b0));
	EXPECT_EQ("0 id + 0 e23 + 0 e31 + 0 e12 + abiex e01 + abiey e02 + abiez e03 + 0 e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + 0 e123 + 0 e032 + 0 e013 + 0 e021", to_string(a2e + b1));
	EXPECT_EQ("bbiEx e23 + bbiEy e31 + bbiEz e12 + abiex e01 + abiey e02 + abiez e03", to_string(a2e + b2E));
	EXPECT_EQ("(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03", to_string(a2e + b2e));
	EXPECT_EQ("0 id + 0 e23 + 0 e31 + 0 e12 + abiex e01 + abiey e02 + abiez e03 + 0 e0123 + "
		"0 e0 + 0 e1 + 0 e2 + 0 e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a2e + b3));
	EXPECT_EQ("abiex e01 + abiey e02 + abiez e03 + be0123 e0123", to_string(a2e + b4));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + abiex e01 + abiey e02 + abiez e03 + 0 e0123", to_string(a2e + b02));
	EXPECT_EQ("bbiEx e23 + bbiEy e31 + bbiEz e12 + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03", to_string(a2e + b22));
	EXPECT_EQ("(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + be0123 e0123", to_string(a2e + b24));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + be0123 e0123", to_string(a2e + b024));
	EXPECT_EQ("0 id + 0 e23 + 0 e31 + 0 e12 + abiex e01 + abiey e02 + abiez e03 + 0 e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a2e + b13));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + be0123 e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a2e + bm));

	EXPECT_EQ("ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a3 + z));
	EXPECT_EQ("bs id + 0 e23 + 0 e31 + 0 e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123 + "
		"0 e0 + 0 e1 + 0 e2 + 0 e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a3 + b0));
	EXPECT_EQ("be0 e0 + bvx e1 + bvy e2 + bvz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a3 + b1));
	EXPECT_EQ("0 id + bbiEx e23 + bbiEy e31 + bbiEz e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123 + "
		"0 e0 + 0 e1 + 0 e2 + 0 e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a3 + b2E));
	EXPECT_EQ("0 id + 0 e23 + 0 e31 + 0 e12 + bbiex e01 + bbiey e02 + bbiez e03 + 0 e0123 + "
		"0 e0 + 0 e1 + 0 e2 + 0 e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a3 + b2e));
	EXPECT_EQ("(ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021", to_string(a3 + b3));
	EXPECT_EQ("0 id + 0 e23 + 0 e31 + 0 e12 + 0 e01 + 0 e02 + 0 e03 + be0123 e0123 + "
		"0 e0 + 0 e1 + 0 e2 + 0 e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a3 + b4));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123 + "
		"0 e0 + 0 e1 + 0 e2 + 0 e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a3 + b02));
	EXPECT_EQ("0 id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + 0 e0123 + "
		"0 e0 + 0 e1 + 0 e2 + 0 e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a3 + b22));
	EXPECT_EQ("0 id + 0 e23 + 0 e31 + 0 e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123 + "
		"0 e0 + 0 e1 + 0 e2 + 0 e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a3 + b24));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123 + "
		"0 e0 + 0 e1 + 0 e2 + 0 e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a3 + b024));
	EXPECT_EQ("be0 e0 + bvx e1 + bvy e2 + bvz e3 + (ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021", to_string(a3 + b13));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + (ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021", to_string(a3 + bm));

	EXPECT_EQ("ae0123 e0123", to_string(a4 + z));
	EXPECT_EQ("bs id + 0 e23 + 0 e31 + 0 e12 + 0 e01 + 0 e02 + 0 e03 + ae0123 e0123", to_string(a4 + b0));
	EXPECT_EQ("0 id + 0 e23 + 0 e31 + 0 e12 + 0 e01 + 0 e02 + 0 e03 + ae0123 e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + 0 e123 + 0 e032 + 0 e013 + 0 e021", to_string(a4 + b1));
	EXPECT_EQ("bbiEx e23 + bbiEy e31 + bbiEz e12 + ae0123 e0123", to_string(a4 + b2E));
	EXPECT_EQ("bbiex e01 + bbiey e02 + bbiez e03 + ae0123 e0123", to_string(a4 + b2e));
	EXPECT_EQ("0 id + 0 e23 + 0 e31 + 0 e12 + 0 e01 + 0 e02 + 0 e03 + ae0123 e0123 + "
		"0 e0 + 0 e1 + 0 e2 + 0 e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a4 + b3));
	EXPECT_EQ("(ae0123+be0123) e0123", to_string(a4 + b4));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + 0 e01 + 0 e02 + 0 e03 + ae0123 e0123", to_string(a4 + b02));
	EXPECT_EQ("0 id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + ae0123 e0123", to_string(a4 + b22));
	EXPECT_EQ("bbiex e01 + bbiey e02 + bbiez e03 + (ae0123+be0123) e0123", to_string(a4 + b24));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + (ae0123+be0123) e0123", to_string(a4 + b024));
	EXPECT_EQ("0 id + 0 e23 + 0 e31 + 0 e12 + 0 e01 + 0 e02 + 0 e03 + ae0123 e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a4 + b13));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + (ae0123+be0123) e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a4 + bm));

	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12", to_string(a02 + z));
	EXPECT_EQ("(as+bs) id + abiEx e23 + abiEy e31 + abiEz e12", to_string(a02 + b0));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + 0 e123 + 0 e032 + 0 e013 + 0 e021", to_string(a02 + b1));
	EXPECT_EQ("as id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12", to_string(a02 + b2E));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + 0 e0123", to_string(a02 + b2e));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123 + "
		"0 e0 + 0 e1 + 0 e2 + 0 e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a02 + b3));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + 0 e01 + 0 e02 + 0 e03 + be0123 e0123", to_string(a02 + b4));
	EXPECT_EQ("(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12", to_string(a02 + b02));
	EXPECT_EQ("as id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + bbiex e01 + bbiey e02 + bbiez e03 + 0 e0123", to_string(a02 + b22));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(a02 + b24));
	EXPECT_EQ("(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123", to_string(a02 + b024));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a02 + b13));
	EXPECT_EQ("(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a02 + bm));

	EXPECT_EQ("abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03", to_string(a22 + z));
	EXPECT_EQ("bs id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + 0 e0123", to_string(a22 + b0));
	EXPECT_EQ("0 id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + 0 e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + 0 e123 + 0 e032 + 0 e013 + 0 e021", to_string(a22 + b1));
	EXPECT_EQ("(abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + abiex e01 + abiey e02 + abiez e03", to_string(a22 + b2E));
	EXPECT_EQ("abiEx e23 + abiEy e31 + abiEz e12 + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03", to_string(a22 + b2e));
	EXPECT_EQ("0 id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + 0 e0123 + "
		"0 e0 + 0 e1 + 0 e2 + 0 e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a22 + b3));
	EXPECT_EQ("0 id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + be0123 e0123", to_string(a22 + b4));
	EXPECT_EQ("bs id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + abiex e01 + abiey e02 + abiez e03 + 0 e0123", to_string(a22 + b02));
	EXPECT_EQ("(abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03", to_string(a22 + b22));
	EXPECT_EQ("0 id + abiEx e23 + abiEy e31 + abiEz e12 + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + be0123 e0123", to_string(a22 + b24));
	EXPECT_EQ("bs id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + "
		"(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + be0123 e0123", to_string(a22 + b024));
	EXPECT_EQ("0 id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + 0 e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a22 + b13));
	EXPECT_EQ("bs id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + "
		"(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + be0123 e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a22 + bm));

	EXPECT_EQ("abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(a24 + z));
	EXPECT_EQ("bs id + 0 e23 + 0 e31 + 0 e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(a24 + b0));
	EXPECT_EQ("0 id + 0 e23 + 0 e31 + 0 e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + 0 e123 + 0 e032 + 0 e013 + 0 e021", to_string(a24 + b1));
	EXPECT_EQ("0 id + bbiEx e23 + bbiEy e31 + bbiEz e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(a24 + b2E));
	EXPECT_EQ("(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + ae0123 e0123", to_string(a24 + b2e));
	EXPECT_EQ("0 id + 0 e23 + 0 e31 + 0 e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123 + "
		"0 e0 + 0 e1 + 0 e2 + 0 e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a24 + b3));
	EXPECT_EQ("abiex e01 + abiey e02 + abiez e03 + (ae0123+be0123) e0123", to_string(a24 + b4));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(a24 + b02));
	EXPECT_EQ("0 id + bbiEx e23 + bbiEy e31 + bbiEz e12 + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + ae0123 e0123", to_string(a24 + b22));
	EXPECT_EQ("(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123", to_string(a24 + b24));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123", to_string(a24 + b024));
	EXPECT_EQ("0 id + 0 e23 + 0 e31 + 0 e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a24 + b13));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a24 + bm));

	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(a024 + z));
	EXPECT_EQ("(as+bs) id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(a024 + b0));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + 0 e123 + 0 e032 + 0 e013 + 0 e021", to_string(a024 + b1));
	EXPECT_EQ("as id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(a024 + b2E));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + ae0123 e0123", to_string(a024 + b2e));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123 + "
		"0 e0 + 0 e1 + 0 e2 + 0 e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a024 + b3));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + (ae0123+be0123) e0123", to_string(a024 + b4));
	EXPECT_EQ("(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123", to_string(a024 + b02));
	EXPECT_EQ("as id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + ae0123 e0123", to_string(a024 + b22));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123", to_string(a024 + b24));
	EXPECT_EQ("(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + "
		"(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123", to_string(a024 + b024));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a024 + b13));
	EXPECT_EQ("(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + "
		"(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123 + "
		"be0 e0 + bvx e1 + bvy e2 + bvz e3 + be123 e123 + btriPx e032 + btriPy e013 + btriPz e021", to_string(a024 + bm));

	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a13 + z));
	EXPECT_EQ("bs id + 0 e23 + 0 e31 + 0 e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a13 + b0));
	EXPECT_EQ("(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a13 + b1));
	EXPECT_EQ("0 id + bbiEx e23 + bbiEy e31 + bbiEz e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a13 + b2E));
	EXPECT_EQ("0 id + 0 e23 + 0 e31 + 0 e12 + bbiex e01 + bbiey e02 + bbiez e03 + 0 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a13 + b2e));
	EXPECT_EQ("ae0 e0 + avx e1 + avy e2 + avz e3 + (ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021", to_string(a13 + b3));
	EXPECT_EQ("0 id + 0 e23 + 0 e31 + 0 e12 + 0 e01 + 0 e02 + 0 e03 + be0123 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a13 + b4));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + 0 e01 + 0 e02 + 0 e03 + 0 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a13 + b02));
	EXPECT_EQ("0 id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + 0 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a13 + b22));
	EXPECT_EQ("0 id + 0 e23 + 0 e31 + 0 e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a13 + b24));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(a13 + b024));
	EXPECT_EQ("(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3 + (ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021", to_string(a13 + b13));
	EXPECT_EQ("bs id + bbiEx e23 + bbiEy e31 + bbiEz e12 + bbiex e01 + bbiey e02 + bbiez e03 + be0123 e0123 + "
		"(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3 + (ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021", to_string(a13 + bm));

	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(am + z));
	EXPECT_EQ("(as+bs) id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(am + b0));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123 + "
		"(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(am + b1));
	EXPECT_EQ("as id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(am + b2E));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + ae0123 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(am + b2e));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + (ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021", to_string(am + b3));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + (ae0123+be0123) e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(am + b4));
	EXPECT_EQ("(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(am + b02));
	EXPECT_EQ("as id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + ae0123 e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(am + b22));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + (abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(am + b24));
	EXPECT_EQ("(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + "
		"(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123 + "
		"ae0 e0 + avx e1 + avy e2 + avz e3 + ae123 e123 + atriPx e032 + atriPy e013 + atriPz e021", to_string(am + b024));
	EXPECT_EQ("as id + abiEx e23 + abiEy e31 + abiEz e12 + abiex e01 + abiey e02 + abiez e03 + ae0123 e0123 + "
		"(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3 + (ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021", to_string(am + b13));
	EXPECT_EQ("(as+bs) id + (abiEx+bbiEx) e23 + (abiEy+bbiEy) e31 + (abiEz+bbiEz) e12 + "
		"(abiex+bbiex) e01 + (abiey+bbiey) e02 + (abiez+bbiez) e03 + (ae0123+be0123) e0123 + "
		"(ae0+be0) e0 + (avx+bvx) e1 + (avy+bvy) e2 + (avz+bvz) e3 + "
		"(ae123+be123) e123 + (atriPx+btriPx) e032 + (atriPy+btriPy) e013 + (atriPz+btriPz) e021", to_string(am + bm));
}

TEST(pga_test, operators_multiply) {
	auto z = make_z();
	auto as = make_as(); auto a0 = make_a0(); auto a1 = make_a1(); auto a2E = make_a2E(); auto a2e = make_a2e(); auto a3 = make_a3(); auto a4 = make_a4();
	auto a02E = make_a02E(); auto a02e = make_a02e(); auto a22 = make_a22(); auto a2E4 = make_a2E4(); auto a2e4 = make_a2e4(); auto a024 = make_a024(); auto a13 = make_a13(); auto am = make_am();
	auto bs = make_bs(); auto b0 = make_b0(); auto b1 = make_b1(); auto b2E = make_b2E(); auto b2e = make_b2e(); auto b3 = make_b3(); auto b4 = make_b4();
	auto b02E = make_b02E(); auto b02e = make_b02e(); auto b22 = make_b22(); auto b2E4 = make_b2E4(); auto b2e4 = make_b2e4(); auto b024 = make_b024(); auto b13 = make_b13(); auto bm = make_bm();

	// zero
	EXPECT_EQ("0", to_string(z * z));
	EXPECT_EQ("0", to_string(z * b0));
	EXPECT_EQ("0", to_string(z * b1));
	EXPECT_EQ("0", to_string(z * b2E));
	EXPECT_EQ("0", to_string(z * b2e));
	EXPECT_EQ("0", to_string(z * b3));
	EXPECT_EQ("0", to_string(z * b4));
	EXPECT_EQ("0", to_string(z * b02E));
	EXPECT_EQ("0", to_string(z * b02e));
	EXPECT_EQ("0", to_string(z * b22));
	EXPECT_EQ("0", to_string(z * b2E4));
	EXPECT_EQ("0", to_string(z * b2e4));
	EXPECT_EQ("0", to_string(z * b024));
	EXPECT_EQ("0", to_string(z * b13));
	EXPECT_EQ("0", to_string(z * bm));
	EXPECT_EQ("0", to_string(a0 * z));
	EXPECT_EQ("0", to_string(a1 * z));
	EXPECT_EQ("0", to_string(a2E * z));
	EXPECT_EQ("0", to_string(a2e * z));
	EXPECT_EQ("0", to_string(a3 * z));
	EXPECT_EQ("0", to_string(a4 * z));
	EXPECT_EQ("0", to_string(a02E * z));
	EXPECT_EQ("0", to_string(a02e * z));
	EXPECT_EQ("0", to_string(a22 * z));
	EXPECT_EQ("0", to_string(a2E4 * z));
	EXPECT_EQ("0", to_string(a2e4 * z));
	EXPECT_EQ("0", to_string(a024 * z));
	EXPECT_EQ("0", to_string(am * z));

	// scalar (commutative)
	EXPECT_EQ("(bs*as) id", to_string(as * b0));
	EXPECT_EQ("(be0*as) e0 + (bvx*as) e1 + (bvy*as) e2 + (bvz*as) e3", to_string(as * b1));
	EXPECT_EQ("(bbiEx*as) e23 + (bbiEy*as) e31 + (bbiEz*as) e12", to_string(as * b2E));
	EXPECT_EQ("(bbiex*as) e01 + (bbiey*as) e02 + (bbiez*as) e03", to_string(as * b2e));
	EXPECT_EQ("(be123*as) e123 + (btriPx*as) e032 + (btriPy*as) e013 + (btriPz*as) e021", to_string(as * b3));
	EXPECT_EQ("(be0123*as) e0123", to_string(as * b4));
	EXPECT_EQ("(bs*as) id + (bbiEx*as) e23 + (bbiEy*as) e31 + (bbiEz*as) e12", to_string(as * b02E));
	EXPECT_EQ("(bs*as) id + (bbiex*as) e01 + (bbiey*as) e02 + (bbiez*as) e03", to_string(as * b02e));
	EXPECT_EQ("(bbiEx*as) e23 + (bbiEy*as) e31 + (bbiEz*as) e12 + (bbiex*as) e01 + (bbiey*as) e02 + (bbiez*as) e03", to_string(as * b22));
	EXPECT_EQ("(bbiEx*as) e23 + (bbiEy*as) e31 + (bbiEz*as) e12 + (be0123*as) e0123", to_string(as * b2E4));
	EXPECT_EQ("(bbiex*as) e01 + (bbiey*as) e02 + (bbiez*as) e03 + (be0123*as) e0123", to_string(as * b2e4));
	EXPECT_EQ("(bs*as) id + (bbiEx*as) e23 + (bbiEy*as) e31 + (bbiEz*as) e12 + (bbiex*as) e01 + (bbiey*as) e02 + (bbiez*as) e03 + (be0123*as) e0123", to_string(as * b024));
	EXPECT_EQ("(be0*as) e0 + (bvx*as) e1 + (bvy*as) e2 + (bvz*as) e3 + (be123*as) e123 + (btriPx*as) e032 + (btriPy*as) e013 + (btriPz*as) e021", to_string(as * b13));
	EXPECT_EQ("(bs*as) id + (bbiEx*as) e23 + (bbiEy*as) e31 + (bbiEz*as) e12 + (bbiex*as) e01 + (bbiey*as) e02 + (bbiez*as) e03 + (be0123*as) e0123 + "
		"(be0*as) e0 + (bvx*as) e1 + (bvy*as) e2 + (bvz*as) e3 + (be123*as) e123 + (btriPx*as) e032 + (btriPy*as) e013 + (btriPz*as) e021", to_string(as * bm));

	// primitive
	EXPECT_EQ("(as*bs) id", to_string(a0 * b0));
	EXPECT_EQ("(as*be0) e0 + (as*bvx) e1 + (as*bvy) e2 + (as*bvz) e3", to_string(a0 * b1));
	EXPECT_EQ("(as*bbiEx) e23 + (as*bbiEy) e31 + (as*bbiEz) e12", to_string(a0 * b2E));
	EXPECT_EQ("(as*bbiex) e01 + (as*bbiey) e02 + (as*bbiez) e03", to_string(a0 * b2e));
	EXPECT_EQ("(as*be123) e123 + (as*btriPx) e032 + (as*btriPy) e013 + (as*btriPz) e021", to_string(a0 * b3));
	EXPECT_EQ("(as*be0123) e0123", to_string(a0 * b4));
	EXPECT_EQ("(ae0*bs) e0 + (avx*bs) e1 + (avy*bs) e2 + (avz*bs) e3", to_string(a1 * b0));
	EXPECT_EQ("(((avx*bvx)+(avy*bvy))+(avz*bvz)) id + ((avy*bvz)-(avz*bvy)) e23 + ((avz*bvx)-(avx*bvz)) e31 + ((avx*bvy)-(avy*bvx)) e12 + ((ae0*bvx)-(avx*be0)) e01 + ((ae0*bvy)-(avy*be0)) e02 + ((ae0*bvz)-(avz*be0)) e03 + 0 e0123", to_string(a1 * b1));
	EXPECT_EQ("0 e0 + (-((avy*bbiEz)-(avz*bbiEy))) e1 + (-((avz*bbiEx)-(avx*bbiEz))) e2 + (-((avx*bbiEy)-(avy*bbiEx))) e3 + (((avx*bbiEx)+(avy*bbiEy))+(avz*bbiEz)) e123 + ((-ae0)*bbiEx) e032 + ((-ae0)*bbiEy) e013 + ((-ae0)*bbiEz) e021", to_string(a1 * b2E));
	EXPECT_EQ("(-(((avx*bbiex)+(avy*bbiey))+(avz*bbiez))) e0 + 0 e1 + 0 e2 + 0 e3 + 0 e123 + ((avy*bbiez)-(avz*bbiey)) e032 + ((avz*bbiex)-(avx*bbiez)) e013 + ((avx*bbiey)-(avy*bbiex)) e021", to_string(a1 * b2e));
	EXPECT_EQ("0 id + (avx*be123) e23 + (avy*be123) e31 + (avz*be123) e12 + (-((avy*btriPz)-(avz*btriPy))) e01 + (-((avz*btriPx)-(avx*btriPz))) e02 + (-((avx*btriPy)-(avy*btriPx))) e03 + ((ae0*be123)+(((avx*btriPx)+(avy*btriPy))+(avz*btriPz))) e0123", to_string(a1 * b3));
	EXPECT_EQ("0 e123 + (avx*be0123) e032 + (avy*be0123) e013 + (avz*be0123) e021", to_string(a1 * b4));
	EXPECT_EQ("(abiEx*bs) e23 + (abiEy*bs) e31 + (abiEz*bs) e12", to_string(a2E * b0));
	EXPECT_EQ("0 e0 + (-((abiEy*bvz)-(abiEz*bvy))) e1 + (-((abiEz*bvx)-(abiEx*bvz))) e2 + (-((abiEx*bvy)-(abiEy*bvx))) e3 + (((abiEx*bvx)+(abiEy*bvy))+(abiEz*bvz)) e123 + (abiEx*(-be0)) e032 + (abiEy*(-be0)) e013 + (abiEz*(-be0)) e021", to_string(a2E * b1));
	EXPECT_EQ("(-(((abiEx*bbiEx)+(abiEy*bbiEy))+(abiEz*bbiEz))) id + (-((abiEy*bbiEz)-(abiEz*bbiEy))) e23 + (-((abiEz*bbiEx)-(abiEx*bbiEz))) e31 + (-((abiEx*bbiEy)-(abiEy*bbiEx))) e12", to_string(a2E * b2E));
	EXPECT_EQ("(-((abiEy*bbiez)-(abiEz*bbiey))) e01 + (-((abiEz*bbiex)-(abiEx*bbiez))) e02 + (-((abiEx*bbiey)-(abiEy*bbiex))) e03 + (((abiEx*bbiex)+(abiEy*bbiey))+(abiEz*bbiez)) e0123", to_string(a2E * b2e));
	EXPECT_EQ("(((abiEx*btriPx)+(abiEy*btriPy))+(abiEz*btriPz)) e0 + (abiEx*(-be123)) e1 + (abiEy*(-be123)) e2 + (abiEz*(-be123)) e3 + 0 e123 + (-((abiEy*btriPz)-(abiEz*btriPy))) e032 + (-((abiEz*btriPx)-(abiEx*btriPz))) e013 + (-((abiEx*btriPy)-(abiEy*btriPx))) e021", to_string(a2E * b3));
	EXPECT_EQ("(abiEx*(-be0123)) e01 + (abiEy*(-be0123)) e02 + (abiEz*(-be0123)) e03", to_string(a2E * b4));
	EXPECT_EQ("(abiex*bs) e01 + (abiey*bs) e02 + (abiez*bs) e03", to_string(a2e * b0));
	EXPECT_EQ("(((abiex*bvx)+(abiey*bvy))+(abiez*bvz)) e0 + 0 e1 + 0 e2 + 0 e3 + 0 e123 + (-((abiey*bvz)-(abiez*bvy))) e032 + (-((abiez*bvx)-(abiex*bvz))) e013 + (-((abiex*bvy)-(abiey*bvx))) e021", to_string(a2e * b1));
	EXPECT_EQ("(-((abiey*bbiEz)-(abiez*bbiEy))) e01 + (-((abiez*bbiEx)-(abiex*bbiEz))) e02 + (-((abiex*bbiEy)-(abiey*bbiEx))) e03 + (((abiex*bbiEx)+(abiey*bbiEy))+(abiez*bbiEz)) e0123", to_string(a2e * b2E));
	EXPECT_EQ("0", to_string(a2e * b2e));
	EXPECT_EQ("0 e123 + (abiex*(-be123)) e032 + (abiey*(-be123)) e013 + (abiez*(-be123)) e021", to_string(a2e * b3));
	EXPECT_EQ("0", to_string(a2e * b4));
	EXPECT_EQ("(ae123*bs) e123 + (atriPx*bs) e032 + (atriPy*bs) e013 + (atriPz*bs) e021", to_string(a3 * b0));
	EXPECT_EQ("0 id + (ae123*bvx) e23 + (ae123*bvy) e31 + (ae123*bvz) e12 + ((atriPy*bvz)-(atriPz*bvy)) e01 + ((atriPz*bvx)-(atriPx*bvz)) e02 + ((atriPx*bvy)-(atriPy*bvx)) e03 + (((-ae123)*be0)-(((atriPx*bvx)+(atriPy*bvy))+(atriPz*bvz))) e0123", to_string(a3 * b1));
	EXPECT_EQ("(((atriPx*bbiEx)+(atriPy*bbiEy))+(atriPz*bbiEz)) e0 + ((-ae123)*bbiEx) e1 + ((-ae123)*bbiEy) e2 + ((-ae123)*bbiEz) e3 + 0 e123 + (-((atriPy*bbiEz)-(atriPz*bbiEy))) e032 + (-((atriPz*bbiEx)-(atriPx*bbiEz))) e013 + (-((atriPx*bbiEy)-(atriPy*bbiEx))) e021", to_string(a3 * b2E));
	EXPECT_EQ("0 e123 + (ae123*bbiex) e032 + (ae123*bbiey) e013 + (ae123*bbiez) e021", to_string(a3 * b2e));
	EXPECT_EQ("((-ae123)*be123) id + ((atriPx*be123)-(ae123*btriPx)) e01 + ((atriPy*be123)-(ae123*btriPy)) e02 + ((atriPz*be123)-(ae123*btriPz)) e03", to_string(a3 * b3));
	EXPECT_EQ("(ae123*be0123) e0 + 0 e1 + 0 e2 + 0 e3", to_string(a3 * b4));
	EXPECT_EQ("(ae0123*bs) e0123", to_string(a4 * b0));
	EXPECT_EQ("0 e123 + ((-ae0123)*bvx) e032 + ((-ae0123)*bvy) e013 + ((-ae0123)*bvz) e021", to_string(a4 * b1));
	EXPECT_EQ("((-ae0123)*bbiEx) e01 + ((-ae0123)*bbiEy) e02 + ((-ae0123)*bbiEz) e03", to_string(a4 * b2E));
	EXPECT_EQ("0", to_string(a4 * b2e));
	EXPECT_EQ("((-ae0123)*be123) e0 + 0 e1 + 0 e2 + 0 e3", to_string(a4 * b3));
	EXPECT_EQ("0", to_string(a4 * b4));

	// multivector
	EXPECT_EQ(
		"((((as*bs)+(-(((abiEx*bbiEx)+(abiEy*bbiEy))+(abiEz*bbiEz))))+(0+0))+(((((avx*bvx)+(avy*bvy))+(avz*bvz))+0)+(0+((-ae123)*be123)))) id + "
		"((((abiEx*bs)+((as*bbiEx)+(-((abiEy*bbiEz)-(abiEz*bbiEy)))))+(0+0))+((((avy*bvz)-(avz*bvy))+(ae123*bvx))+(avx*be123))) e23 + "
		"((((abiEy*bs)+((as*bbiEy)+(-((abiEz*bbiEx)-(abiEx*bbiEz)))))+(0+0))+((((avz*bvx)-(avx*bvz))+(ae123*bvy))+(avy*be123))) e31 + "
		"((((abiEz*bs)+((as*bbiEz)+(-((abiEx*bbiEy)-(abiEy*bbiEx)))))+(0+0))+((((avx*bvy)-(avy*bvx))+(ae123*bvz))+(avz*be123))) e12 + "
		"((((abiex*bs)+((-((abiey*bbiEz)-(abiez*bbiEy)))+((-ae0123)*bbiEx)))+(((as*bbiex)+(-((abiEy*bbiez)-(abiEz*bbiey))))+(abiEx*(-be0123))))+((((ae0*bvx)-(avx*be0))+((atriPy*bvz)-(atriPz*bvy)))+((-((avy*btriPz)-(avz*btriPy)))+((atriPx*be123)-(ae123*btriPx))))) e01 + "
		"((((abiey*bs)+((-((abiez*bbiEx)-(abiex*bbiEz)))+((-ae0123)*bbiEy)))+(((as*bbiey)+(-((abiEz*bbiex)-(abiEx*bbiez))))+(abiEy*(-be0123))))+((((ae0*bvy)-(avy*be0))+((atriPz*bvx)-(atriPx*bvz)))+((-((avz*btriPx)-(avx*btriPz)))+((atriPy*be123)-(ae123*btriPy))))) e02 + "
		"((((abiez*bs)+((-((abiex*bbiEy)-(abiey*bbiEx)))+((-ae0123)*bbiEz)))+(((as*bbiez)+(-((abiEx*bbiey)-(abiEy*bbiex))))+(abiEz*(-be0123))))+((((ae0*bvz)-(avz*be0))+((atriPx*bvy)-(atriPy*bvx)))+((-((avx*btriPy)-(avy*btriPx)))+((atriPz*be123)-(ae123*btriPz))))) e03 + "
		"((((ae0123*bs)+(((abiex*bbiEx)+(abiey*bbiEy))+(abiez*bbiEz)))+((((abiEx*bbiex)+(abiEy*bbiey))+(abiEz*bbiez))+(as*be0123)))+((0+(((-ae123)*be0)-(((atriPx*bvx)+(atriPy*bvy))+(atriPz*bvz))))+((ae0*be123)+(((avx*btriPx)+(avy*btriPy))+(avz*btriPz))))) e0123 + "
		"((((ae0*bs)+(0+(((atriPx*bbiEx)+(atriPy*bbiEy))+(atriPz*bbiEz))))+((-(((avx*bbiex)+(avy*bbiey))+(avz*bbiez)))+(ae123*be0123)))+((((as*be0)+0)+(((abiex*bvx)+(abiey*bvy))+(abiez*bvz)))+((((abiEx*btriPx)+(abiEy*btriPy))+(abiEz*btriPz))+((-ae0123)*be123)))) e0 + "
		"((((avx*bs)+((-((avy*bbiEz)-(avz*bbiEy)))+((-ae123)*bbiEx)))+(0+0))+((((as*bvx)+(-((abiEy*bvz)-(abiEz*bvy))))+0)+((abiEx*(-be123))+0))) e1 + "
		"((((avy*bs)+((-((avz*bbiEx)-(avx*bbiEz)))+((-ae123)*bbiEy)))+(0+0))+((((as*bvy)+(-((abiEz*bvx)-(abiEx*bvz))))+0)+((abiEy*(-be123))+0))) e2 + "
		"((((avz*bs)+((-((avx*bbiEy)-(avy*bbiEx)))+((-ae123)*bbiEz)))+(0+0))+((((as*bvz)+(-((abiEx*bvy)-(abiEy*bvx))))+0)+((abiEz*(-be123))+0))) e3 + "
		"((((ae123*bs)+((((avx*bbiEx)+(avy*bbiEy))+(avz*bbiEz))+0))+((0+0)+0))+(((((abiEx*bvx)+(abiEy*bvy))+(abiEz*bvz))+(0+0))+(((as*be123)+0)+0))) e123 + "
		"((((atriPx*bs)+(((-ae0)*bbiEx)+(-((atriPy*bbiEz)-(atriPz*bbiEy)))))+((((avy*bbiez)-(avz*bbiey))+(ae123*bbiex))+(avx*be0123)))+(((abiEx*(-be0))+((-((abiey*bvz)-(abiez*bvy)))+((-ae0123)*bvx)))+(((as*btriPx)+(-((abiEy*btriPz)-(abiEz*btriPy))))+(abiex*(-be123))))) e032 + "
		"((((atriPy*bs)+(((-ae0)*bbiEy)+(-((atriPz*bbiEx)-(atriPx*bbiEz)))))+((((avz*bbiex)-(avx*bbiez))+(ae123*bbiey))+(avy*be0123)))+(((abiEy*(-be0))+((-((abiez*bvx)-(abiex*bvz)))+((-ae0123)*bvy)))+(((as*btriPy)+(-((abiEz*btriPx)-(abiEx*btriPz))))+(abiey*(-be123))))) e013 + "
		"((((atriPz*bs)+(((-ae0)*bbiEz)+(-((atriPx*bbiEy)-(atriPy*bbiEx)))))+((((avx*bbiey)-(avy*bbiex))+(ae123*bbiez))+(avz*be0123)))+(((abiEz*(-be0))+((-((abiex*bvy)-(abiey*bvx)))+((-ae0123)*bvz)))+(((as*btriPz)+(-((abiEx*btriPy)-(abiEy*btriPx))))+(abiez*(-be123))))) e021",
		to_string(am * bm));
}

TEST(pga_test, operators_wedge) {
	auto z = make_z();
	auto as = make_as(); auto a0 = make_a0(); auto a1 = make_a1(); auto a2E = make_a2E(); auto a2e = make_a2e(); auto a3 = make_a3(); auto a4 = make_a4();
	auto a02E = make_a02E(); auto a02e = make_a02e(); auto a22 = make_a22(); auto a2E4 = make_a2E4(); auto a2e4 = make_a2e4(); auto a024 = make_a024(); auto a13 = make_a13(); auto am = make_am();
	auto bs = make_bs(); auto b0 = make_b0(); auto b1 = make_b1(); auto b2E = make_b2E(); auto b2e = make_b2e(); auto b3 = make_b3(); auto b4 = make_b4();
	auto b02E = make_b02E(); auto b02e = make_b02e(); auto b22 = make_b22(); auto b2E4 = make_b2E4(); auto b2e4 = make_b2e4(); auto b024 = make_b024(); auto b13 = make_b13(); auto bm = make_bm();

	// zero
	EXPECT_EQ("0", to_string(z ^ z));
	EXPECT_EQ("0", to_string(z ^ b0));
	EXPECT_EQ("0", to_string(z ^ b1));
	EXPECT_EQ("0", to_string(z ^ b2E));
	EXPECT_EQ("0", to_string(z ^ b2e));
	EXPECT_EQ("0", to_string(z ^ b3));
	EXPECT_EQ("0", to_string(z ^ b4));
	EXPECT_EQ("0", to_string(z ^ b02E));
	EXPECT_EQ("0", to_string(z ^ b02e));
	EXPECT_EQ("0", to_string(z ^ b22));
	EXPECT_EQ("0", to_string(z ^ b2E4));
	EXPECT_EQ("0", to_string(z ^ b2e4));
	EXPECT_EQ("0", to_string(z ^ b024));
	EXPECT_EQ("0", to_string(z ^ b13));
	EXPECT_EQ("0", to_string(z ^ bm));
	EXPECT_EQ("0", to_string(a0 ^ z));
	EXPECT_EQ("0", to_string(a1 ^ z));
	EXPECT_EQ("0", to_string(a2E ^ z));
	EXPECT_EQ("0", to_string(a2e ^ z));
	EXPECT_EQ("0", to_string(a3 ^ z));
	EXPECT_EQ("0", to_string(a4 ^ z));
	EXPECT_EQ("0", to_string(a02E ^ z));
	EXPECT_EQ("0", to_string(a02e ^ z));
	EXPECT_EQ("0", to_string(a22 ^ z));
	EXPECT_EQ("0", to_string(a2E4 ^ z));
	EXPECT_EQ("0", to_string(a2e4 ^ z));
	EXPECT_EQ("0", to_string(a024 ^ z));
	EXPECT_EQ("0", to_string(am ^ z));

	// primitive
	EXPECT_EQ("(as*bs) id", to_string(a0 ^ b0));
	EXPECT_EQ("(as*be0) e0 + (as*bvx) e1 + (as*bvy) e2 + (as*bvz) e3", to_string(a0 ^ b1));
	EXPECT_EQ("(as*bbiEx) e23 + (as*bbiEy) e31 + (as*bbiEz) e12", to_string(a0 ^ b2E));
	EXPECT_EQ("(as*bbiex) e01 + (as*bbiey) e02 + (as*bbiez) e03", to_string(a0 ^ b2e));
	EXPECT_EQ("(as*be123) e123 + (as*btriPx) e032 + (as*btriPy) e013 + (as*btriPz) e021", to_string(a0 ^ b3));
	EXPECT_EQ("(as*be0123) e0123", to_string(a0 ^ b4));

	EXPECT_EQ("(ae0*bs) e0 + (avx*bs) e1 + (avy*bs) e2 + (avz*bs) e3", to_string(a1 ^ b0));
	EXPECT_EQ("((avy*bvz)-(avz*bvy)) e23 + ((avz*bvx)-(avx*bvz)) e31 + ((avx*bvy)-(avy*bvx)) e12 + ((ae0*bvx)-(avx*be0)) e01 + ((ae0*bvy)-(avy*be0)) e02 + ((ae0*bvz)-(avz*be0)) e03", to_string(a1 ^ b1));
	EXPECT_EQ("(((avx*bbiEx)+(avy*bbiEy))+(avz*bbiEz)) e123 + ((-ae0)*bbiEx) e032 + ((-ae0)*bbiEy) e013 + ((-ae0)*bbiEz) e021", to_string(a1 ^ b2E));
	EXPECT_EQ("0 e123 + ((avy*bbiez)-(avz*bbiey)) e032 + ((avz*bbiex)-(avx*bbiez)) e013 + ((avx*bbiey)-(avy*bbiex)) e021", to_string(a1 ^ b2e));
	EXPECT_EQ("((ae0*be123)+(((avx*btriPx)+(avy*btriPy))+(avz*btriPz))) e0123", to_string(a1 ^ b3));
	EXPECT_EQ("0", to_string(a1 ^ b4));

	EXPECT_EQ("(abiEx*bs) e23 + (abiEy*bs) e31 + (abiEz*bs) e12", to_string(a2E ^ b0));
	EXPECT_EQ("(((abiEx*bvx)+(abiEy*bvy))+(abiEz*bvz)) e123 + (abiEx*(-be0)) e032 + (abiEy*(-be0)) e013 + (abiEz*(-be0)) e021", to_string(a2E ^ b1));
	EXPECT_EQ("0", to_string(a2E ^ b2E));
	EXPECT_EQ("(((abiEx*bbiex)+(abiEy*bbiey))+(abiEz*bbiez)) e0123", to_string(a2E ^ b2e));
	EXPECT_EQ("0", to_string(a2E ^ b3));
	EXPECT_EQ("0", to_string(a2E ^ b4));

	EXPECT_EQ("(abiex*bs) e01 + (abiey*bs) e02 + (abiez*bs) e03", to_string(a2e ^ b0));
	EXPECT_EQ("0 e123 + (-((abiey*bvz)-(abiez*bvy))) e032 + (-((abiez*bvx)-(abiex*bvz))) e013 + (-((abiex*bvy)-(abiey*bvx))) e021", to_string(a2e ^ b1));
	EXPECT_EQ("(((abiex*bbiEx)+(abiey*bbiEy))+(abiez*bbiEz)) e0123", to_string(a2e ^ b2E));
	EXPECT_EQ("0", to_string(a2e ^ b2e));
	EXPECT_EQ("0", to_string(a2e ^ b3));
	EXPECT_EQ("0", to_string(a2e ^ b4));

	EXPECT_EQ("(ae123*bs) e123 + (atriPx*bs) e032 + (atriPy*bs) e013 + (atriPz*bs) e021", to_string(a3 ^ b0));
	EXPECT_EQ("(((-ae123)*be0)-(((atriPx*bvx)+(atriPy*bvy))+(atriPz*bvz))) e0123", to_string(a3 ^ b1));
	EXPECT_EQ("0", to_string(a3 ^ b2E));
	EXPECT_EQ("0", to_string(a3 ^ b2e));
	EXPECT_EQ("0", to_string(a3 ^ b3));
	EXPECT_EQ("0", to_string(a3 ^ b4));

	EXPECT_EQ("(ae0123*bs) e0123", to_string(a4 ^ b0));
	EXPECT_EQ("0", to_string(a4 ^ b1));
	EXPECT_EQ("0", to_string(a4 ^ b2E));
	EXPECT_EQ("0", to_string(a4 ^ b2e));
	EXPECT_EQ("0", to_string(a4 ^ b3));
	EXPECT_EQ("0", to_string(a4 ^ b4));

	// multivector
	EXPECT_EQ(
		"((((as*bs)+0)+0)+(0+0)) id + "
		"((((abiEx*bs)+(as*bbiEx))+0)+(((avy*bvz)-(avz*bvy))+0)) e23 + "
		"((((abiEy*bs)+(as*bbiEy))+0)+(((avz*bvx)-(avx*bvz))+0)) e31 + "
		"((((abiEz*bs)+(as*bbiEz))+0)+(((avx*bvy)-(avy*bvx))+0)) e12 + "
		"((((abiex*bs)+0)+(as*bbiex))+(((ae0*bvx)-(avx*be0))+0)) e01 + "
		"((((abiey*bs)+0)+(as*bbiey))+(((ae0*bvy)-(avy*be0))+0)) e02 + "
		"((((abiez*bs)+0)+(as*bbiez))+(((ae0*bvz)-(avz*be0))+0)) e03 + "
		"((((ae0123*bs)+(((abiex*bbiEx)+(abiey*bbiEy))+(abiez*bbiEz)))+((((abiEx*bbiex)+(abiEy*bbiey))+(abiEz*bbiez))+(as*be0123)))+((((-ae123)*be0)-(((atriPx*bvx)+(atriPy*bvy))+(atriPz*bvz)))+((ae0*be123)+(((avx*btriPx)+(avy*btriPy))+(avz*btriPz))))) e0123 + "
		"((((ae0*bs)+0)+0)+((as*be0)+0)) e0 + "
		"((((avx*bs)+0)+0)+((as*bvx)+0)) e1 + "
		"((((avy*bs)+0)+0)+((as*bvy)+0)) e2 + "
		"((((avz*bs)+0)+0)+((as*bvz)+0)) e3 + "
		"((((ae123*bs)+(((avx*bbiEx)+(avy*bbiEy))+(avz*bbiEz)))+0)+(((((abiEx*bvx)+(abiEy*bvy))+(abiEz*bvz))+0)+(as*be123))) e123 + "
		"((((atriPx*bs)+((-ae0)*bbiEx))+((avy*bbiez)-(avz*bbiey)))+(((abiEx*(-be0))+(-((abiey*bvz)-(abiez*bvy))))+(as*btriPx))) e032 + "
		"((((atriPy*bs)+((-ae0)*bbiEy))+((avz*bbiex)-(avx*bbiez)))+(((abiEy*(-be0))+(-((abiez*bvx)-(abiex*bvz))))+(as*btriPy))) e013 + "
		"((((atriPz*bs)+((-ae0)*bbiEz))+((avx*bbiey)-(avy*bbiex)))+(((abiEz*(-be0))+(-((abiex*bvy)-(abiey*bvx))))+(as*btriPz))) e021",
		to_string(am ^ bm));
}

TEST(pga_test, operators_dot) {
	auto z = make_z();
	auto as = make_as(); auto a0 = make_a0(); auto a1 = make_a1(); auto a2E = make_a2E(); auto a2e = make_a2e(); auto a3 = make_a3(); auto a4 = make_a4();
	auto a02E = make_a02E(); auto a02e = make_a02e(); auto a22 = make_a22(); auto a2E4 = make_a2E4(); auto a2e4 = make_a2e4(); auto a024 = make_a024(); auto a13 = make_a13(); auto am = make_am();
	auto bs = make_bs(); auto b0 = make_b0(); auto b1 = make_b1(); auto b2E = make_b2E(); auto b2e = make_b2e(); auto b3 = make_b3(); auto b4 = make_b4();
	auto b02E = make_b02E(); auto b02e = make_b02e(); auto b22 = make_b22(); auto b2E4 = make_b2E4(); auto b2e4 = make_b2e4(); auto b024 = make_b024(); auto b13 = make_b13(); auto bm = make_bm();

	// zero
	EXPECT_EQ("0", to_string(z & z));
	EXPECT_EQ("0", to_string(z & b0));
	EXPECT_EQ("0", to_string(z & b1));
	EXPECT_EQ("0", to_string(z & b2E));
	EXPECT_EQ("0", to_string(z & b2e));
	EXPECT_EQ("0", to_string(z & b3));
	EXPECT_EQ("0", to_string(z & b4));
	EXPECT_EQ("0", to_string(z & b02E));
	EXPECT_EQ("0", to_string(z & b02e));
	EXPECT_EQ("0", to_string(z & b22));
	EXPECT_EQ("0", to_string(z & b2E4));
	EXPECT_EQ("0", to_string(z & b2e4));
	EXPECT_EQ("0", to_string(z & b024));
	EXPECT_EQ("0", to_string(z & b13));
	EXPECT_EQ("0", to_string(z & bm));
	EXPECT_EQ("0", to_string(a0 & z));
	EXPECT_EQ("0", to_string(a1 & z));
	EXPECT_EQ("0", to_string(a2E & z));
	EXPECT_EQ("0", to_string(a2e & z));
	EXPECT_EQ("0", to_string(a3 & z));
	EXPECT_EQ("0", to_string(a4 & z));
	EXPECT_EQ("0", to_string(a02E & z));
	EXPECT_EQ("0", to_string(a02e & z));
	EXPECT_EQ("0", to_string(a22 & z));
	EXPECT_EQ("0", to_string(a2E4 & z));
	EXPECT_EQ("0", to_string(a2e4 & z));
	EXPECT_EQ("0", to_string(a024 & z));
	EXPECT_EQ("0", to_string(am & z));

	// primitive
	EXPECT_EQ("(as*bs) id", to_string(a0 & b0));
	EXPECT_EQ("(as*be0) e0 + (as*bvx) e1 + (as*bvy) e2 + (as*bvz) e3", to_string(a0 & b1));
	EXPECT_EQ("(as*bbiEx) e23 + (as*bbiEy) e31 + (as*bbiEz) e12", to_string(a0 & b2E));
	EXPECT_EQ("(as*bbiex) e01 + (as*bbiey) e02 + (as*bbiez) e03", to_string(a0 & b2e));
	EXPECT_EQ("(as*be123) e123 + (as*btriPx) e032 + (as*btriPy) e013 + (as*btriPz) e021", to_string(a0 & b3));
	EXPECT_EQ("(as*be0123) e0123", to_string(a0 & b4));
	EXPECT_EQ("(ae0*bs) e0 + (avx*bs) e1 + (avy*bs) e2 + (avz*bs) e3", to_string(a1 & b0));
	EXPECT_EQ("(((avx*bvx)+(avy*bvy))+(avz*bvz)) id", to_string(a1 & b1));
	EXPECT_EQ("0 e0 + (-((avy*bbiEz)-(avz*bbiEy))) e1 + (-((avz*bbiEx)-(avx*bbiEz))) e2 + (-((avx*bbiEy)-(avy*bbiEx))) e3", to_string(a1 & b2E));
	EXPECT_EQ("(-(((avx*bbiex)+(avy*bbiey))+(avz*bbiez))) e0 + 0 e1 + 0 e2 + 0 e3", to_string(a1 & b2e));
	EXPECT_EQ("(avx*be123) e23 + (avy*be123) e31 + (avz*be123) e12 + (-((avy*btriPz)-(avz*btriPy))) e01 + (-((avz*btriPx)-(avx*btriPz))) e02 + (-((avx*btriPy)-(avy*btriPx))) e03", to_string(a1 & b3));
	EXPECT_EQ("0 e123 + (avx*be0123) e032 + (avy*be0123) e013 + (avz*be0123) e021", to_string(a1 & b4));
	EXPECT_EQ("(abiEx*bs) e23 + (abiEy*bs) e31 + (abiEz*bs) e12", to_string(a2E & b0));
	EXPECT_EQ("0 e0 + (-((abiEy*bvz)-(abiEz*bvy))) e1 + (-((abiEz*bvx)-(abiEx*bvz))) e2 + (-((abiEx*bvy)-(abiEy*bvx))) e3", to_string(a2E & b1));
	EXPECT_EQ("(-(((abiEx*bbiEx)+(abiEy*bbiEy))+(abiEz*bbiEz))) id", to_string(a2E & b2E));
	EXPECT_EQ("0", to_string(a2E & b2e));
	EXPECT_EQ("(((abiEx*btriPx)+(abiEy*btriPy))+(abiEz*btriPz)) e0 + (abiEx*(-be123)) e1 + (abiEy*(-be123)) e2 + (abiEz*(-be123)) e3", to_string(a2E & b3));
	EXPECT_EQ("(abiEx*(-be0123)) e01 + (abiEy*(-be0123)) e02 + (abiEz*(-be0123)) e03", to_string(a2E & b4));
	EXPECT_EQ("(abiex*bs) e01 + (abiey*bs) e02 + (abiez*bs) e03", to_string(a2e & b0));
	EXPECT_EQ("(((abiex*bvx)+(abiey*bvy))+(abiez*bvz)) e0 + 0 e1 + 0 e2 + 0 e3", to_string(a2e & b1));
	EXPECT_EQ("0", to_string(a2e & b2E));
	EXPECT_EQ("0", to_string(a2e & b2e));
	EXPECT_EQ("0", to_string(a2e & b3));
	EXPECT_EQ("0", to_string(a2e & b4));
	EXPECT_EQ("(ae123*bs) e123 + (atriPx*bs) e032 + (atriPy*bs) e013 + (atriPz*bs) e021", to_string(a3 & b0));
	EXPECT_EQ("(ae123*bvx) e23 + (ae123*bvy) e31 + (ae123*bvz) e12 + ((atriPy*bvz)-(atriPz*bvy)) e01 + ((atriPz*bvx)-(atriPx*bvz)) e02 + ((atriPx*bvy)-(atriPy*bvx)) e03", to_string(a3 & b1));
	EXPECT_EQ("(((atriPx*bbiEx)+(atriPy*bbiEy))+(atriPz*bbiEz)) e0 + ((-ae123)*bbiEx) e1 + ((-ae123)*bbiEy) e2 + ((-ae123)*bbiEz) e3", to_string(a3 & b2E));
	EXPECT_EQ("0", to_string(a3 & b2e));
	EXPECT_EQ("((-ae123)*be123) id", to_string(a3 & b3));
	EXPECT_EQ("(ae123*be0123) e0 + 0 e1 + 0 e2 + 0 e3", to_string(a3 & b4));
	EXPECT_EQ("(ae0123*bs) e0123", to_string(a4 & b0));
	EXPECT_EQ("0 e123 + ((-ae0123)*bvx) e032 + ((-ae0123)*bvy) e013 + ((-ae0123)*bvz) e021", to_string(a4 & b1));
	EXPECT_EQ("((-ae0123)*bbiEx) e01 + ((-ae0123)*bbiEy) e02 + ((-ae0123)*bbiEz) e03", to_string(a4 & b2E));
	EXPECT_EQ("0", to_string(a4 & b2e));
	EXPECT_EQ("((-ae0123)*be123) e0 + 0 e1 + 0 e2 + 0 e3", to_string(a4 & b3));
	EXPECT_EQ("0", to_string(a4 & b4));

	// multivector
	EXPECT_EQ(
		"((((as*bs)+(-(((abiEx*bbiEx)+(abiEy*bbiEy))+(abiEz*bbiEz))))+(0+0))+((((avx*bvx)+(avy*bvy))+(avz*bvz))+((-ae123)*be123))) id + "
		"((((abiEx*bs)+(as*bbiEx))+(0+0))+((ae123*bvx)+(avx*be123))) e23 + "
		"((((abiEy*bs)+(as*bbiEy))+(0+0))+((ae123*bvy)+(avy*be123))) e31 + "
		"((((abiEz*bs)+(as*bbiEz))+(0+0))+((ae123*bvz)+(avz*be123))) e12 + "
		"((((abiex*bs)+((-ae0123)*bbiEx))+((as*bbiex)+(abiEx*(-be0123))))+(((atriPy*bvz)-(atriPz*bvy))+(-((avy*btriPz)-(avz*btriPy))))) e01 + "
		"((((abiey*bs)+((-ae0123)*bbiEy))+((as*bbiey)+(abiEy*(-be0123))))+(((atriPz*bvx)-(atriPx*bvz))+(-((avz*btriPx)-(avx*btriPz))))) e02 + "
		"((((abiez*bs)+((-ae0123)*bbiEz))+((as*bbiez)+(abiEz*(-be0123))))+(((atriPx*bvy)-(atriPy*bvx))+(-((avx*btriPy)-(avy*btriPx))))) e03 + "
		"((((ae0123*bs)+0)+(0+(as*be0123)))+(0+0)) e0123 + "
		"((((ae0*bs)+(0+(((atriPx*bbiEx)+(atriPy*bbiEy))+(atriPz*bbiEz))))+((-(((avx*bbiex)+(avy*bbiey))+(avz*bbiez)))+(ae123*be0123)))+((((as*be0)+0)+(((abiex*bvx)+(abiey*bvy))+(abiez*bvz)))+((((abiEx*btriPx)+(abiEy*btriPy))+(abiEz*btriPz))+((-ae0123)*be123)))) e0 + "
		"((((avx*bs)+((-((avy*bbiEz)-(avz*bbiEy)))+((-ae123)*bbiEx)))+(0+0))+((((as*bvx)+(-((abiEy*bvz)-(abiEz*bvy))))+0)+((abiEx*(-be123))+0))) e1 + "
		"((((avy*bs)+((-((avz*bbiEx)-(avx*bbiEz)))+((-ae123)*bbiEy)))+(0+0))+((((as*bvy)+(-((abiEz*bvx)-(abiEx*bvz))))+0)+((abiEy*(-be123))+0))) e2 + "
		"((((avz*bs)+((-((avx*bbiEy)-(avy*bbiEx)))+((-ae123)*bbiEz)))+(0+0))+((((as*bvz)+(-((abiEx*bvy)-(abiEy*bvx))))+0)+((abiEz*(-be123))+0))) e3 + "
		"((((ae123*bs)+0)+(0+0))+(0+(as*be123))) e123 + "
		"((((atriPx*bs)+0)+(0+(avx*be0123)))+(((-ae0123)*bvx)+(as*btriPx))) e032 + "
		"((((atriPy*bs)+0)+(0+(avy*be0123)))+(((-ae0123)*bvy)+(as*btriPy))) e013 + "
		"((((atriPz*bs)+0)+(0+(avz*be0123)))+(((-ae0123)*bvz)+(as*btriPz))) e021",
		to_string(am & bm));
}

TEST(pga_test, operators_join) {
	auto z = make_z();
	auto as = make_as(); auto a0 = make_a0(); auto a1 = make_a1(); auto a2E = make_a2E(); auto a2e = make_a2e(); auto a3 = make_a3(); auto a4 = make_a4();
	auto a02E = make_a02E(); auto a02e = make_a02e(); auto a22 = make_a22(); auto a2E4 = make_a2E4(); auto a2e4 = make_a2e4(); auto a024 = make_a024(); auto a13 = make_a13(); auto am = make_am();
	auto bs = make_bs(); auto b0 = make_b0(); auto b1 = make_b1(); auto b2E = make_b2E(); auto b2e = make_b2e(); auto b3 = make_b3(); auto b4 = make_b4();
	auto b02E = make_b02E(); auto b02e = make_b02e(); auto b22 = make_b22(); auto b2E4 = make_b2E4(); auto b2e4 = make_b2e4(); auto b024 = make_b024(); auto b13 = make_b13(); auto bm = make_bm();

	// zero
	EXPECT_EQ("0", to_string(z | z));
	EXPECT_EQ("0", to_string(z | b0));
	EXPECT_EQ("0", to_string(z | b1));
	EXPECT_EQ("0", to_string(z | b2E));
	EXPECT_EQ("0", to_string(z | b2e));
	EXPECT_EQ("0", to_string(z | b3));
	EXPECT_EQ("0", to_string(z | b4));
	EXPECT_EQ("0", to_string(z | b02E));
	EXPECT_EQ("0", to_string(z | b02e));
	EXPECT_EQ("0", to_string(z | b22));
	EXPECT_EQ("0", to_string(z | b2E4));
	EXPECT_EQ("0", to_string(z | b2e4));
	EXPECT_EQ("0", to_string(z | b024));
	EXPECT_EQ("0", to_string(z | b13));
	EXPECT_EQ("0", to_string(z | bm));
	EXPECT_EQ("0", to_string(a0 | z));
	EXPECT_EQ("0", to_string(a1 | z));
	EXPECT_EQ("0", to_string(a2E | z));
	EXPECT_EQ("0", to_string(a2e | z));
	EXPECT_EQ("0", to_string(a3 | z));
	EXPECT_EQ("0", to_string(a4 | z));
	EXPECT_EQ("0", to_string(a02E | z));
	EXPECT_EQ("0", to_string(a02e | z));
	EXPECT_EQ("0", to_string(a22 | z));
	EXPECT_EQ("0", to_string(a2E4 | z));
	EXPECT_EQ("0", to_string(a2e4 | z));
	EXPECT_EQ("0", to_string(a024 | z));
	EXPECT_EQ("0", to_string(am | z));

	// primitive
	EXPECT_EQ("0", to_string(a0 | b0));
	EXPECT_EQ("0", to_string(a0 | b1));
	EXPECT_EQ("0", to_string(a0 | b2E));
	EXPECT_EQ("0", to_string(a0 | b2e));
	EXPECT_EQ("0", to_string(a0 | b3));
	EXPECT_EQ("(as*be0123) id", to_string(a0 | b4));
	EXPECT_EQ("0", to_string(a1 | b0));
	EXPECT_EQ("0", to_string(a1 | b1));
	EXPECT_EQ("0", to_string(a1 | b2E));
	EXPECT_EQ("0", to_string(a1 | b2e));
	EXPECT_EQ("(((-ae0)*be123)-(((avx*btriPx)+(avy*btriPy))+(avz*btriPz))) id", to_string(a1 | b3));
	EXPECT_EQ("(ae0*be0123) e0 + (avx*be0123) e1 + (avy*be0123) e2 + (avz*be0123) e3", to_string(a1 | b4));
	EXPECT_EQ("0", to_string(a2E | b0));
	EXPECT_EQ("0", to_string(a2E | b1));
	EXPECT_EQ("0", to_string(a2E | b2E));
	EXPECT_EQ("(((abiEx*bbiex)+(abiEy*bbiey))+(abiEz*bbiez)) id", to_string(a2E | b2e));
	EXPECT_EQ("0 e0 + (-((abiEy*btriPz)-(abiEz*btriPy))) e1 + (-((abiEz*btriPx)-(abiEx*btriPz))) e2 + (-((abiEx*btriPy)-(abiEy*btriPx))) e3", to_string(a2E | b3));
	EXPECT_EQ("(abiEx*be0123) e23 + (abiEy*be0123) e31 + (abiEz*be0123) e12", to_string(a2E | b4));
	EXPECT_EQ("0", to_string(a2e | b0));
	EXPECT_EQ("0", to_string(a2e | b1));
	EXPECT_EQ("(((abiex*bbiEx)+(abiey*bbiEy))+(abiez*bbiEz)) id", to_string(a2e | b2E));
	EXPECT_EQ("0", to_string(a2e | b2e));
	EXPECT_EQ("(((abiex*btriPx)+(abiey*btriPy))+(abiez*btriPz)) e0 + (abiex*(-be123)) e1 + (abiey*(-be123)) e2 + (abiez*(-be123)) e3", to_string(a2e | b3));
	EXPECT_EQ("(abiex*be0123) e01 + (abiey*be0123) e02 + (abiez*be0123) e03", to_string(a2e | b4));
	EXPECT_EQ("0", to_string(a3 | b0));
	EXPECT_EQ("((ae123*be0)+(((atriPx*bvx)+(atriPy*bvy))+(atriPz*bvz))) id", to_string(a3 | b1));
	EXPECT_EQ("0 e0 + ((atriPy*bbiEz)-(atriPz*bbiEy)) e1 + ((atriPz*bbiEx)-(atriPx*bbiEz)) e2 + ((atriPx*bbiEy)-(atriPy*bbiEx)) e3", to_string(a3 | b2E));
	EXPECT_EQ("(((atriPx*bbiex)+(atriPy*bbiey))+(atriPz*bbiez)) e0 + ((-ae123)*bbiex) e1 + ((-ae123)*bbiey) e2 + ((-ae123)*bbiez) e3", to_string(a3 | b2e));
	EXPECT_EQ("((atriPx*(-be123))+(ae123*btriPx)) e23 + ((atriPy*(-be123))+(ae123*btriPy)) e31 + ((atriPz*(-be123))+(ae123*btriPz)) e12 + ((atriPy*btriPz)-(atriPz*btriPy)) e01 + ((atriPz*btriPx)-(atriPx*btriPz)) e02 + ((atriPx*btriPy)-(atriPy*btriPx)) e03", to_string(a3 | b3));
	EXPECT_EQ("(ae123*be0123) e123 + (atriPx*be0123) e032 + (atriPy*be0123) e013 + (atriPz*be0123) e021", to_string(a3 | b4));
	EXPECT_EQ("(ae0123*bs) id", to_string(a4 | b0));
	EXPECT_EQ("(ae0123*be0) e0 + (ae0123*bvx) e1 + (ae0123*bvy) e2 + (ae0123*bvz) e3", to_string(a4 | b1));
	EXPECT_EQ("(ae0123*bbiEx) e23 + (ae0123*bbiEy) e31 + (ae0123*bbiEz) e12", to_string(a4 | b2E));
	EXPECT_EQ("(ae0123*bbiex) e01 + (ae0123*bbiey) e02 + (ae0123*bbiez) e03", to_string(a4 | b2e));
	EXPECT_EQ("(ae0123*be123) e123 + (ae0123*btriPx) e032 + (ae0123*btriPy) e013 + (ae0123*btriPz) e021", to_string(a4 | b3));
	EXPECT_EQ("(ae0123*be0123) e0123", to_string(a4 | b4));

	// multivector
	EXPECT_EQ(
		"((((ae0123*bs)+(((abiex*bbiEx)+(abiey*bbiEy))+(abiez*bbiEz)))+((((abiEx*bbiex)+(abiEy*bbiey))+(abiEz*bbiez))+(as*be0123)))+(((ae123*be0)+(((atriPx*bvx)+(atriPy*bvy))+(atriPz*bvz)))+(((-ae0)*be123)-(((avx*btriPx)+(avy*btriPy))+(avz*btriPz))))) id + "
		"(((ae0123*bbiEx)+(0+(abiEx*be0123)))+(0+((atriPx*(-be123))+(ae123*btriPx)))) e23 + "
		"(((ae0123*bbiEy)+(0+(abiEy*be0123)))+(0+((atriPy*(-be123))+(ae123*btriPy)))) e31 + "
		"(((ae0123*bbiEz)+(0+(abiEz*be0123)))+(0+((atriPz*(-be123))+(ae123*btriPz)))) e12 + "
		"((0+((ae0123*bbiex)+(abiex*be0123)))+(0+((atriPy*btriPz)-(atriPz*btriPy)))) e01 + "
		"((0+((ae0123*bbiey)+(abiey*be0123)))+(0+((atriPz*btriPx)-(atriPx*btriPz)))) e02 + "
		"((0+((ae0123*bbiez)+(abiez*be0123)))+(0+((atriPx*btriPy)-(atriPy*btriPx)))) e03 + "
		"((0+(0+(ae0123*be0123)))+(0+0)) e0123 + "
		"((0+((((atriPx*bbiex)+(atriPy*bbiey))+(atriPz*bbiez))+(ae0*be0123)))+((ae0123*be0)+(0+(((abiex*btriPx)+(abiey*btriPy))+(abiez*btriPz))))) e0 + "
		"((((atriPy*bbiEz)-(atriPz*bbiEy))+(((-ae123)*bbiex)+(avx*be0123)))+((ae0123*bvx)+((-((abiEy*btriPz)-(abiEz*btriPy)))+(abiex*(-be123))))) e1 + "
		"((((atriPz*bbiEx)-(atriPx*bbiEz))+(((-ae123)*bbiey)+(avy*be0123)))+((ae0123*bvy)+((-((abiEz*btriPx)-(abiEx*btriPz)))+(abiey*(-be123))))) e2 + "
		"((((atriPx*bbiEy)-(atriPy*bbiEx))+(((-ae123)*bbiez)+(avz*be0123)))+((ae0123*bvz)+((-((abiEx*btriPy)-(abiEy*btriPx)))+(abiez*(-be123))))) e3 + "
		"((0+(0+(ae123*be0123)))+(0+(ae0123*be123))) e123 + "
		"((0+(0+(atriPx*be0123)))+(0+(ae0123*btriPx))) e032 + "
		"((0+(0+(atriPy*be0123)))+(0+(ae0123*btriPy))) e013 + "
		"((0+(0+(atriPz*be0123)))+(0+(ae0123*btriPz))) e021",
		to_string(am | bm));
}

TEST(pga_test, primitives) {
	auto pa = pga::plane<double>({ 2, 5, 14 }, 4);
	auto pb = pga::plane<double>({ -2, 6, 9 }, 3);
	auto pc = pga::plane<double>({ -3, 4, 12 }, 2);
	auto Pa = pga::point<double>({ 4, 1, 8 });
	auto Pb = pga::point<double>({ 2, -6, 3 });
	auto Pc = pga::point<double>({ -10, 11, 2 });
	auto le = pga::line<double>({ 2, -5, 5 }, { 2, 6, 3 });

	EXPECT_EQ("4 e0 + 2 e1 + 5 e2 + 14 e3", to_string(pa));
	EXPECT_EQ("1 e123 + -4 e032 + -1 e013 + -8 e021", to_string(Pa));
	EXPECT_EQ("2 e23 + -5 e31 + 5 e12 + -45 e01 + 4 e02 + 22 e03", to_string(le));
	
	auto lab = pa ^ pb;
	EXPECT_EQ("-39 e23 + -46 e31 + 22 e12 + -14 e01 + 9 e02 + -6 e03", to_string(lab));
	
	auto Pabc = pa ^ pb ^ pc;
	EXPECT_EQ("197 e123 + -54 e032 + -94 e013 + -15 e021", to_string(Pabc));

	auto pabc = Pa | Pb | Pc;
	EXPECT_EQ("518 e0 + -92 e1 + -58 e2 + 118 e3", to_string(pabc));

	auto l = pa ^ pabc;
	EXPECT_EQ("1402 e23 + -1524 e31 + 344 e12 + -1404 e01 + -2822 e02 + -6780 e03", to_string(l));
	
	auto pabc2 = l | Pa / 121;
	EXPECT_EQ("518 e0 + -92 e1 + -58 e2 + 118 e3", to_string(pabc2));
}

TEST(pga_test, operations) {
	auto t = pga::translator<double>({ 3, -5, 7 });
	auto r1 = pga::rotor<double>({ 9, -12, 20 }, 1.28700221758657);
	auto r2 = pga::rotor<double>({ 9, -12, 20 }, 0.8, 0.6);
	auto P = pga::point<double>({ 4, 1, 8 });
	auto Pt = t * P * t.rev();
	auto Pr1 = r1 * P * r1.rev();
	auto Pr2 = r2 * P * r2.rev();
	EXPECT_EQ("1 e123 + -7 e032 + 4 e013 + -15 e021", to_string(Pt));
	EXPECT_EQ("0 e0 + 0 e1 + 0 e2 + 0 e3 + 15625 e123 + 22292 e032 + 30569 e013 + -135440 e021", to_string(Pr1 * 15625));
	EXPECT_EQ("0 e0 + 0 e1 + 0 e2 + 0 e3 + 15625 e123 + 22292 e032 + 30569 e013 + -135440 e021", to_string(Pr2 * 15625));
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
