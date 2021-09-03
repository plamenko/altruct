#include "altruct/structure/math/clifford3.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

template<typename T>
void test_comparison(bool eq, bool lt, const T& lhs, const T& rhs) {
    ASSERT_FALSE(eq && lt);
    EXPECT_EQ(eq, lhs == rhs);
    EXPECT_EQ(!eq, lhs != rhs);
    EXPECT_EQ(lt, lhs < rhs);
    EXPECT_EQ(!(lt || eq), lhs > rhs);
    EXPECT_EQ((lt || eq), lhs <= rhs);
    EXPECT_EQ(!lt, lhs >= rhs);
    EXPECT_EQ(eq, rhs == lhs);
    EXPECT_EQ(!eq, rhs != lhs);
    EXPECT_EQ(lt, rhs > lhs);
    EXPECT_EQ(!(lt || eq), rhs < lhs);
    EXPECT_EQ((lt || eq), rhs >= lhs);
    EXPECT_EQ(!lt, rhs <= lhs);
}

using rot3 = cl3::rotor<double>;
using vec3 = cl3::vector<double>;
using mvec3 = cl3::multivector<double>;

// arithmetic operator coverage:
//
// "#": "assumed";  "/": "ill-defined";  "!": "unnecessary, use rotor(S)";  "+": "handled"; 
//
// * S R V M   / S R V M   + S R V M   - S R V M   *=S R V M   /=S R V M   +=S R V M   -=S R V M
// S # + + +   S # + + +   S # ! ! !   S # ! ! !   S # / / /   S # / / /   S # / / /   S # / / /
// R + + + +   R + + + +   R + + + +   R + + + +   R + + / /   R + + / /   R + + / /   R + + / /
// V + + + +   V + + + +   V ! + + +   V ! + + +   V + + / /   V + + / /   V / / + /   V / / + /
// M + + + +   M + + + +   M + + + +   M + + + +   M + + + +   M + + + +   M + + + +   M + + + +

//------------------------------------------------------------------------------------------------/
TEST(clifford3_test, rotor_constructor) {
    rot3 r1;
    EXPECT_EQ(1.0, r1.s);
    EXPECT_EQ(0.0, r1.yz);
    EXPECT_EQ(0.0, r1.zx);
    EXPECT_EQ(0.0, r1.xy);
    rot3 r2(5.0);
    EXPECT_EQ(5.0, r2.s);
    EXPECT_EQ(0.0, r2.yz);
    EXPECT_EQ(0.0, r2.zx);
    EXPECT_EQ(0.0, r2.xy);
    rot3 r3(5.0, 4.0, 3.0, 2.0);
    EXPECT_EQ(5.0, r3.s);
    EXPECT_EQ(4.0, r3.yz);
    EXPECT_EQ(3.0, r3.zx);
    EXPECT_EQ(2.0, r3.xy);
    rot3 r4(r3);
    EXPECT_EQ(5.0, r4.s);
    EXPECT_EQ(4.0, r4.yz);
    EXPECT_EQ(3.0, r4.zx);
    EXPECT_EQ(2.0, r4.xy);
}

TEST(clifford3_test, rotor_operators_comparison) {
    test_comparison(true, false, rot3(2, 5, 7, 8), rot3(2, 5, 7, 8));
    test_comparison(false, false, rot3(5, 1, 2, 3), rot3(4, 7, 8, 9));
    test_comparison(false, false, rot3(5, 4, 1, 2), rot3(5, 3, 8, 9));
    test_comparison(false, false, rot3(5, 4, 3, 1), rot3(5, 4, 2, 9));
    test_comparison(false, false, rot3(5, 4, 3, 2), rot3(5, 4, 3, 1));
    test_comparison(true, false, rot3(5, 4, 3, 2), rot3(5, 4, 3, 2));
}

void ROT_EXPECT_NEAR(const rot3& r1, const rot3& r2, double eps = 1e-10) {
    EXPECT_NEAR(r1.s, r2.s, eps);
    EXPECT_NEAR(r1.yz, r2.yz, eps);
    EXPECT_NEAR(r1.zx, r2.zx, eps);
    EXPECT_NEAR(r1.xy, r2.xy, eps);
}

TEST(clifford3_test, rotor_operators_arithmetic) {
    ROT_EXPECT_NEAR(rot3(30, 40, 50, 100), 10. * rot3(3, 4, 5, 10));
    ROT_EXPECT_NEAR(rot3(30, -40, -50, -100), 1500. / rot3(3, 4, 5, 10));
    ROT_EXPECT_NEAR(rot3(30, 40, 50, 100), rot3(3, 4, 5, 10) * 10.);
    ROT_EXPECT_NEAR(rot3(3, 4, 5, 10), rot3(30, 40, 50, 100) / 10.);
    ROT_EXPECT_NEAR(rot3(-34, 23, 86, 63), rot3(7, 5, 3, 2) * rot3(3, 4, 5, 10));
    ROT_EXPECT_NEAR(rot3(76, 7, -68, -51), rot3(1050, 750, 450, 300) / rot3(3, 4, 5, 10));

    ROT_EXPECT_NEAR(rot3(13, 4, 5, 10), rot3(3, 4, 5, 10) + 10.);
    ROT_EXPECT_NEAR(rot3(-7, 4, 5, 10), rot3(3, 4, 5, 10) - 10.);
    ROT_EXPECT_NEAR(rot3(10, 9, 8, 12), rot3(7, 5, 3, 2) + rot3(3, 4, 5, 10));
    ROT_EXPECT_NEAR(rot3(4, 1, -2, -8), rot3(7, 5, 3, 2) - rot3(3, 4, 5, 10));
    ROT_EXPECT_NEAR(rot3(-3, 4, -5, 10), -rot3(3, -4, 5, -10));
}

TEST(clifford3_test, rotor_operators_inplace) {
    rot3 r;
    r = rot3(3, 4, 5, 10);
    ROT_EXPECT_NEAR(rot3(30, 40, 50, 100), r *= 10.);
    r = rot3(30, 40, 50, 100);
    ROT_EXPECT_NEAR(rot3(3, 4, 5, 10), r /= 10.);
    r = rot3(7, 5, 3, 2);
    ROT_EXPECT_NEAR(rot3(-34, 23, 86, 63), r *= rot3(3, 4, 5, 10));
    r = rot3(1050, 750, 450, 300);
    ROT_EXPECT_NEAR(rot3(76, 7, -68, -51), r /= rot3(3, 4, 5, 10));

    r = rot3(3, 4, 5, 10);
    ROT_EXPECT_NEAR(rot3(13, 4, 5, 10), r += 10.);
    r = rot3(3, 4, 5, 10);
    ROT_EXPECT_NEAR(rot3(-7, 4, 5, 10), r -= 10.);
    r = rot3(7, 5, 3, 2);
    ROT_EXPECT_NEAR(rot3(10, 9, 8, 12), r += rot3(3, 4, 5, 10));
    r = rot3(7, 5, 3, 2);
    ROT_EXPECT_NEAR(rot3(4, 1, -2, -8), r -= rot3(3, 4, 5, 10));
}

TEST(clifford3_test, rotor_operators_inplace_self) {
    rot3 r;
    r = rot3(3, 4, 5, 10);
    ROT_EXPECT_NEAR(rot3(-132, 24, 30, 60), r *= r);
    r = rot3(3, 4, 5, 10);
    ROT_EXPECT_NEAR(rot3(1, 0, 0, 0), r /= r);
    r = rot3(3, 4, 5, 10);
    ROT_EXPECT_NEAR(rot3(6, 8, 10, 20), r += r);
    r = rot3(3, 4, 5, 10);
    ROT_EXPECT_NEAR(rot3(0, 0, 0, 0), r -= r);
}

TEST(clifford3_test, rotor_functions) {
    const rot3 r = rot3(3, 4, 5, 10);
    ROT_EXPECT_NEAR(rot3(3, -4, -5, -10), r.conj());
    EXPECT_NEAR(150., r.abs2(), 1e-10);
    ROT_EXPECT_NEAR(rot3(3, -4, -5, -10) / 150., r.inv());
}

TEST(clifford3_test, rotor_specializations) {
    const rot3 r = rot3(3, 4, 5, 10);
    ROT_EXPECT_NEAR(rot3(12, 0, 0, 0), castOf<rot3>(12));
    ROT_EXPECT_NEAR(rot3(12, 0, 0, 0), castOf<rot3>(r, 12));
    ROT_EXPECT_NEAR(rot3(0, 0, 0, 0), zeroOf(r));
    ROT_EXPECT_NEAR(rot3(1, 0, 0, 0), identityOf(r));
    ROT_EXPECT_NEAR(rot3(3, -4, -5, -10), conjugateT<rot3>::of(r));
}

//------------------------------------------------------------------------------------------------/
TEST(clifford3_test, vector_constructor) {
    vec3 v1;
    EXPECT_EQ(0.0, v1.x);
    EXPECT_EQ(0.0, v1.y);
    EXPECT_EQ(0.0, v1.z);
    EXPECT_EQ(0.0, v1.w);
    vec3 v2(5.0, 4.0, 3.0);
    EXPECT_EQ(5.0, v2.x);
    EXPECT_EQ(4.0, v2.y);
    EXPECT_EQ(3.0, v2.z);
    EXPECT_EQ(0.0, v2.w);
    vec3 v3(5.0, 4.0, 3.0, 2.0);
    EXPECT_EQ(5.0, v3.x);
    EXPECT_EQ(4.0, v3.y);
    EXPECT_EQ(3.0, v3.z);
    EXPECT_EQ(2.0, v3.w);
    vec3 v4(v3);
    EXPECT_EQ(5.0, v4.x);
    EXPECT_EQ(4.0, v4.y);
    EXPECT_EQ(3.0, v4.z);
    EXPECT_EQ(2.0, v4.w);
}

TEST(clifford3_test, vector_operators_comparison) {
    test_comparison(true, false, vec3(2, 5, 7, 8), vec3(2, 5, 7, 8));
    test_comparison(false, false, vec3(5, 1, 2, 3), vec3(4, 7, 8, 9));
    test_comparison(false, false, vec3(5, 4, 1, 2), vec3(5, 3, 8, 9));
    test_comparison(false, false, vec3(5, 4, 3, 1), vec3(5, 4, 2, 9));
    test_comparison(false, false, vec3(5, 4, 3, 2), vec3(5, 4, 3, 1));
    test_comparison(true, false, vec3(5, 4, 3, 2), vec3(5, 4, 3, 2));
}

void VEC_EXPECT_NEAR(const vec3& v1, const vec3& v2, double eps = 1e-10) {
    EXPECT_NEAR(v1.x, v2.x, eps);
    EXPECT_NEAR(v1.y, v2.y, eps);
    EXPECT_NEAR(v1.z, v2.z, eps);
    EXPECT_NEAR(v1.w, v2.w, eps);
}

TEST(clifford3_test, vector_operators_arithmetic) {
    ROT_EXPECT_NEAR(rot3(36, 89, 32, 53), vec3(7, 5, 3, 2) * vec3(3, 4, 5, 10));
    ROT_EXPECT_NEAR(rot3(76, -51, -68, -7), vec3(1050, 750, 450, 300) / vec3(3, 4, 5, 10));

    VEC_EXPECT_NEAR(vec3(30, 40, 50, 100), 10. * vec3(3, 4, 5, 10));
    VEC_EXPECT_NEAR(vec3(30, 40, 50, -100), 1500. / vec3(3, 4, 5, 10));
    VEC_EXPECT_NEAR(vec3(30, 40, 50, 100), vec3(3, 4, 5, 10) * 10.);
    VEC_EXPECT_NEAR(vec3(3, 4, 5, 10), vec3(30, 40, 50, 100) / 10.);
    VEC_EXPECT_NEAR(vec3(-22, 63, -26, 89), vec3(7, 5, 3, 2) * rot3(3, 4, 5, 10));
    VEC_EXPECT_NEAR(vec3(64, -33, 44, -77), vec3(1050, 750, 450, 300) / rot3(3, 4, 5, 10));
    VEC_EXPECT_NEAR(vec3(-36, 17, 4, 107), rot3(7, 5, 3, 2) * vec3(3, 4, 5, 10));
    VEC_EXPECT_NEAR(vec3(64, 77, 44, -33), rot3(1050, 750, 450, 300) / vec3(3, 4, 5, 10));

    VEC_EXPECT_NEAR(vec3(10, 9, 8, 12), vec3(7, 5, 3, 2) + vec3(3, 4, 5, 10));
    VEC_EXPECT_NEAR(vec3(4, 1, -2, -8), vec3(7, 5, 3, 2) - vec3(3, 4, 5, 10));
    VEC_EXPECT_NEAR(vec3(-3, 4, -5, 10), -vec3(3, -4, 5, -10));
}

TEST(clifford3_test, vector_operators_inplace) {
    vec3 v;
    v = vec3(3, 4, 5, 10);
    VEC_EXPECT_NEAR(vec3(30, 40, 50, 100), v *= 10.);
    v = vec3(30, 40, 50, 100);
    VEC_EXPECT_NEAR(vec3(3, 4, 5, 10), v /= 10.);
    v = vec3(7, 5, 3, 2);
    VEC_EXPECT_NEAR(vec3(-22, 63, -26, 89), v *= rot3(3, 4, 5, 10));
    v = vec3(1050, 750, 450, 300);
    VEC_EXPECT_NEAR(vec3(64, -33, 44, -77), v /= rot3(3, 4, 5, 10));

    v = vec3(7, 5, 3, 2);
    VEC_EXPECT_NEAR(vec3(10, 9, 8, 12), v += vec3(3, 4, 5, 10));
    v = vec3(7, 5, 3, 2);
    VEC_EXPECT_NEAR(vec3(4, 1, -2, -8), v -= vec3(3, 4, 5, 10));
}

TEST(clifford3_test, vector_operators_inplace_self) {
    vec3 v;
    v = vec3(3, 4, 5, 10);
    VEC_EXPECT_NEAR(vec3(6, 8, 10, 20), v += v);
    v = vec3(3, 4, 5, 10);
    VEC_EXPECT_NEAR(vec3(0, 0, 0, 0), v -= v);
}

TEST(clifford3_test, vector_functions) {
    const vec3 v = vec3(3, 4, 5, 10);
    VEC_EXPECT_NEAR(vec3(3, 4, 5, -10), v.conj());
    EXPECT_NEAR(150., v.abs2(), 1e-10);
    VEC_EXPECT_NEAR(vec3(3, 4, 5, -10) / 150., v.inv());
}

TEST(clifford3_test, vector_specializations) {
    const vec3 v = vec3(3, 4, 5, 10);
    VEC_EXPECT_NEAR(vec3(0, 0, 0, 12), castOf<vec3>(12));
    VEC_EXPECT_NEAR(vec3(0, 0, 0, 12), castOf<vec3>(v, 12));
    VEC_EXPECT_NEAR(vec3(0, 0, 0, 0), zeroOf(v));
    VEC_EXPECT_NEAR(vec3(0, 0, 0, 0), identityOf(v)); // there is no identity of type vector
    VEC_EXPECT_NEAR(vec3(3, 4, 5, -10), conjugateT<vec3>::of(v));
}

//------------------------------------------------------------------------------------------------/
TEST(clifford3_test, multivector_constructor) {
    mvec3 m1;
    EXPECT_EQ(rot3(1, 0, 0, 0), m1.r);
    EXPECT_EQ(vec3(0, 0, 0, 0), m1.v);
    mvec3 m2(rot3(2, 3, 5, 7));
    EXPECT_EQ(rot3(2, 3, 5, 7), m2.r);
    EXPECT_EQ(vec3(0, 0, 0, 0), m2.v);
    mvec3 m3(vec3(3, 4, 5, 10));
    EXPECT_EQ(rot3(0, 0, 0, 0), m3.r);
    EXPECT_EQ(vec3(3, 4, 5, 10), m3.v);
    mvec3 m4(m3);
    EXPECT_EQ(rot3(0, 0, 0, 0), m4.r);
    EXPECT_EQ(vec3(3, 4, 5, 10), m4.v);
}

TEST(clifford3_test, multivector_operators_comparison) {
    const rot3 r1(3, 4, 5, 10);
    const vec3 v1(2, 5, 7, 8);
    const rot3 r2(7, 5, 3, 2);
    const vec3 v2(4, 1, 9, 3);
    test_comparison(true, false, mvec3(r1, v1), mvec3(r1, v1));
    test_comparison(false, false, mvec3(r2, v1), mvec3(r1, v2));
    test_comparison(false, false, mvec3(r1, v2), mvec3(r1, v1));
}

void MVEC_EXPECT_NEAR(const mvec3& m1, const mvec3& m2, double eps = 1e-10) {
    ROT_EXPECT_NEAR(m1.r, m2.r, eps);
    VEC_EXPECT_NEAR(m1.v, m2.v, eps);
}

TEST(clifford3_test, multivector_operators_arithmetic) {
    const rot3 r1(3, 4, 5, 10);
    const vec3 v1(2, 5, 7, 8);
    const rot3 r2(7, 5, 3, 2);
    const vec3 v2(4, 1, 9, 3);
    const mvec3 m1(r1, v1);
    const mvec3 m2(r2, v2);

    MVEC_EXPECT_NEAR(mvec3(r1, v1), r1 + v1);
    MVEC_EXPECT_NEAR(mvec3(r1, -v1), r1 - v1);
    MVEC_EXPECT_NEAR(mvec3(r1, v1), v1 + r1);
    MVEC_EXPECT_NEAR(mvec3(-r1, v1), v1 - r1);
    
    MVEC_EXPECT_NEAR(mvec3(10. * r2, 10. * v2), 10. * m2);
    MVEC_EXPECT_NEAR(mvec3(rot3(-13, -3, 1, -16), vec3(14, 7, 13, 11)), 100. / m2);
    MVEC_EXPECT_NEAR(mvec3(r1 * 10., v1 * 10), m1 * 10.);
    MVEC_EXPECT_NEAR(mvec3(r1 / 10., v1 / 10), m1 / 10.);
    MVEC_EXPECT_NEAR(mvec3(r1 * r2, r1 * v2), r1 * m2);
    MVEC_EXPECT_NEAR(mvec3(r1 * m2.inv()), r1 / m2);
    MVEC_EXPECT_NEAR(mvec3(r1 * r2, v1 * r2), m1 * r2);
    MVEC_EXPECT_NEAR(mvec3(m1 * r2.inv()), m1 / r2);
    MVEC_EXPECT_NEAR(mvec3(v1 * v2, v1 * r2), v1 * m2);
    MVEC_EXPECT_NEAR(mvec3(v1 * m2.inv()), v1 / m2);
    MVEC_EXPECT_NEAR(mvec3(v1 * v2, r1 * v2), m1 * v2);
    MVEC_EXPECT_NEAR(mvec3(m1 * v2.inv()), m1 / v2);
    MVEC_EXPECT_NEAR(mvec3(r1 * r2 + v1 * v2, v1 * r2 + r1 * v2), m1 * m2);
    MVEC_EXPECT_NEAR(mvec3(m1 * m2.inv()), m1 / m2);

    MVEC_EXPECT_NEAR(mvec3(r1 + 10., v1), m1 + 10.);
    MVEC_EXPECT_NEAR(mvec3(r1 - 10., v1), m1 - 10.);
    MVEC_EXPECT_NEAR(mvec3(r1 + r2, v2), r1 + m2);
    MVEC_EXPECT_NEAR(mvec3(r1 - r2, -v2), r1 - m2);
    MVEC_EXPECT_NEAR(mvec3(r1 + r2, v1), m1 + r2);
    MVEC_EXPECT_NEAR(mvec3(r1 - r2, v1), m1 - r2);
    MVEC_EXPECT_NEAR(mvec3(r2, v1 + v2), v1 + m2);
    MVEC_EXPECT_NEAR(mvec3(-r2, v1 - v2), v1 - m2);
    MVEC_EXPECT_NEAR(mvec3(r1, v1 + v2), m1 + v2);
    MVEC_EXPECT_NEAR(mvec3(r1, v1 - v2), m1 - v2);
    MVEC_EXPECT_NEAR(mvec3(r1 + r2, v1 + v2), m1 + m2);
    MVEC_EXPECT_NEAR(mvec3(r1 - r2, v1 - v2), m1 - m2);
    MVEC_EXPECT_NEAR(mvec3(-r2, -v2), -m2);
}

TEST(clifford3_test, multivector_operators_inplace) {
    const rot3 r1(3, 4, 5, 10);
    const vec3 v1(2, 5, 7, 8);
    const rot3 r2(7, 5, 3, 2);
    const vec3 v2(4, 1, 9, 3);
    const mvec3 m1(r1, v1);
    const mvec3 m2(r2, v2);

    mvec3 m;
    m = m1;
    MVEC_EXPECT_NEAR(mvec3(r1 * 10., v1 * 10.), m *= 10.);
    m = m1;
    MVEC_EXPECT_NEAR(mvec3(r1 / 10., v1 / 10.), m /= 10.);
    m = m1;
    MVEC_EXPECT_NEAR(mvec3(r1 * r2, v1 * r2), m *= r2);
    m = m1;
    MVEC_EXPECT_NEAR(mvec3(m1 * r2.inv()), m /= r2);
    m = m1;
    MVEC_EXPECT_NEAR(mvec3(v1 * v2, r1 * v2), m *= v2);
    m = m1;
    MVEC_EXPECT_NEAR(mvec3(m1 * v2.inv()), m /= v2);
    m = m1;
    MVEC_EXPECT_NEAR(mvec3(r1 * r2 + v1 * v2, v1 * r2 + r1 * v2), m *= m2);
    m = m1;
    MVEC_EXPECT_NEAR(mvec3(m1 * m2.inv()), m /= m2);

    m = m1;
    MVEC_EXPECT_NEAR(mvec3(r1 + 10., v1), m += 10.);
    m = m1;
    MVEC_EXPECT_NEAR(mvec3(r1 - 10., v1), m -= 10.);
    m = m1;
    MVEC_EXPECT_NEAR(mvec3(r1 + r2, v1), m += r2);
    m = m1;
    MVEC_EXPECT_NEAR(mvec3(r1 - r2, v1), m -= r2);
    m = m1;
    MVEC_EXPECT_NEAR(mvec3(r1, v1 + v2), m += v2);
    m = m1;
    MVEC_EXPECT_NEAR(mvec3(r1, v1 - v2), m -= v2);
    m = m1;
    MVEC_EXPECT_NEAR(mvec3(r1 + r2, v1 + v2), m += m2);
    m = m1;
    MVEC_EXPECT_NEAR(mvec3(r1 - r2, v1 - v2), m -= m2);
}

TEST(clifford3_test, multivector_operators_inplace_self) {
    const rot3 r1(3, 4, 5, 10);
    const vec3 v1(2, 5, 7, 8);
    const mvec3 m1(r1, v1);
    mvec3 m;
    m = m1;
    MVEC_EXPECT_NEAR(mvec3(r1 * r1 + v1 * v1, r1 * v1 + v1 * r1), m *= m);
    m = m1;
    MVEC_EXPECT_NEAR(mvec3(1), m /= m);
    m = m1;
    MVEC_EXPECT_NEAR(mvec3(r1 + r1, v1 + v1), m += m);
    m = m1;
    MVEC_EXPECT_NEAR(mvec3(0), m -= m);
}

TEST(clifford3_test, multivector_functions) {
    const mvec3 m(rot3(7, 5, 3, 2), vec3(4, 1, 9, 3));
    MVEC_EXPECT_NEAR(mvec3(rot3(7, -5, -3, -2), vec3(4, 1, 9, -3)), m.conj());
    MVEC_EXPECT_NEAR(mvec3(rot3(-13, -3, 1, -16), vec3(14, 7, 13, 11)) / 100., m.inv());
}

TEST(clifford3_test, multivector_specializations) {
    const mvec3 m(rot3(7, 5, 3, 2), vec3(4, 1, 9, 3));
    MVEC_EXPECT_NEAR(mvec3(12), castOf<mvec3>(12));
    MVEC_EXPECT_NEAR(mvec3(12), castOf<mvec3>(m, 12));
    MVEC_EXPECT_NEAR(mvec3(0), zeroOf(m));
    MVEC_EXPECT_NEAR(mvec3(1), identityOf(m));
    MVEC_EXPECT_NEAR(mvec3(rot3(7, -5, -3, -2), vec3(4, 1, 9, -3)), conjugateT<mvec3>::of(m));
}
