#include "altruct/structure/math/quadratic.h"
#include "altruct/structure/math/modulo.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

typedef moduloX<int> modx;
typedef quadratic<modx, 5> quad;
typedef quadraticX<modx> quadx;

template<>
modx quadratic_members<modx, 5, quadratic_storage::STATIC>::_D = modx(5, 1009);

template<typename Q>
vector<int> to_vec(const Q& q) {
    return{ q.a.v, q.a.M(), q.b.v, q.b.M(), q.D().v, q.D().M() };
}

TEST(quadratic_modx_test, constructor) {
    quad q1;
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 5, 1009}), to_vec(q1));
    quad q2(modx(10, 1009));
    EXPECT_EQ((vector<int>{10, 1009, 0, 1009, 5, 1009}), to_vec(q2));
    quad q3(modx(+2, 1009), modx(-5, 1009));
    EXPECT_EQ((vector<int>{2, 1009, 1004, 1009, 5, 1009}), to_vec(q3));
    quad q4(modx(+2, 1009), modx(-5, 1009), modx(7, 1009)); // ignore D if static
    EXPECT_EQ((vector<int>{2, 1009, 1004, 1009, 5, 1009}), to_vec(q4));
    quad q5(q4);
    EXPECT_EQ((vector<int>{2, 1009, 1004, 1009, 5, 1009}), to_vec(q5));
}

TEST(quadratic_modx_test, constructor_x) {
    quadx q1;
    EXPECT_EQ((vector<int>{0, 1, 0, 1, -1, 1}), to_vec(q1));
    quadx q2(modx(10, 1009));
    EXPECT_EQ((vector<int>{10, 1009, 0, 1009, -1, 1}), to_vec(q2));
    quadx q3(modx(+2, 1009), modx(-5, 1009));
    EXPECT_EQ((vector<int>{2, 1009, 1004, 1009, -1, 1}), to_vec(q3));
    quadx q4(modx(+2, 1009), modx(-5, 1009), modx(7, 1009)); // ignore D if static
    EXPECT_EQ((vector<int>{2, 1009, 1004, 1009, 7, 1009}), to_vec(q4));
    quadx q5(q4);
    EXPECT_EQ((vector<int>{2, 1009, 1004, 1009, 7, 1009}), to_vec(q5));
}

template<typename T>
void test_comparison(bool eq, bool lt, const T& lhs, const T& rhs) {
    ASSERT_FALSE(eq && lt);
    EXPECT_EQ(eq, lhs == rhs);
    EXPECT_EQ(!eq, lhs != rhs);
    EXPECT_EQ(lt, lhs < rhs);
    EXPECT_EQ(!(lt || eq), lhs > rhs);
    EXPECT_EQ((lt || eq), lhs <= rhs);
    EXPECT_EQ(!lt, lhs >= rhs);
}

TEST(quadratic_modx_test, operators_comparison) {
    test_comparison(true, false, quad(modx(2, 1009), modx(5, 1009)), quad(modx(2, 1009), modx(5, 1009)));
    test_comparison(false, false, quad(modx(2, 1009), modx(5, 1009)), quad(modx(2, 1009), modx(3, 1009)));
    test_comparison(false, true, quad(modx(2, 1009), modx(5, 1009)), quad(modx(2, 1009), modx(7, 1009)));
    test_comparison(false, true, quad(modx(2, 1009), modx(5, 1009)), quad(modx(4, 1009), modx(5, 1009)));
    test_comparison(false, true, quad(modx(2, 1009), modx(5, 1009)), quad(modx(4, 1009), modx(3, 1009)));
    test_comparison(false, true, quad(modx(2, 1009), modx(5, 1009)), quad(modx(4, 1009), modx(7, 1009)));
    test_comparison(false, false, quad(modx(2, 1009), modx(5, 1009)), quad(modx(1, 1009), modx(5, 1009)));
    test_comparison(false, false, quad(modx(2, 1009), modx(5, 1009)), quad(modx(1, 1009), modx(3, 1009)));
    test_comparison(false, false, quad(modx(2, 1009), modx(5, 1009)), quad(modx(1, 1009), modx(7, 1009)));
}

TEST(quadratic_modx_test, operators_comparison_x) {
    test_comparison(true, false, quadx(modx(2, 1009), modx(5, 1009), modx(5, 1009)), quadx(modx(2, 1009), modx(5, 1009), modx(5, 1009)));
    test_comparison(false, false, quadx(modx(2, 1009), modx(5, 1009), modx(5, 1009)), quadx(modx(2, 1009), modx(3, 1009), modx(5, 1009)));
    test_comparison(false, true, quadx(modx(2, 1009), modx(5, 1009), modx(5, 1009)), quadx(modx(2, 1009), modx(7, 1009), modx(5, 1009)));
    test_comparison(false, true, quadx(modx(2, 1009), modx(5, 1009), modx(5, 1009)), quadx(modx(4, 1009), modx(5, 1009), modx(5, 1009)));
    test_comparison(false, true, quadx(modx(2, 1009), modx(5, 1009), modx(5, 1009)), quadx(modx(4, 1009), modx(3, 1009), modx(5, 1009)));
    test_comparison(false, true, quadx(modx(2, 1009), modx(5, 1009), modx(5, 1009)), quadx(modx(4, 1009), modx(7, 1009), modx(5, 1009)));
    test_comparison(false, false, quadx(modx(2, 1009), modx(5, 1009), modx(5, 1009)), quadx(modx(1, 1009), modx(5, 1009), modx(5, 1009)));
    test_comparison(false, false, quadx(modx(2, 1009), modx(5, 1009), modx(5, 1009)), quadx(modx(1, 1009), modx(3, 1009), modx(5, 1009)));
    test_comparison(false, false, quadx(modx(2, 1009), modx(5, 1009), modx(5, 1009)), quadx(modx(1, 1009), modx(7, 1009), modx(5, 1009)));
}

TEST(quadratic_modx_test, operators_arithmetic) {
    const quad q1(modx(2, 1009), modx(-5, 1009));
    const quad q2(modx(3, 1009), modx(4, 1009));
    const quad q3(modx(3, 1009), modx(-2, 1009));
    EXPECT_EQ((vector<int>{5, 1009, 1008, 1009, 5, 1009}), to_vec(q1 + q2));
    EXPECT_EQ((vector<int>{1008, 1009, 1000, 1009, 5, 1009}), to_vec(q1 - q2));
    EXPECT_EQ((vector<int>{1007, 1009, 5, 1009, 5, 1009}), to_vec(-q1));
    EXPECT_EQ((vector<int>{915, 1009, 1002, 1009, 5, 1009}), to_vec(q1 * q2));
    EXPECT_EQ((vector<int>{4, 1009, 1, 1009, 5, 1009}), to_vec(q1 / q3));
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 5, 1009}), to_vec(q1 % q2));
    EXPECT_EQ((vector<int>{5, 1009, 1008, 1009, 5, 1009}), to_vec(q2 + q1));
    EXPECT_EQ((vector<int>{1, 1009, 9, 1009, 5, 1009}), to_vec(q2 - q1));
    EXPECT_EQ((vector<int>{1006, 1009, 1005, 1009, 5, 1009}), to_vec(-q2));
    EXPECT_EQ((vector<int>{915, 1009, 1002, 1009, 5, 1009}), to_vec(q2 * q1));
    EXPECT_EQ((vector<int>{1003, 1009, 15, 1009, 5, 1009}), to_vec(q1 * modx(-3, 1009)));
    EXPECT_EQ((vector<int>{1, 1009, 502, 1009, 5, 1009}), to_vec(q1 / modx(2, 1009)));
}

TEST(quadratic_modx_test, operators_arithmetic_x) {
    const quadx q1(modx(2, 1009), modx(-5, 1009), modx(5, 1009));
    const quadx q2(modx(3, 1009), modx(4, 1009), modx(5, 1009));
    const quadx q3(modx(3, 1009), modx(-2, 1009), modx(5, 1009));
    EXPECT_EQ((vector<int>{5, 1009, 1008, 1009, 5, 1009}), to_vec(q1 + q2));
    EXPECT_EQ((vector<int>{1008, 1009, 1000, 1009, 5, 1009}), to_vec(q1 - q2));
    EXPECT_EQ((vector<int>{1007, 1009, 5, 1009, 5, 1009}), to_vec(-q1));
    EXPECT_EQ((vector<int>{915, 1009, 1002, 1009, 5, 1009}), to_vec(q1 * q2));
    EXPECT_EQ((vector<int>{4, 1009, 1, 1009, 5, 1009}), to_vec(q1 / q3));
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 5, 1009}), to_vec(q1 % q2));
    EXPECT_EQ((vector<int>{5, 1009, 1008, 1009, 5, 1009}), to_vec(q2 + q1));
    EXPECT_EQ((vector<int>{1, 1009, 9, 1009, 5, 1009}), to_vec(q2 - q1));
    EXPECT_EQ((vector<int>{1006, 1009, 1005, 1009, 5, 1009}), to_vec(-q2));
    EXPECT_EQ((vector<int>{915, 1009, 1002, 1009, 5, 1009}), to_vec(q2 * q1));
    EXPECT_EQ((vector<int>{1003, 1009, 15, 1009, 5, 1009}), to_vec(q1 * modx(-3, 1009)));
    EXPECT_EQ((vector<int>{1, 1009, 502, 1009, 5, 1009}), to_vec(q1 / modx(2, 1009)));
}

TEST(quadratic_modx_test, operators_inplace) {
    const quad q1(modx(2, 1009), modx(-5, 1009));
    const quad q2(modx(3, 1009), modx(4, 1009));
    const quad q3(modx(3, 1009), modx(-2, 1009));
    quad qr;
    qr = q1; qr += q2;
    EXPECT_EQ((vector<int>{5, 1009, 1008, 1009, 5, 1009}), to_vec(qr));
    qr = q1; qr -= q2;
    EXPECT_EQ((vector<int>{1008, 1009, 1000, 1009, 5, 1009}), to_vec(qr));
    qr = q1; qr *= q2;
    EXPECT_EQ((vector<int>{915, 1009, 1002, 1009, 5, 1009}), to_vec(qr));
    qr = q1; qr /= q3;
    EXPECT_EQ((vector<int>{4, 1009, 1, 1009, 5, 1009}), to_vec(qr));
    qr = q1; qr %= q2;
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 5, 1009}), to_vec(qr));
    qr = q2; qr += q1;
    EXPECT_EQ((vector<int>{5, 1009, 1008, 1009, 5, 1009}), to_vec(qr));
    qr = q2; qr -= q1;
    EXPECT_EQ((vector<int>{1, 1009, 9, 1009, 5, 1009}), to_vec(qr));
    qr = q2; qr *= q1;
    EXPECT_EQ((vector<int>{915, 1009, 1002, 1009, 5, 1009}), to_vec(qr));
    qr = q1; qr *= modx(-3, 1009);
    EXPECT_EQ((vector<int>{1003, 1009, 15, 1009, 5, 1009}), to_vec(qr));
    qr = q1; qr /= modx(2, 1009);
    EXPECT_EQ((vector<int>{1, 1009, 502, 1009, 5, 1009}), to_vec(qr));
}

TEST(quadratic_modx_test, operators_inplace_x) {
    const quadx q1(modx(2, 1009), modx(-5, 1009), modx(5, 1009));
    const quadx q2(modx(3, 1009), modx(4, 1009), modx(5, 1009));
    const quadx q3(modx(3, 1009), modx(-2, 1009), modx(5, 1009));
    quadx qr;
    qr = q1; qr += q2;
    EXPECT_EQ((vector<int>{5, 1009, 1008, 1009, 5, 1009}), to_vec(qr));
    qr = q1; qr -= q2;
    EXPECT_EQ((vector<int>{1008, 1009, 1000, 1009, 5, 1009}), to_vec(qr));
    qr = q1; qr *= q2;
    EXPECT_EQ((vector<int>{915, 1009, 1002, 1009, 5, 1009}), to_vec(qr));
    qr = q1; qr /= q3;
    EXPECT_EQ((vector<int>{4, 1009, 1, 1009, 5, 1009}), to_vec(qr));
    qr = q1; qr %= q2;
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 5, 1009}), to_vec(qr));
    qr = q2; qr += q1;
    EXPECT_EQ((vector<int>{5, 1009, 1008, 1009, 5, 1009}), to_vec(qr));
    qr = q2; qr -= q1;
    EXPECT_EQ((vector<int>{1, 1009, 9, 1009, 5, 1009}), to_vec(qr));
    qr = q2; qr *= q1;
    EXPECT_EQ((vector<int>{915, 1009, 1002, 1009, 5, 1009}), to_vec(qr));
    qr = q1; qr *= modx(-3, 1009);
    EXPECT_EQ((vector<int>{1003, 1009, 15, 1009, 5, 1009}), to_vec(qr));
    qr = q1; qr /= modx(2, 1009);
    EXPECT_EQ((vector<int>{1, 1009, 502, 1009, 5, 1009}), to_vec(qr));
}

TEST(quadratic_modx_test, operators_inplace_self) {
    const quad q1(modx(2, 1009), modx(-5, 1009));
    quad qr;
    qr = q1; qr += qr;
    EXPECT_EQ((vector<int>{4, 1009, 999, 1009, 5, 1009}), to_vec(qr));
    qr = q1; qr -= qr;
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 5, 1009}), to_vec(qr));
    qr = q1; qr *= qr;
    EXPECT_EQ((vector<int>{129, 1009, 989, 1009, 5, 1009}), to_vec(qr));
    qr = q1; qr /= qr;
    EXPECT_EQ((vector<int>{1, 1009, 0, 1009, 5, 1009}), to_vec(qr));
    qr = q1; qr %= qr;
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 5, 1009}), to_vec(qr));
}

TEST(quadratic_modx_test, operators_inplace_self_x) {
    const quadx q1(modx(2, 1009), modx(-5, 1009), modx(5, 1009));
    quadx qr;
    qr = q1; qr += qr;
    EXPECT_EQ((vector<int>{4, 1009, 999, 1009, 5, 1009}), to_vec(qr));
    qr = q1; qr -= qr;
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 5, 1009}), to_vec(qr));
    qr = q1; qr *= qr;
    EXPECT_EQ((vector<int>{129, 1009, 989, 1009, 5, 1009}), to_vec(qr));
    qr = q1; qr /= qr;
    EXPECT_EQ((vector<int>{1, 1009, 0, 1009, 5, 1009}), to_vec(qr));
    qr = q1; qr %= qr;
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 5, 1009}), to_vec(qr));
}

TEST(quadratic_modx_test, conjugate) {
    const quad q1(modx(2, 1009), modx(-5, 1009));
    const quad q2(modx(2, 1009), modx(3, 1009));
    EXPECT_EQ((vector<int>{2, 1009, 5, 1009, 5, 1009}), to_vec(q1.conjugate()));
    EXPECT_EQ((vector<int>{2, 1009, 1006, 1009, 5, 1009}), to_vec(q2.conjugate()));
    EXPECT_EQ((vector<int>{2, 1009, 1006, 1009, 5, 1009}), to_vec(conjugateT<quad>::of(q2)));
}

TEST(quadratic_modx_test, conjugate_x) {
    const quadx q1(modx(2, 1009), modx(-5, 1009), modx(5, 1009));
    const quadx q2(modx(2, 1009), modx(3, 1009), modx(5, 1009));
    EXPECT_EQ((vector<int>{2, 1009, 5, 1009, 5, 1009}), to_vec(q1.conjugate()));
    EXPECT_EQ((vector<int>{2, 1009, 1006, 1009, 5, 1009}), to_vec(q2.conjugate()));
    EXPECT_EQ((vector<int>{2, 1009, 1006, 1009, 5, 1009}), to_vec(conjugateT<quadx>::of(q2)));
}

TEST(quadratic_modx_test, norm) {
    const quad q1(modx(2, 1009), modx(-5, 1009));
    const quad q2(modx(3, 1009), modx(4, 1009));
    EXPECT_EQ(modx(-121, 1009), q1.norm());
    EXPECT_EQ(modx(-71, 1009), q2.norm());
}

TEST(quadratic_modx_test, norm_x) {
    const quadx q1(modx(2, 1009), modx(-5, 1009), modx(5, 1009));
    const quadx q2(modx(3, 1009), modx(4, 1009), modx(5, 1009));
    EXPECT_EQ(modx(-121, 1009), q1.norm());
    EXPECT_EQ(modx(-71, 1009), q2.norm());
}

TEST(quadratic_modx_test, casts) {
    const quad q(modx(2, 1009), modx(-5, 1009));
    const quad e0 = zeroT<quad>::of(q);
    const quad e1 = identityT<quad>::of(q);
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 5, 1009}), to_vec(e0));
    EXPECT_EQ((vector<int>{1, 1009, 0, 1009, 5, 1009}), to_vec(e1));
    const quad q3 = castOf<quad>(3);
    EXPECT_EQ((vector<int>{3, 1, 0, 1, 5, 1009}), to_vec(q3));
    const quad q4 = castOf(q, 4);
    EXPECT_EQ((vector<int>{4, 1009, 0, 1009, 5, 1009}), to_vec(q4));
    const quad q6 = castOf(q, q4);
    EXPECT_EQ((vector<int>{4, 1009, 0, 1009, 5, 1009}), to_vec(q6));
    const quad q7 = castOf<quad>(q4);
    EXPECT_EQ((vector<int>{4, 1009, 0, 1009, 5, 1009}), to_vec(q7));
}

TEST(quadratic_modx_test, casts_x) {
    const quadx z(modx(2, 1009), modx(-5, 1009), modx(7, 1009));
    const quadx z0 = zeroT<quadx>::of(z);
    const quadx z1 = identityT<quadx>::of(z);
    EXPECT_EQ((vector<int>{0, 1009, 0, 1009, 7, 1009}), to_vec(z0));
    EXPECT_EQ((vector<int>{1, 1009, 0, 1009, 7, 1009}), to_vec(z1));
    const quadx z5 = castOf(z1, 5);
    EXPECT_EQ((vector<int>{5, 1009, 0, 1009, 7, 1009}), to_vec(z5));
    const quadx z6 = castOf(z1, z);
    EXPECT_EQ((vector<int>{2, 1009, 1004, 1009, 7, 1009}), to_vec(z6));
    const quadx z7 = castOf<quadx>(z6);
    EXPECT_EQ((vector<int>{2, 1009, 1004, 1009, 7, 1009}), to_vec(z7));
}
