#include "altruct/structure/math/permutation.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

typedef permutation<int> perm;

TEST(permutation_test, constructor) {
    perm p1;
    EXPECT_EQ(0, p1.n);
    EXPECT_EQ(perm::cycles_t{}, p1.cycles);

    perm p2(10);
    EXPECT_EQ(10, p2.n);
    EXPECT_EQ(perm::cycles_t{}, p2.cycles);

    perm p3(perm::cycles_t{ { 0, 7, 9 }, { 2, 8 }, { 3, 6, 4, 5 } }, 10);
    EXPECT_EQ(10, p3.n);
    EXPECT_EQ((perm::cycles_t{ { 0, 7, 9 }, { 2, 8 }, { 3, 6, 4, 5 } }), p3.cycles);

    perm p4(perm::line_t{ 2, 1, 3, 0, 4, 6, 5 });
    EXPECT_EQ(7, p4.n);
    EXPECT_EQ((perm::cycles_t{ { 0, 2, 3 }, { 5, 6 } }), p4.cycles);

    perm p5(perm::transpositions_t{ { 3, 5 }, { 3, 6 }, { 4, 2 } }, 10);
    EXPECT_EQ(10, p5.n);
    EXPECT_EQ((perm::cycles_t{ { 2, 4 }, { 3, 6, 5 } }), p5.cycles);
}

void test_permutation_comparison(const perm& p1, const perm& p2, bool eq, bool lt) {
    EXPECT_EQ(eq, p1 == p2);
    EXPECT_EQ(!eq, p1 != p2);
    EXPECT_EQ(lt, p1 < p2);
    EXPECT_EQ(!lt && !eq, p1 > p2);
    EXPECT_EQ(lt || eq, p1 <= p2);
    EXPECT_EQ(!lt, p1 >= p2);
}

TEST(permutation_test, operators_comparison) {
    const perm p1(perm::line_t{ 4, 2, 1, 6, 0, 5, 3 });
    const perm p2(perm::line_t{ 2, 1, 3, 0, 4, 6, 5 });
    const perm p3(perm::line_t{ 2, 1, 3, 0, 4, 6, 5, 7, 8 });
    test_permutation_comparison(p1, p1, true, false);
    test_permutation_comparison(p2, p2, true, false);
    test_permutation_comparison(p3, p3, true, false);
    test_permutation_comparison(p1, p2, false, false);
    test_permutation_comparison(p2, p1, false, true);
    test_permutation_comparison(p2, p3, false, true);
    test_permutation_comparison(p3, p2, false, false);
    test_permutation_comparison(p1, p3, false, true);
    test_permutation_comparison(p3, p1, false, false);
}

TEST(permutation_test, operators) {
    const perm p0;
    const perm p1(perm::line_t{ 4, 2, 1, 6, 0, 5, 3 });
    const perm p2(perm::line_t{ 2, 1, 3, 0, 4, 6, 5, 7, 8 });

    // multiplication by identity
    EXPECT_EQ(p0, p0 * p0);
    EXPECT_EQ(p1, p0 * p1);
    EXPECT_EQ(p1, p1 * p0);
    EXPECT_EQ(p2, p0 * p2);
    EXPECT_EQ(p2, p2 * p0);

    // multiplication
    EXPECT_EQ(perm(7), p1 * p1);
    EXPECT_EQ(perm(perm::cycles_t{ { 0, 4, 2, 1, 3, 5, 6 } }, 9), p1 * p2);
    EXPECT_EQ(perm(perm::cycles_t{ { 0, 1, 2, 6, 5, 3, 4 } }, 9), p2 * p1);
    EXPECT_EQ(perm(perm::cycles_t{ { 0, 3, 2 } }, 9), p2 * p2);

    // product works as composition
    auto l = perm::identity_line(10);
    random_shuffle(l.begin(), l.end());
    auto la = l; p2.apply_to(la); p1.apply_to(la);
    auto lb = l; (p1 * p2).apply_to(lb);
    EXPECT_EQ(la, lb);

    // division by identity
    EXPECT_EQ(p1.cycles, (p1 / p0).cycles);
    EXPECT_EQ(p2.cycles, (p2 / p0).cycles);

    // division by itself
    EXPECT_EQ(p0.cycles, (p0 / p0).cycles);
    EXPECT_EQ(p0.cycles, (p1 / p1).cycles);
    EXPECT_EQ(p0.cycles, (p2 / p2).cycles);

    // dvision
    auto p12 = p1 * p2;
    EXPECT_EQ(p1.cycles, (p12 / p2).cycles);
    auto p21 = p2 * p1;
    EXPECT_EQ(p2.cycles, (p21 / p1).cycles);

    // inplace
    perm pr;
    pr = p1; pr *= p2;
    EXPECT_EQ(p1 * p2, pr);
    pr = p1; pr *= pr;
    EXPECT_EQ(p1 * p1, pr);
    pr = p1; pr /= p2;
    EXPECT_EQ(p1 / p2, pr);
    pr = p1; pr /= pr;
    EXPECT_EQ(perm(7), pr);
}

TEST(permutation_test, power) {
    const perm p(perm::line_t{ 4, 2, 1, 6, 0, 5, 3 });
    const perm pi = p.inv();
    EXPECT_EQ(perm(perm::cycles_t{ { 0, 4 }, { 1, 2 }, { 3, 6 } }, 7), pi);
    EXPECT_EQ(perm(7), p * pi);
    EXPECT_EQ(perm(7), pi * p);
    EXPECT_EQ(perm(7), p.pow(0));
    EXPECT_EQ(p, p.pow(1));
    EXPECT_EQ(p * p, p.pow(2));
    EXPECT_EQ(p * p * p, p.pow(3));
    EXPECT_EQ(pi, p.pow(-1));
    EXPECT_EQ(pi * pi, p.pow(-2));
    EXPECT_EQ(pi * pi * pi, p.pow(-3));
}

perm split_to_cycles(const std::vector<int>& a, const std::vector<int>& lengths) {
    int i = 0;
    perm::cycles_t vc;
    for (auto l : lengths) {
        perm::cycle_t c;
        while (l-- > 0) c.push_back(a[i++]);
        vc.push_back(c);
    }
    return perm(vc, (int)a.size());
}

TEST(permutation_test, root) {
    const auto a = perm::identity_line(100);
    const perm p14 = split_to_cycles(a, { 2, 8, 20, 9, 49 }).pow(14);
    auto p14_14 = p14.root(14);
    EXPECT_EQ(p14.to_line(), p14_14.pow(14).to_line());
    auto p14_7 = p14.root(7);
    EXPECT_EQ(p14.to_line(), p14_7.pow(7).to_line());
    auto p14_2 = p14.root(2);
    EXPECT_EQ(p14.to_line(), p14_2.pow(2).to_line());

    auto p14_14_0 = p14.root(14, 0);
    EXPECT_EQ(0, p14_14_0.to_transpositions().size() % 2);
    EXPECT_EQ(p14.to_line(), p14_14_0.pow(14).to_line());
    auto p14_7_0 = p14.root(7, 0);
    EXPECT_EQ(0, p14_7_0.to_transpositions().size() % 2);
    EXPECT_EQ(p14.to_line(), p14_7_0.pow(7).to_line());
    auto p14_2_0 = p14.root(2, 0);
    EXPECT_EQ(0, p14_2_0.to_transpositions().size() % 2);
    EXPECT_EQ(p14.to_line(), p14_2_0.pow(2).to_line());

    auto p14_14_1 = p14.root(14, 1);
    EXPECT_EQ(1, p14_14_1.to_transpositions().size() % 2);
    EXPECT_EQ(p14.to_line(), p14_14_1.pow(14).to_line());
    auto p14_7_1 = p14.root(7, 1);
    EXPECT_EQ(perm(), p14_7_1) << "there should be no such root";
    auto p14_2_1 = p14.root(2, 1);
    EXPECT_EQ(1, p14_2_1.to_transpositions().size() % 2);
    EXPECT_EQ(p14.to_line(), p14_2_1.pow(2).to_line());
}

TEST(permutation_test, conversion) {
    const perm p(perm::cycles_t{{ 0, 2, 3 }, { 5, 6 }}, 7);
    EXPECT_EQ((perm::cycles_t{ { 0, 2, 3 }, { 5, 6 } }), p.to_cycles());
    EXPECT_EQ((perm::cycles_t{ { 0, 2, 3 }, { 5, 6 }, { 1 }, { 4 } }), p.to_all_cycles());
    EXPECT_EQ((perm::line_t{ 2, 1, 3, 0, 4, 6, 5 }), p.to_line());
    EXPECT_EQ((perm::transpositions_t{ { 0, 2 }, { 2, 3 }, { 5, 6 } }), p.to_transpositions());
}

TEST(permutation_test, static_helpers) {
    const auto cycles = perm::cycles_t{ { 0, 2, 3 }, { 5, 6 } };
    const auto line = perm::line_t{ 2, 1, 3, 0, 4, 6, 5 };
    const auto line10 = perm::line_t{ 2, 1, 3, 0, 4, 6, 5, 7, 8, 9 };
    const auto transpositions = perm::transpositions_t{ { 0, 2 }, { 2, 3 }, { 5, 6 } };

    // static cycle helpers
    EXPECT_EQ((perm::cycles_t{ { 0 }, { 1 }, { 2 }, { 3 }, { 4 } }), perm::all_cycles(perm::cycles_t(), 5));
    EXPECT_EQ((perm::cycles_t{ { 0, 2, 3 }, { 5, 6 }, { 1 }, { 4 }, { 7 }, { 8 }, { 9 } }), perm::all_cycles(cycles, 10));

    // static line helpers
    EXPECT_EQ((perm::line_t{ 0, 1, 2, 3, 4, 5, 6 }), perm::identity_line(7));
    auto temp_line = line;
    EXPECT_EQ((perm::line_t{ 2, 1, 3, 0, 4, 6, 5, 7, 8, 9 }), perm::expand_line(temp_line, 10));

    // static conversions
    EXPECT_EQ(line, perm::cycles_to_line(cycles, 7));
    EXPECT_EQ(line10, perm::cycles_to_line(cycles, 10));
    EXPECT_EQ(cycles, perm::line_to_cycles(line));
    EXPECT_EQ(line, perm::transpositions_to_line(transpositions, 7));
    EXPECT_EQ(line10, perm::transpositions_to_line(transpositions, 10));
    EXPECT_EQ(transpositions, perm::line_to_transpositions(line));
    EXPECT_EQ(cycles, perm::transpositions_to_cycles(transpositions));
    EXPECT_EQ(transpositions, perm::cycles_to_transpositions(cycles));

    // static size helpers
    EXPECT_EQ(0, perm::size(perm::cycles_t()));
    EXPECT_EQ(7, perm::size(cycles));
    EXPECT_EQ(0, perm::size(perm::transpositions_t()));
    EXPECT_EQ(7, perm::size(transpositions));
}

TEST(permutation_test, identity) {
    const perm p(perm::line_t{ 4, 2, 1, 6, 0, 5, 3 });
    EXPECT_EQ(perm(7), identityT<perm>::of(p));
}
