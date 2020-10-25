#include "altruct/algorithm/math/fractions.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

typedef fraction<int> frac;

TEST(fractions_test, farey_sequence_prev_inc) {
    vector<frac> vf;
    int n = 8;
    frac qp = frac(-1, 0), q = frac(0, 1), qn;
    for (; q < frac(1, 1); qp = q, q = qn) {
        vf.push_back(q);
        qn = farey_neighbour(n, qp, q);
    }
    vf.push_back(q);
    EXPECT_EQ((vector<frac>{
        { 0, 1 }, { 1, 8 }, { 1, 7 }, { 1, 6 }, { 1, 5 }, { 1, 4 }, { 2, 7 }, { 1, 3 }, { 3, 8 }, { 2, 5 }, { 3, 7 }, { 1, 2 },
        { 4, 7 }, { 3, 5 }, { 5, 8 }, { 2, 3 }, { 5, 7 }, { 3, 4 }, { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 8 }, { 1, 1 }}), vf);
}

TEST(fractions_test, farey_sequence_prev_dec) {
    vector<frac> vf;
    int n = 5;
    frac qp = frac(+1, 0), q = frac(1, 1), qn;
    for (; q > frac(0, 1); qp = q, q = qn) {
        vf.push_back(q);
        qn = farey_neighbour(n, qp, q);
    }
    vf.push_back(q);
    EXPECT_EQ((vector<frac>{{1, 1}, { 4, 5 }, { 3, 4 }, { 2, 3 }, { 3, 5 }, { 1, 2 }, { 2, 5 }, { 1, 3 }, { 1, 4 }, { 1, 5 }, { 0, 1 }}), vf);
}

TEST(fractions_test, farey_sequence_inc) {
    vector<frac> vf;
    for (frac q = frac(0, 1); q <= frac(1, 1); q = farey_neighbour<int>(5, { -1, 0 }, q)) {
        vf.push_back(q);
    }
    EXPECT_EQ((vector<frac>{{ 0, 1 }, { 1, 5 }, { 1, 4 }, { 1, 3 }, { 2, 5 }, { 1, 2 }, { 3, 5 }, { 2, 3 }, { 3, 4 }, { 4, 5 }, { 1, 1 }}), vf);
}

TEST(fractions_test, farey_sequence_dec) {
    vector<frac> vf;
    for (frac q = frac(1, 1); q >= frac(0, 1); q = farey_neighbour<int>(5, { +1, 0 }, q)) {
        vf.push_back(q);
    }
    EXPECT_EQ((vector<frac>{{ 1, 1 }, { 4, 5 }, { 3, 4 }, { 2, 3 }, { 3, 5 }, { 1, 2 }, { 2, 5 }, { 1, 3 }, { 1, 4 }, { 1, 5 }, { 0, 1 }}), vf);
}

TEST(fractions_test, repeating_decimal) {
    EXPECT_EQ(make_pair(std::vector<int>{5}, size_t(0)), repeating_decimal(10, 1, 2));
    EXPECT_EQ(make_pair(std::vector<int>{3}, size_t(1)), repeating_decimal(10, 1, 3));
    EXPECT_EQ(make_pair(std::vector<int>{2, 5}, size_t(0)), repeating_decimal(10, 1, 4));
    EXPECT_EQ(make_pair(std::vector<int>{0, 1}, size_t(0)), repeating_decimal(10, 1, 100));
    EXPECT_EQ(make_pair(std::vector<int>{1, 4, 2, 8, 5, 7}, size_t(6)), repeating_decimal(10, 1, 7));
    EXPECT_EQ(make_pair(std::vector<int>{8, 3, 1, 4, 2, 8, 5, 7}, size_t(6)), repeating_decimal(10, 291, 350));
}

TEST(fractions_test, rational_digit) {
    EXPECT_EQ(2, rational_digit(0, 10, 1, 4));
    EXPECT_EQ(5, rational_digit(1, 10, 1, 4));
    EXPECT_EQ(0, rational_digit(2, 10, 1, 4));
    EXPECT_EQ(0, rational_digit(3, 10, 1, 4));
    EXPECT_EQ(8, rational_digit(0, 10, 291, 350));
    EXPECT_EQ(3, rational_digit(1, 10, 291, 350));
    EXPECT_EQ(1, rational_digit(2, 10, 291, 350));
    EXPECT_EQ(4, rational_digit(3, 10, 291, 350));
    EXPECT_EQ(2, rational_digit(4, 10, 291, 350));
    EXPECT_EQ(8, rational_digit(5, 10, 291, 350));
    EXPECT_EQ(5, rational_digit(6, 10, 291, 350));
    EXPECT_EQ(7, rational_digit(7, 10, 291, 350));
    EXPECT_EQ(1, rational_digit(8, 10, 291, 350));
    EXPECT_EQ(4, rational_digit(9, 10, 291, 350));
    EXPECT_EQ(2, rational_digit(10, 10, 291, 350));
    EXPECT_EQ(8, rational_digit(11, 10, 291, 350));
    EXPECT_EQ(5, rational_digit(12, 10, 291, 350));
    EXPECT_EQ(7, rational_digit(13, 10, 291, 350));
}

TEST(fractions_test, rational_digits) {
    EXPECT_EQ((std::vector<int>{ 3, 4, 5, 6, 7, 8, 9, 0, 5, 6, 7, 8, 9, 0, 5, 6, 7, 8, 9, 0, 5, 6, 7, 8, 9, 0, 5, 6, 7, 8 }), rational_digits<int64_t>(2, 30, 10, 77160416, 624999375));
}
    
//cout << from_digits<int>(integer_digits(12345678, 10), 10) << endl;
