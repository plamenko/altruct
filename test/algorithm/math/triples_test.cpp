#include "algorithm/collections/collections.h"
#include "algorithm/math/primes.h"
#include "algorithm/math/triples.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;
using namespace altruct::collections;

TEST(triples_test, pythagorean_triples) {
    EXPECT_EQ((vector<triple<int>>{{ 3, 4, 5 }, { 5, 12, 13 }, { 7, 24, 25 }, { 8, 15, 17 }, { 9, 40, 41 }, { 12, 35, 37 }, { 20, 21, 29 }}), sorted(pythagorean_triples(50, true)));
    EXPECT_EQ((vector<triple<int>>{{ 3, 4, 5 }, { 5, 12, 13 }, { 6, 8, 10 }, { 7, 24, 25 }, { 8, 15, 17 }, { 9, 12, 15 }, { 9, 40, 41 }, { 10, 24, 26 }, { 12, 16, 20 }, { 12, 35, 37 }, { 14, 48, 50 }, { 15, 20, 25 }, { 15, 36, 39 }, { 16, 30, 34 }, { 18, 24, 30 }, { 20, 21, 29 }, { 21, 28, 35 }, { 24, 32, 40 }, { 27, 36, 45 }, { 30, 40, 50 }}), sorted(pythagorean_triples(50, false)));
}

TEST(triples_test, pythagorean_triples_fixed_leg) {
    EXPECT_EQ((vector<triple<int>>{{ 60, 11, 61 }, { 60, 25, 65 }, { 60, 32, 68 }, { 60, 45, 75 }, { 60, 63, 87 }, { 60, 80, 100 }, { 60, 91, 109 }, { 60, 144, 156 }, { 60, 175, 185 }, { 60, 221, 229 }, { 60, 297, 303 }, { 60, 448, 452 }, { 60, 899, 901 }}), sorted(pythagorean_triples_fixed_leg(60, factor_integer_slow(60))));
}

TEST(triples_test, eisenstein_triples) {
    EXPECT_EQ((vector<triple<int>>{{ 1, 1, 1 }, { 3, 7, 8 }, { 5, 7, 8 }, { 5, 19, 21 }, { 7, 13, 15 }, { 7, 37, 40 }, { 8, 13, 15 }, { 11, 31, 35 }, { 13, 43, 48 }, { 16, 19, 21 }, { 24, 31, 35 }, { 33, 37, 40 }, { 35, 43, 48 }}), sorted(eisenstein_triples60(50, true)));
    EXPECT_EQ((vector<triple<int>>{{ 1, 1, 1 }, { 2, 2, 2 }, { 3, 3, 3 }, { 3, 7, 8 }, { 4, 4, 4 }, { 5, 5, 5 }, { 5, 7, 8 }, { 6, 6, 6 }, { 6, 14, 16 }, { 7, 7, 7 }, { 7, 13, 15 }, { 8, 8, 8 }, { 8, 13, 15 }, { 9, 9, 9 }, { 10, 10, 10 }, { 10, 14, 16 }, { 11, 11, 11 }, { 12, 12, 12 }, { 13, 13, 13 }, { 14, 14, 14 }, { 15, 15, 15 }, { 16, 16, 16 }, { 17, 17, 17 }, { 18, 18, 18 }, { 19, 19, 19 }, { 20, 20, 20 }}), sorted(eisenstein_triples60(20, false)));
}
