#include "algorithm/math/convergents.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

TEST(convergents_test, sqrt_convergent) {
    EXPECT_EQ((pair<int64_t, int64_t>(2140758220993, 1513744654945)), (sqrt_convergent<int64_t>(2, 1000000000000)));
    EXPECT_EQ((pair<int64_t, int64_t>(2140758220993, 1513744654945)), (sqrt_convergent<int64_t>(2, 1513744654945 - 1)));
    EXPECT_EQ((pair<int64_t, int64_t>(7454517039243, 1099108574456)), (sqrt_convergent<int64_t>(46, 1000000000000)));
    EXPECT_EQ((pair<int64_t, int64_t>(105503093353351, 9512893564020)), (sqrt_convergent<int64_t>(123, 1000000000000)));
}

TEST(convergents_test, continued_fraction) {
    EXPECT_EQ((vector<int64_t>{1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}), (continued_fraction<int64_t>(2140758220993, 1513744654945)));
    EXPECT_EQ((vector<int64_t>{6, 1, 3, 1, 1, 2, 6, 2, 1, 1, 3, 1, 12, 1, 3, 1, 1, 2, 6, 2, 1, 1, 3, 1, 12, 1, 3, 1, 1, 2, 6, 3}), (continued_fraction<int64_t>(7454517039243, 1099108574456)));
    EXPECT_EQ((vector<int64_t>{11, 11, 22, 11, 22, 11, 22, 11, 22, 11, 22, 11}), (continued_fraction<int64_t>(105503093353351, 9512893564020)));
}

TEST(convergents_test, convergents) {
    // only convergents
    EXPECT_EQ((vector<pair<int64_t, int64_t>>{
        { 1, 1 }, { 3, 2 }, { 7, 5 }, { 17, 12 }, { 41, 29 }, { 99, 70 }, { 239, 169 }, { 577, 408 },
        { 1393, 985 }, { 3363, 2378 }, { 8119, 5741 }, { 19601, 13860 }, { 47321, 33461 }, {114243, 80782} }),
        (convergents<int64_t>({ 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 })));
    // all the best approximations
    EXPECT_EQ((vector<pair<int64_t, int64_t>>{
        { 1, 1 }, { 2, 1 }, { 3, 2 }, { 4, 3 }, { 7, 5 }, { 10, 7 }, { 17, 12 }, { 24, 17 }, { 41, 29 }, { 58, 41 }, { 99, 70 },
        { 140, 99 }, { 239, 169 }, { 338, 239 }, { 577, 408 }, { 816, 577 }, { 1393, 985 }, { 1970, 1393 }, { 3363, 2378 }, { 4756, 3363 },
        { 8119, 5741 }, { 11482, 8119 }, { 19601, 13860 }, { 27720, 19601 }, { 47321, 33461 }, { 66922, 47321 }, { 114243, 80782 } }),
        (convergents<int64_t>({ 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 }, 1000000000)));
    // only convergents
    EXPECT_EQ((vector<pair<int64_t, int64_t>>{
        { 6, 1 }, { 7, 1 }, { 27, 4 }, { 34, 5 }, { 61, 9 }, { 156, 23 }, { 997, 147 }, { 2150, 317 }, { 3147, 464 },
        { 5297, 781 }, { 19038, 2807 }, { 24335, 3588 }, { 311058, 45863 }, { 335393, 49451 }, { 1317237, 194216 },
        { 1652630, 243667 }, { 2969867, 437883 }, { 7592364, 1119433 }, { 48524051, 7154481 }, { 104640466, 15428395 },
        { 153164517, 22582876 }, { 257804983, 38011271 }, { 926579466, 136616689 }, { 1184384449, 174627960 }, { 15139192854, 2232152209 } }),
        (convergents<int64_t>({ 6, 1, 3, 1, 1, 2, 6, 2, 1, 1, 3, 1, 12, 1, 3, 1, 1, 2, 6, 2, 1, 1, 3, 1, 12 })));
}

TEST(convergents_test, line_closest_lattice_point) {
    int u = 10;
    for (int a = -u; a <= u; a++) {
        for (int b = -u; b <= u; b++) {
            for (int c = -u; c <= u; c++) {
                for (int l = -u; l <= u; l++) {
                    for (int r = l; r <= u; r++) {
                        int x0 = line_closest_lattice_point(a, b, c, l, r);
                        int y0 = (b == 0) ? 0 : div_round(a * x0 + c, -b);
                        int d0 = absT(a * x0 + b * y0 + c);
                        for (int x1 = l; x1 <= r; x1++) {
                            int y1 = (b == 0) ? 0 : div_round(a * x1 + c, -b);
                            int d1 = absT(a * x1 + b * y1 + c);
                            if (d1 < d0) EXPECT_FALSE(d1 < d0) << "ERROR: (" << a << " " << b << " " << c << ") (" << l << " " << r << "): " << d1 << " < " << d0;
                        }
                    }
                }
            }
        }
    }
}

TEST(convergents_test, minimize_floor_ladder) {
    int u = 5;
    for (int a = -u; a <= u; a++) {
        for (int b = -u; b <= u; b++) {
            for (int c = -u; c <= u; c++) {
                for (int d = -u; d <= u; d++) {
                    for (int e = -u; e <= u; e++) {
                        if (e == 0) continue;
                        for (int l = -u; l <= u; l++) {
                            for (int r = l; r <= u; r++) {
                                int x0 = minimize_floor_ladder(a, b, c, d, e, l, r);
                                int y0 = div_floor(c * x0 + d, e);
                                int s0 = a * x0 + b * y0;
                                for (int x1 = l; x1 <= r; x1++) {
                                    int y1 = div_floor(c * x1 + d, e);
                                    int s1 = a * x1 + b * y1;
                                    if (s1 < s0) EXPECT_FALSE(s1 < s0) << "ERROR: (" << a << " " << b << " " << c << " " << d << " " << e << ") (" << l << " " << r << "): " << s1 << " < " << s0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
