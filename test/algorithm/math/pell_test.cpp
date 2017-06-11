#include "altruct/algorithm/collections/collections.h"
#include "altruct/algorithm/math/pell.h"
#include "altruct/algorithm/math/factorization.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;
using namespace altruct::collections;

namespace {
typedef quadraticX<int64_t> quadx;

vector<quadx> to_quadx(int D, const vector<pair<int64_t, int64_t>>& v) {
    return transform(v, [&](pair<int64_t, int64_t> q) {
    	return quadx(q.first, q.second, D);
    });
}

vector<pair<int64_t, int64_t>> zip(const vector<int64_t>&v1, const vector<int64_t>&v2) {
    EXPECT_EQ(v1.size(), v2.size());
    vector<pair<int64_t, int64_t>> v(v1.size());
    for (int i = 0; i < (int)v.size(); i++) {
        v[i] = { v1[i], v2[i] };
    }
    return v;
}
}

TEST(pell_test, PQa_13) {
    if (true) {
        int64_t D = 13, P0 = 0, Q0 = 1;
        std::vector<int64_t> P{ P0 }, Q{ Q0 }, a, A, B, G;
        EXPECT_EQ(5, pell_PQa<int64_t>(D, P, Q, a, A, B, G));
        EXPECT_EQ((vector<int64_t>{0, 3, 1, 2, 1, 3, 3}), P);
        EXPECT_EQ((vector<int64_t>{1, 4, 3, 3, 4, 1, 4}), Q);
        EXPECT_EQ((vector<int64_t>{3, 1, 1, 1, 1, 6}), a);
        EXPECT_EQ((vector<int64_t>{3, 4, 7, 11, 18, 119}), A);
        EXPECT_EQ((vector<int64_t>{1, 1, 2, 3, 5, 33}), B);
        EXPECT_EQ((vector<int64_t>{3, 4, 7, 11, 18, 119}), G);
        EXPECT_EQ(5, pell_PQa<int64_t>(D, P, Q, a, A, B, G));
        EXPECT_EQ((vector<int64_t>{0, 3, 1, 2, 1, 3, 3, 1, 2, 1, 3, 3}), P);
        EXPECT_EQ((vector<int64_t>{1, 4, 3, 3, 4, 1, 4, 3, 3, 4, 1, 4}), Q);
        EXPECT_EQ((vector<int64_t>{3, 1, 1, 1, 1, 6, 1, 1, 1, 1, 6}), a);
        EXPECT_EQ((vector<int64_t>{3, 4, 7, 11, 18, 119, 137, 256, 393, 649, 4287}), A);
        EXPECT_EQ((vector<int64_t>{1, 1, 2, 3, 5, 33, 38, 71, 109, 180, 1189}), B);
        EXPECT_EQ((vector<int64_t>{3, 4, 7, 11, 18, 119, 137, 256, 393, 649, 4287}), G);
    }
    if (true) {
        int64_t D = 13, P0 = 1, Q0 = 2;
        std::vector<int64_t> P{ P0 }, Q{ Q0 }, a, A, B, G;
        EXPECT_EQ(1, pell_PQa<int64_t>(D, P, Q, a, A, B, G));
        EXPECT_EQ((vector<int64_t>{1, 3, 3}), P);
        EXPECT_EQ((vector<int64_t>{2, 2, 2}), Q);
        EXPECT_EQ((vector<int64_t>{2, 3}), a);
        EXPECT_EQ((vector<int64_t>{2, 7}), A);
        EXPECT_EQ((vector<int64_t>{1, 3}), B);
        EXPECT_EQ((vector<int64_t>{3, 11}), G);
        EXPECT_EQ(1, pell_PQa<int64_t>(D, P, Q, a, A, B, G));
        EXPECT_EQ(1, pell_PQa<int64_t>(D, P, Q, a, A, B, G));
        EXPECT_EQ(1, pell_PQa<int64_t>(D, P, Q, a, A, B, G));
        EXPECT_EQ(1, pell_PQa<int64_t>(D, P, Q, a, A, B, G));
        EXPECT_EQ((vector<int64_t>{1, 3, 3, 3, 3, 3, 3}), P);
        EXPECT_EQ((vector<int64_t>{2, 2, 2, 2, 2, 2, 2}), Q);
        EXPECT_EQ((vector<int64_t>{2, 3, 3, 3, 3, 3}), a);
        EXPECT_EQ((vector<int64_t>{2, 7, 23, 76, 251, 829}), A);
        EXPECT_EQ((vector<int64_t>{1, 3, 10, 33, 109, 360}), B);
        EXPECT_EQ((vector<int64_t>{3, 11, 36, 119, 393, 1298}), G);
    }
    if (true) {
        int64_t D = 13, P0 = 11, Q0 = 108;
        std::vector<int64_t> P{ P0 }, Q{ Q0 }, a, A, B, G;
        EXPECT_EQ(5, pell_PQa<int64_t>(D, P, Q, a, A, B, G));
        EXPECT_EQ((vector<int64_t>{11, -11, 4, 2, 1, 3, 3, 1, 2}), P);
        EXPECT_EQ((vector<int64_t>{108, -1, 3, 3, 4, 1, 4, 3, 3}), Q);
        EXPECT_EQ((vector<int64_t>{0, 7, 2, 1, 1, 6, 1, 1}), a);
        EXPECT_EQ((vector<int64_t>{0, 1, 2, 3, 5, 33, 38, 71}), A);
        EXPECT_EQ((vector<int64_t>{1, 7, 15, 22, 37, 244, 281, 525}), B);
        EXPECT_EQ((vector<int64_t>{-11, 31, 51, 82, 133, 880, 1013, 1893}), G);
        EXPECT_EQ(5, pell_PQa<int64_t>(D, P, Q, a, A, B, G));
        EXPECT_EQ((vector<int64_t>{11, -11, 4, 2, 1, 3, 3, 1, 2, 1, 3, 3, 1, 2}), P);
        EXPECT_EQ((vector<int64_t>{108, -1, 3, 3, 4, 1, 4, 3, 3, 4, 1, 4, 3, 3}), Q);
        EXPECT_EQ((vector<int64_t>{0, 7, 2, 1, 1, 6, 1, 1, 1, 1, 6, 1, 1}), a);
        EXPECT_EQ((vector<int64_t>{0, 1, 2, 3, 5, 33, 38, 71, 109, 180, 1189, 1369, 2558}), A);
        EXPECT_EQ((vector<int64_t>{1, 7, 15, 22, 37, 244, 281, 525, 806, 1331, 8792, 10123, 18915}), B);
        EXPECT_EQ((vector<int64_t>{-11, 31, 51, 82, 133, 880, 1013, 1893, 2906, 4799, 31700, 36499, 68199}), G);
    }
    if (true) {
        int64_t D = 13, P0 = 43, Q0 = 108;
        std::vector<int64_t> P{ P0 }, Q{ Q0 }, a, A, B, G;
        EXPECT_EQ(5, pell_PQa<int64_t>(D, P, Q, a, A, B, G));
        EXPECT_EQ((vector<int64_t>{43, -43, 9, 3, 3, 1, 2, 1, 3}), P);
        EXPECT_EQ((vector<int64_t>{108, -17, 4, 1, 4, 3, 3, 4, 1}), Q);
        EXPECT_EQ((vector<int64_t>{0, 2, 3, 6, 1, 1, 1, 1}), a);
        EXPECT_EQ((vector<int64_t>{0, 1, 3, 19, 22, 41, 63, 104}), A);
        EXPECT_EQ((vector<int64_t>{1, 2, 7, 44, 51, 95, 146, 241}), B);
        EXPECT_EQ((vector<int64_t>{-43, 22, 23, 160, 183, 343, 526, 869}), G);
    }
}

TEST(pell_test, pell1) {
    int64_t x0 = 0, y0 = 0;
    // x^2 - 13 y^2 = +1
    EXPECT_EQ(5, pell1<int64_t>(13, +1, x0, y0));
    EXPECT_EQ(649, x0);
    EXPECT_EQ(180, y0);
    // x^2 - 13 y^2 = -1
    EXPECT_EQ(5, pell1<int64_t>(13, -1, x0, y0));
    EXPECT_EQ(18, x0);
    EXPECT_EQ(5, y0);
    // x^2 - 3 y^2 = +1
    EXPECT_EQ(2, pell1<int64_t>(3, +1, x0, y0));
    EXPECT_EQ(2, x0);
    EXPECT_EQ(1, y0);
    // x^2 - 3 y^2 = -1
    EXPECT_EQ(0, pell1<int64_t>(3, -1, x0, y0));
    // x^2 - 23 y^2 = +1
    EXPECT_EQ(4, pell1<int64_t>(23, +1, x0, y0));
    EXPECT_EQ(24, x0);
    EXPECT_EQ(5, y0);
    // x^2 - 23 y^2 = -1
    EXPECT_EQ(0, pell1<int64_t>(23, -1, x0, y0));
}

TEST(pell_test, pellS) {
    auto calc = [](int D, int N) {
        vector<int64_t> xc0, yc0;
        pellS<int64_t>(D, N, xc0, yc0);
        return sorted(zip(xc0, yc0));
    };

    // x^2 - 13 y^2 = +1
    EXPECT_EQ((vector<pair<int64_t, int64_t>>{{ 649, 180 }}), calc(13, +1));
    // x^2 - 13 y^2 = -1
    EXPECT_EQ((vector<pair<int64_t, int64_t>>{{ 18, 5 }}), calc(13, -1));
    // x^2 - 3 y^2 = +1
    EXPECT_EQ((vector<pair<int64_t, int64_t>>{{ 2, 1 }}), calc(3, +1));
    // x^2 - 3 y^2 = -1
    EXPECT_EQ((vector<pair<int64_t, int64_t>>{}), calc(3, -1));
    // x^2 - 23 y^2 = +1
    EXPECT_EQ((vector<pair<int64_t, int64_t>>{{ 24, 5 }}), calc(23, +1));
    // x^2 - 23 y^2 = -1
    EXPECT_EQ((vector<pair<int64_t, int64_t>>{}), calc(23, -1));

    // x^2 - 157 y^2 = +12
    EXPECT_EQ((vector<pair<int64_t, int64_t>>{{ 13, 1 }, { 10663, 851 }, { 579160, 46222 }, { 483790960, 38610722 }, { 26277068347, 2097138361 }, { 21950079635497, 1751807067011 }}), calc(157, +12));
    // x^2 - 157 y^2 = -12
    EXPECT_EQ((vector<pair<int64_t, int64_t>>{{ 50, 4 }, { 2719, 217 }, { 2271269, 181267 }, { 123363799, 9845503 }, { 103049745749, 8224265053 }, { 5597138921710, 446700316396 }}), calc(157, -12));
}

TEST(pell_test, pell) {
    auto calc = [](int D, int N) {
        vector<int64_t> xc0, yc0;
        pell<int64_t>(D, N, factor_integer_slow(abs(N)), xc0, yc0);
        return sorted(zip(xc0, yc0));
    };

    // x^2 - 13 y^2 = +1
    EXPECT_EQ((vector<pair<int64_t, int64_t>>{{ 649, 180 }}), calc(13, +1));
    // x^2 - 13 y^2 = -1
    EXPECT_EQ((vector<pair<int64_t, int64_t>>{{ 18, 5 }}), calc(13, -1));
    // x^2 - 3 y^2 = +1
    EXPECT_EQ((vector<pair<int64_t, int64_t>>{{ 2, 1 }}), calc(3, +1));
    // x^2 - 3 y^2 = -1
    EXPECT_EQ((vector<pair<int64_t, int64_t>>{}), calc(3, -1));
    // x^2 - 23 y^2 = +1
    EXPECT_EQ((vector<pair<int64_t, int64_t>>{{ 24, 5 }}), calc(23, +1));
    // x^2 - 23 y^2 = -1
    EXPECT_EQ((vector<pair<int64_t, int64_t>>{}), calc(23, -1));

    // x^2 - 157 y^2 = +12
    EXPECT_EQ((vector<pair<int64_t, int64_t>>{{ 13, 1 }, { 10663, 851 }, { 579160, 46222 }, { 483790960, 38610722 }, { 26277068347, 2097138361 }, { 21950079635497, 1751807067011 }}), calc(157, +12));
    // x^2 - 157 y^2 = -12
    EXPECT_EQ((vector<pair<int64_t, int64_t>>{{ 50, 4 }, { 2719, 217 }, { 2271269, 181267 }, { 123363799, 9845503 }, { 103049745749, 8224265053 }, { 5597138921710, 446700316396 }}), calc(157, -12));
    
    // x^2 - 13 y^2 = 108
    EXPECT_EQ((vector<pair<int64_t, int64_t>>{{ -15, 3 }, { -11, 1 }, { 11, 1 }, { 15, 3 }, { 24, 6 }, { 41, 11 }, { 80, 22 }, { 141, 39 }, { 249, 69 }, { 440, 122 }, { 869, 241 }, { 1536, 426 }}), calc(13, 108));
}

TEST(pell_test, pell1_solutions) {
    // x^2 - 13 y^2 = +1
    auto vs1 = sorted(pell1<int64_t>(13, +1, 0, 0, 1));
    EXPECT_EQ(to_quadx(13, { { 649, 180 } }), vs1);
    // x^2 - 13 y^2 = -1
    auto vs2 = sorted(pell1<int64_t>(13, -1, 0, 0, 1));
    EXPECT_EQ(to_quadx(13, { { 18, 5 } }), vs2);
    
    // x^2 - 13 y^2 = +1
    auto vs3 = sorted(pell1<int64_t>(13, +1, 0, 0, 4));
    EXPECT_EQ(to_quadx(13, { { 649, 180 }, { 842401, 233640 }, { 1093435849, 303264540 }, { 1419278889601, 393637139280 } }), vs3);
    // x^2 - 13 y^2 = -1
    auto vs4 = sorted(pell1<int64_t>(13, -1, 0, 0, 4));
    EXPECT_EQ(to_quadx(13, { { 18, 5 }, { 23382, 6485 }, { 30349818, 8417525 }, { 39394040382, 10925940965 } }), vs4);
}

TEST(pell_test, pell_solutions) {
    // x^2 - 13 y^2 = +1
    auto vs3 = sorted(pell<int64_t>(13, +1, factor_integer_slow(1), 0, 0, 4));
    EXPECT_EQ(to_quadx(13, { { 649, 180 }, { 842401, 233640 }, { 1093435849, 303264540 }, { 1419278889601, 393637139280 } }), vs3);
    // x^2 - 13 y^2 = -1
    auto vs4 = sorted(pell<int64_t>(13, -1, factor_integer_slow(1), 0, 0, 4));
    EXPECT_EQ(to_quadx(13, { { 18, 5 }, { 23382, 6485 }, { 30349818, 8417525 }, { 39394040382, 10925940965 } }), vs4);
    
    // x^2 - 13 y^2 = 108
    auto vs12 = sorted(pell<int64_t>(13, 108, factor_integer_slow(108), 0, 0, 12));
    EXPECT_EQ(to_quadx(13, { { 11, 1 }, { 11, 1 }, { 15, 3 }, { 15, 3 }, { 24, 6 }, { 41, 11 }, { 80, 22 }, { 141, 39 }, { 249, 69 }, { 440, 122 }, { 869, 241 }, { 1536, 426 } }), vs12);

    auto vs10K = sorted(pell<int64_t>(13, 108, factor_integer_slow(108), 10000, 0, 0));
    EXPECT_EQ(to_quadx(13, { { 11, 1 }, { 11, 1 }, { 15, 3 }, { 15, 3 }, { 24, 6 }, { 41, 11 }, { 80, 22 }, { 141, 39 }, { 249, 69 }, { 440, 122 }, { 869, 241 }, { 1536, 426 }, { 2715, 753 }, { 4799, 1331 }, { 9479, 2629 } }), vs10K);

    auto vs50K = sorted(pell<int64_t>(13, 108, factor_integer_slow(108), 0, 50000, 0));
    EXPECT_EQ(to_quadx(13, { { 11, 1 }, { 11, 1 }, { 15, 3 }, { 15, 3 }, { 24, 6 }, { 41, 11 }, { 80, 22 }, { 141, 39 }, { 249, 69 }, { 440, 122 }, { 869, 241 }, { 1536, 426 }, { 2715, 753 }, { 4799, 1331 }, { 9479, 2629 }, { 16755, 4647 }, { 29616, 8214 }, { 52349, 14519 }, { 103400, 28678 } }), vs50K);
}
