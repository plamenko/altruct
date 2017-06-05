#include "altruct/algorithm/math/sequences.h"
#include "altruct/structure/math/modulo.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

typedef modulo<int, 107> mod;

namespace {
template<typename F, typename I, typename R = typename std::result_of<F(I)>::type>
vector<R> seq(F f, I n) {
    vector<R> v;
    for (I i = 0; i < n; i++) {
        v.push_back(f(i));
    }
    return v;
}
}

TEST(sequences_test, seqeunces) {
    EXPECT_EQ((vector<mod> {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}), seq(sequences::delta<mod, int64_t>, 20));
    EXPECT_EQ((vector<mod> {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}), seq(sequences::dirichlet_id<mod, int64_t>, 20));
    EXPECT_EQ((vector<mod> {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}), seq(sequences::zero<mod, int64_t>, 20));
    EXPECT_EQ((vector<mod> {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}), seq(sequences::one<mod, int64_t>, 20));
    EXPECT_EQ((vector<mod> {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}), seq(sequences::identity<mod, int64_t>, 20));
    EXPECT_EQ((vector<mod> {0, 1, 4, 9, 16, 25, 36, 49, 64, 81, 100, 14, 37, 62, 89, 11, 42, 75, 3, 40}), seq(sequences::square<mod, int64_t>, 20));
    EXPECT_EQ((vector<mod> {0, 1, 8, 27, 64, 18, 2, 22, 84, 87, 37, 47, 16, 57, 69, 58, 30, 98, 54, 11}), seq(sequences::cube<mod, int64_t>, 20));
    EXPECT_EQ((vector<mod> {0, 1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66, 78, 91, 105, 13, 29, 46, 64, 83}), seq(sequences::triangular<mod, int64_t>, 20));
    EXPECT_EQ((vector<mod> {0, 1, 4, 10, 20, 35, 56, 84, 13, 58, 6, 72, 43, 27, 25, 38, 67, 6, 70, 46}), seq(sequences::tetrahedral<mod, int64_t>, 20));
    EXPECT_EQ((vector<mod> {0, 1, 5, 14, 30, 55, 91, 33, 97, 71, 64, 78, 8, 70, 52, 63, 105, 73, 76, 9}), seq(sequences::pyramidal<mod, int64_t>, 20));
    EXPECT_EQ((vector<mod> {0, 1, 6, 19, 44, 85, 39, 17, 23, 61, 28, 35, 86, 78, 15, 8, 61, 71, 42, 85}), seq(sequences::octahedral<mod, int64_t>, 20));
    EXPECT_EQ((vector<mod> {0, 1, 20, 84, 6, 27, 67, 46, 98, 36, 101, 106, 78, 44, 31, 66, 69, 67, 87, 49}), seq(sequences::dodecahedral<mod, int64_t>, 20));
    EXPECT_EQ((vector<mod> {0, 1, 12, 48, 17, 41, 28, 100, 58, 24, 13, 40, 13, 54, 71, 79, 93, 21, 92, 0}), seq(sequences::icosahedral<mod, int64_t>, 20));
    EXPECT_EQ((vector<mod> {0, 1, 3, 5, 8, 10, 14, 16, 20, 23, 27, 29, 35, 37, 41, 45, 50, 52, 58, 60}), seq(sequences::sum_sigma0<mod, int64_t>, 20));
    EXPECT_EQ((vector<mod> {0, 1, 4, 8, 15, 21, 33, 41, 56, 69, 87, 99, 20, 34, 58, 82, 6, 24, 63, 83}), seq(sequences::sum_sigma1<mod, int64_t>, 20));
    EXPECT_EQ((vector<mod> {0, 1, 6, 16, 37, 63, 6, 56, 34, 18, 41, 56, 52, 8, 44, 90, 3, 79, 106, 40}), seq(sequences::sum_sigma2<mod, int64_t>, 20));
}
