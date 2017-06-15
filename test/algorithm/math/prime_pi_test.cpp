#include "altruct/algorithm/math/prime_pi.h"
#include "altruct/structure/math/prime_holder.h"

#include <algorithm>

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

typedef long long ll;

TEST(prime_pi_test, prime_pi) {
    prime_holder prim(5000000);
    vector<ll> pi2 = { 0, 1, 2, 4, 6, 11, 18, 31, 54, 97, 172, 309, 564, 1028, 1900, 3512, 6542, 12251, 23000, 43390, 82025 };
    for (int k = 0; k < pi2.size(); k++) {
        prime_pi_PHI_tbl<ll, int>().clear();
        EXPECT_EQ(pi2[k], prime_pi_deprecated(powT(2LL, k), prim.pi().data(), prim.p().data())) << "2^" << k;
    }
    vector<ll> pi10 = { 0, 4, 25, 168, 1229, 9592, 78498, 664579 };
    for (int k = 0; k < pi10.size(); k++) {
        prime_pi_PHI_tbl<ll, int>().clear();
        EXPECT_EQ(pi10[k], prime_pi_deprecated(powT(10LL, k), prim.pi().data(), prim.p().data())) << "10^" << k;
    }
    for (int m = 0; m < 1000; m++) {
        EXPECT_EQ(prim.pi(m), prime_pi_deprecated(m, prim.pi().data(), prim.p().data())) << "pi(" << m << ")";
    }
}
