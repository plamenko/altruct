#include "altruct/algorithm/math/prime_counting.h"
#include "altruct/algorithm/math/primes.h"
#include "altruct/structure/math/modulo.h"

#include "gtest/gtest.h"

#include <functional>

using namespace std;
using namespace altruct::math;
using namespace altruct::container;

namespace {
typedef moduloX<int> modx;
}

TEST(prime_counting_test, prime_power_sum_sqrt_modx) {
    vector<char> vq(500);
    primes(nullptr, vq.data(), (int)vq.size());
    for (int z = 0; z <= 3; z++) {
        vector<modx> vps;
        modx c = { 0, 1009 };
        for (int n = 0; n < vq.size(); n++) {
            vps.push_back(c += powT(modx(n, 1009), z) * vq[n]);
        }
        for (int n = 1; n < vq.size(); n++) {
            auto mps = prime_power_sum_sqrt(z, n, modx(1, 1009));
            vector<modx> ve, va;
            for (int k = 1; k <= n; k++) {
                ve.push_back(vps[n / k]);
                va.push_back(mps[n / k]);
            }
            EXPECT_EQ(ve, va) << "unexpected prime_sum_sqrt result at n = " << n << " z = " << z;
        }
    }
}

TEST(prime_counting_test, prime_sum_modx) {
    vector<char> vq(1000);
    primes(nullptr, vq.data(), (int)vq.size());
    vector<modx> ve, va;
    modx c = { 0, 1009 };
    for (int n = 0; n < vq.size(); n++) {
        ve.push_back(c += {n * vq[n], 1009});
        va.push_back(prime_sum(n, modx(1, 1009)));
    }
    EXPECT_EQ(ve, va);
}

TEST(prime_counting_test, prime_sum) {
    vector<char> vq(1000);
    primes(nullptr, vq.data(), (int)vq.size());
    vector<int64_t> ve, va;
    int c = 0;
    for (int n = 0; n < vq.size(); n++) {
        ve.push_back(c += n * vq[n]);
        va.push_back(prime_sum(n, 1LL));
    }
    EXPECT_EQ(ve, va);
}

TEST(prime_counting_test, prime_pi_sqrt) {
    vector<char> vq(1000);
    primes(nullptr, vq.data(), (int)vq.size());
    vector<int64_t> vpi;
    int c = 0;
    for (int n = 0; n < vq.size(); n++) {
        vpi.push_back(c += vq[n]);
    }
    for (int n = 1; n < vq.size(); n++) {
        auto mpi = prime_pi_sqrt(n);
        vector<int64_t> ve, va;
        for (int k = 1; k <= n; k++) {
            ve.push_back(vpi[n / k]);
            va.push_back(mpi[n / k]);
        }
        EXPECT_EQ(ve, va) << "unexpected prime_pi_sqrt result at " << n;
    }
}

TEST(prime_counting_test, prime_pi) {
    vector<char> vq(1000);
    primes(nullptr, vq.data(), (int)vq.size());
    vector<int64_t> ve, va;
    int c = 0;
    for (int n = 0; n < vq.size(); n++) {
        ve.push_back(c += vq[n]);
        va.push_back(prime_pi(n));
    }
    EXPECT_EQ(ve, va);
}

TEST(prime_counting_test, prime_pi13) {
    vector<char> vq(1000);
    primes(nullptr, vq.data(), (int)vq.size());
    vector<int64_t> ve1, va1, ve3, va3;
    int c1 = 0, c3 = 0;
    for (int n = 0; n < vq.size(); n++) {
        ve1.push_back(c1 += vq[n] && (n % 4 == 1));
        va1.push_back(prime_pi1(n));
        ve3.push_back(c3 += vq[n] && (n % 4 == 3));
        va3.push_back(prime_pi3(n));
    }
    EXPECT_EQ(ve1, va1);
    EXPECT_EQ(ve3, va3);
}
