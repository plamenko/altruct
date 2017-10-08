#include "altruct/concurrency/concurrency.h"
#include "altruct/algorithm/math/primes.h"
#include "altruct/algorithm/math/reduce.h"

#include "gtest/gtest.h"

#include <string>
#include <set>
#include <vector>
#include <ppl.h>

using namespace std;
using namespace concurrency;
using namespace altruct::concurrency;

typedef long long ll;

TEST(concurrency_test, lock_macro) {
    ll r = 0;
    mutex r_mutex;
    parallel_for(0, 10000, [&](ll i) {
        ll q = i * i;
        LOCK(r_mutex) { r += q; }
    });
    EXPECT_EQ(333283335000LL, r);
}

TEST(concurrency_test, parallel_execute) {
    struct pi_worker_provider {
        struct pi_worker {
            vector<int> vp;
            vector<char> vq;
            pi_worker(int64_t U) : vp(altruct::math::isqrt(U) + 1) {
                vp.resize(altruct::math::primes(vp.data(), nullptr, (int)vp.size()));
            }
            int64_t execute_job(const std::pair<int64_t, int64_t>& job) {
                vq.resize(job.second - job.first); // reuses vq instead of allocating each time
                altruct::math::segmented_q(vq.data(), job.first, job.second, vp.data(), (int)vp.size());
                return altruct::math::reduce_sum(vq, int64_t(0));
            }
        };
        int64_t U; // max value job can take, primes are preprocessed up to sqrt(U)
        pi_worker_provider(int64_t U) : U(U) {}
        pi_worker create_worker() { return pi_worker(U); }
    };
    
    int N = 10007, B = 100;
    add_result_collector<int64_t> rc;
    range_job_provider<int64_t> jp(0, N + 1, B);
    pi_worker_provider wp(N + 1);
    parallel_execute(rc, jp, wp, 4);
    EXPECT_EQ(1230, rc.result); // pi(10007) = 1230
}