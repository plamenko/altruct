#include "altruct/concurrency/concurrency.h"

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
