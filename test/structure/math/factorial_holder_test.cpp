#include "altruct/structure/math/factorial_holder.h"
#include "altruct/structure/math/modulo.h"

#include "gtest/gtest.h"

#include <vector>

using namespace std;
using namespace altruct::math;

typedef moduloX<int> modx;

TEST(factorial_holder_test, factorial_holder_modx17) {
    factorial_holder<modx> fh(17, modx(1, 17));

    EXPECT_EQ(17, fh.size());

    const vector<modx> expected_fact = { 1, 1, 2, 6, 7, 1, 6, 8, 13, 15, 14, 1, 12, 3, 8, 1, 16 };
    const vector<modx> expected_ifact = { 1, 1, 9, 3, 5, 1, 3, 15, 4, 8, 11, 1, 10, 6, 15, 1, 16 };
    const vector<modx> expected_inv = { 0, 1, 9, 6, 13, 7, 3, 5, 15, 2, 12, 14, 10, 4, 11, 8, 16 };
    EXPECT_EQ(expected_fact, fh.fact());
    EXPECT_EQ(expected_ifact, fh.ifact());
    EXPECT_EQ(expected_inv, fh.inv());

    vector<modx> actual_fact; for (int i = 0; i < fh.size(); i++) actual_fact.push_back(fh.fact(i));
    vector<modx> actual_ifact; for (int i = 0; i < fh.size(); i++) actual_ifact.push_back(fh.ifact(i));
    vector<modx> actual_inv; for (int i = 0; i < fh.size(); i++) actual_inv.push_back(fh.inv(i));
    EXPECT_EQ(expected_fact, actual_fact);
    EXPECT_EQ(expected_ifact, actual_ifact);
    EXPECT_EQ(expected_inv, actual_inv);

    vector<vector<modx>> bin(17, vector<modx>(17, modx(0, 17)));
    for (int n = 0; n < 17; n++) {
        bin[n][0] = modx(1, 17);
        EXPECT_EQ(modx(0, 17), fh.bin(n, -1)) << n << " " << -1;
        EXPECT_EQ(modx(1, 17), fh.bin(n, 0)) << n << " " << 0;
        for (int k = 1; k <= n; k++) {
            bin[n][k] = bin[n - 1][k - 1] + bin[n - 1][k];
            EXPECT_EQ(bin[n][k], fh.bin(n, k)) << n << " " << k;
        }
        EXPECT_EQ(modx(0, 17), fh.bin(n, n+1)) << n << " " << n+1;
    }
}
