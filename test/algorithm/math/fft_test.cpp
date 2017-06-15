#include "altruct/algorithm/math/fft.h"
#include "altruct/structure/math/modulo.h"

#include <algorithm>
#include <vector>

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

/**
 * {     prime, root, order}
 *
 * {1012924417,  198,  2^21}
 * {1004535809, 4172,  2^21}
 * { 985661441,  210,  2^22}
 * { 998244353,   31,  2^23}
 * { 897581057,   45,  2^23}
 * { 754974721,  362,  2^24}
 * { 469762049,   30,  2^26}
 *
 * {     12289,   41,  2^12}
 */

typedef modulo<int, 12289> mod;

TEST(fft_test, perf) {
    return; // skip perf tests

    typedef modulo<int, 754974721> mod;
    const int n = 1 << 16;
    mod root = powT(mod(362), (1 << 24) / n);
    vector<mod> a(n), t(n);

    for (int i = 0; i < n; i++) a[i] = i * i;
    auto T0 = clock();
    for (int i = 0; i < 100; i++)
        fft_rec(&t[0], &a[0], n, root);
    cout << "fft_rec: " << clock() - T0 << " ms" << endl;

    for (int i = 0; i < n; i++) a[i] = i * i;
    auto T1 = clock();
    for (int i = 0; i < 100; i++)
        fft(&a[0], n, root);
    cout << "fft: " << clock() - T1 << " ms" << endl;
}

TEST(fft_test, fft_rec_inverse) {
    const int n = 16;
    mod A[n], a[n] = { 671, 9230, 3302, 4764, 6135, 7750, 9881, 1189, 411, 8144 };
    fft_rec(A, a, n, powT(mod(41), (1 << 12) / n));
    EXPECT_EQ((vector<mod>{ 2321, 2621, 3262, 4649, 3137, 4957, 7242, 3878, 1612, 11833, 6116, 150, 9509, 964, 35, 9895 }), vector<mod>(A, A + n));
    fft_rec(a, A, n, powT(mod(1) / mod(41), (1 << 12) / n));
    for (auto& v : a) v /= n;
    EXPECT_EQ((vector<mod>{ 671, 9230, 3302, 4764, 6135, 7750, 9881, 1189, 411, 8144, 0, 0, 0, 0, 0, 0 }), vector<mod>(a, a + n));
}

TEST(fft_test, fft_inverse) {
    const int n = 16;
    mod a[n] = { 671, 9230, 3302, 4764, 6135, 7750, 9881, 1189, 411, 8144 };
    fft(a, n, powT(mod(41), (1 << 12) / n));
    EXPECT_EQ((vector<mod>{ 2321, 2621, 3262, 4649, 3137, 4957, 7242, 3878, 1612, 11833, 6116, 150, 9509, 964, 35, 9895 }), vector<mod>(a, a + n));
    fft(a, n, powT(mod(1) / mod(41), (1 << 12) / n));
    for (auto& v : a) v /= n;
    EXPECT_EQ((vector<mod>{ 671, 9230, 3302, 4764, 6135, 7750, 9881, 1189, 411, 8144, 0, 0, 0, 0, 0, 0 }), vector<mod>(a, a + n));
}

TEST(fft_test, fft_cyclic_convolution) {
    const int n = 16;
    mod u[n] = { 671, 9230, 3302, 4764, 6135, 7750, 9881, 1189, 411, 8144 };
    mod v[n] = { 8468, 3944, 4798, 6405, 8016, 8884, 1006, 54, 7066, 3531 };

    mod e[n] = { 0 };
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            e[k] += u[i] * v[modT(k - i, n)];
        }
    }

    mod a[n] = { 0 };
    fft_cyclic_convolution(a, u, v, n, mod(41), 1 << 12);
    EXPECT_EQ(vector<mod>(e, e + n), vector<mod>(a, a + n));
}

TEST(fft_test, fft_convolution) {
    const int n = 16;
    mod u[n] = { 671, 9230, 3302, 4764, 6135, 7750, 9881 };
    mod v[n] = { 8468, 3944, 4798, 6405, 8016, 8884, 1006, 54, 7066, 3531 };

    // necessary condition for cyclic convolution to equal normal convolution:
    // len(u) + len(v) <= n + 1;
    // 7 + 10 <= 16 + 1

    mod e[n] = { 0 };
    for (int k = 0; k < n; k++) {
        for (int i = 0; i <= k; i++) {
            e[k] += u[i] * v[k - i];
        }
    }

    mod a[n] = { 0 };
    fft_cyclic_convolution(a, u, v, n, mod(41), 1 << 12);
    EXPECT_EQ(vector<mod>(e, e + n), vector<mod>(a, a + n));
}

void test_ordinary_convolution(const vector<mod>& u, const vector<mod>& v) {
    vector<mod> eo(u.size() + v.size() - 1);
    for (int k = 0; k < (int)eo.size(); k++) {
        for (int i = 0; i < (int)u.size(); i++) {
            if (0 <= k - i && k - i < v.size()) {
                eo[k] += u[i] * v[k - i];
            }
        }
    }
    auto ao = convolution<mod>(u.begin(), u.end(), v.begin(), v.end(), mod(41), 1 << 12);
    EXPECT_EQ(eo, ao);
}
void test_cyclic_convolution(const vector<mod>& u, const vector<mod>& v) {
    vector<mod> ec(u.size() + v.size() - 1);
    for (int k = 0; k < (int)ec.size(); k++) {
        for (int i = 0; i < (int)u.size(); i++) {
            ec[k] += u[i] * v[modT(k - i, (int)v.size())];
        }
    }
    auto ac = cyclic_convolution<mod>(u.begin(), u.end(), v.begin(), v.end(), mod(41), 1 << 12);
    EXPECT_EQ(ec, ac);
}

void test_convolutions(const vector<mod>& u, const vector<mod>& v) {
    test_ordinary_convolution(u, v);
    test_cyclic_convolution(u, v);
}

TEST(fft_test, convolutions) {
    test_convolutions({ 671, 9230, 3302, 4764, 6135, 7750, 9881, 2930 }, { 8468, 3944, 4798, 6405, 8016, 8884, 1006, 54 });
    test_convolutions({ 671, 9230, 3302, 4764, 6135, 7750, 9881 }, { 8468, 3944, 4798, 6405, 8016, 8884, 1006, 54, 7066, 3531, 12, 3407, 551 });
    test_convolutions({ 8468, 3944, 4798, 6405, 8016, 8884, 1006, 54, 7066, 3531, 12, 3407, 551 }, { 671, 9230, 3302, 4764, 6135, 7750, 9881 });
    test_convolutions({ 671, 9230, 3302 }, { 8468, 3944, 4798, 6405, 8016, 8884, 1006, 54, 7066, 3531, 12, 3407, 551 });
    test_convolutions({ 8468, 3944, 4798, 6405, 8016, 8884, 1006, 54, 7066, 3531, 12, 3407, 551 }, { 671, 9230, 3302 });
    test_convolutions({ 671 }, { 8468, 3944, 4798, 6405, 8016, 8884, 1006, 54, 7066, 3531, 12, 3407, 551 });
    test_convolutions({ 8468, 3944, 4798, 6405, 8016, 8884, 1006, 54, 7066, 3531, 12, 3407, 551 }, { 671 });
}
