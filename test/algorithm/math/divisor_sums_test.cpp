#include "altruct/algorithm/collections/collections.h"
#include "altruct/algorithm/math/divisor_sums.h"
#include "altruct/algorithm/math/mertens.h"
#include "altruct/algorithm/math/primes.h"
#include "altruct/algorithm/math/ranges.h"
#include "altruct/structure/math/modulo.h"

#include "gtest/gtest.h"

#include <functional>

using namespace std;
using namespace altruct::math;
using namespace altruct::container;
using namespace altruct::collections;

namespace {
typedef moduloX<int> modx;
typedef modulo<int, 1000000007> field;

vector<int> primes_table(int n) {
    vector<int> p(n);
    int m = primes(p.data(), nullptr, n);
    p.resize(m);
    return p;
}

vector<int> prime_factor_table(int n) {
    auto p = primes_table(n);
    vector<int> pf(n);
    factor(pf.data(), n, p.data(), (int)p.size());
    return pf;
}

vector<moduloX<int>> to_modx(int M, const vector<int>& v) {
    return transform(v, [&](int a) {
        return moduloX<int>(a, M);
    });
}

template<typename T>
function<T(int)> to_func(const vector<T>& v) {
    return[=](int n){ return v[n]; };
}

const int n = 21;

auto v_e = to_modx(1009, { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });
auto v_1 = to_modx(1009, { 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 });
auto v_id = to_modx(1009, { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 });
auto v_mu = to_modx(1009, { 0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1, 0, -1, 0, -1, 0 });
auto v_d = to_modx(1009, { 0, 1, 2, 2, 3, 2, 4, 2, 4, 3, 4, 2, 6, 2, 4, 4, 5, 2, 6, 2, 6 });
auto v_s = to_modx(1009, { 0, 1, 3, 4, 7, 6, 12, 8, 15, 13, 18, 12, 28, 14, 24, 24, 31, 18, 39, 20, 42 });
auto v_phi = to_modx(1009, { 0, 1, 1, 2, 2, 4, 2, 6, 4, 6, 4, 10, 4, 12, 6, 8, 8, 16, 6, 18, 8 });

auto f_e = [](int n){ return modx(int(n == 1), 1009); };
auto f_1 = [](int n){ return modx(1, 1009); };
auto f_id = [](int n){ return modx(n, 1009); };
auto f_mu = to_func(v_mu);
auto f_d = to_func(v_d);
auto f_s = to_func(v_s);
auto f_phi = to_func(v_phi);
}

TEST(divisor_sums_test, dirichlet_convolution) {
    vector<modx> phi(n); dirichlet_convolution(phi, f_id, f_mu, n);
    EXPECT_EQ(v_phi, phi);
}

TEST(divisor_sums_test, dirichlet_division) {
    vector<modx> phi1(n); dirichlet_division(phi1, f_id, f_1, n);
    EXPECT_EQ(v_phi, phi1);

    vector<modx> phi2(n); dirichlet_division(phi2, f_s, f_d, n);
    EXPECT_EQ(v_phi, phi2);

    auto f = to_func(to_modx(1009, { 0, 6, 34, 66, 156, 160, 408, 294, 680, 648, 1020, 682, 2016, 936, 1904, 2100, 2928, 1564, 4266, 1938, 5160 }));
    auto g = to_func(to_modx(1009, { 0, 2, 6, 12, 20, 30, 42, 56, 72, 90, 110, 132, 156, 182, 210, 240, 272, 306, 342, 380, 420 }));
    vector<modx> h(n); dirichlet_division(h, f, g, n);
    EXPECT_EQ(to_modx(1009, { 0, 3, 8, 15, 24, 35, 48, 63, 80, 99, 120, 143, 168, 195, 224, 255, 288, 323, 360, 399, 440 }), h);
}

TEST(divisor_sums_test, dirichlet_inverse) {
    auto f = [](int n){ return modx(n * (n + 2), 1009); };
    vector<modx> f_inv(n); dirichlet_inverse(f_inv, f, n);
    EXPECT_EQ(to_modx(1009, { 0, 673, 896, 671, 635, 893, 452, 1002, 435, 670, 269, 881, 113, 651, 573, 459, 441, 861, 678, 292, 861 }), f_inv);
}

TEST(divisor_sums_test, calc_multiplicative) {
    int n = 51;
    auto pa = primes_table(n);
    auto id_expected = to_modx(1009, { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50 });
    auto id_actual = to_modx(1009, { 0, 1, 2, 3, 4, 5, 1, 7, 8, 9, 1, 11, 1, 13, 1, 1, 16, 17, 1, 19, 1, 1, 1, 23, 1, 25, 1, 27, 1, 29, 1, 31, 32, 1, 1, 1, 1, 37, 1, 1, 1, 41, 1, 43, 1, 1, 1, 47, 1, 49, 1 });
    calc_multiplicative(id_actual, n, pa.data(), (int)pa.size());
    EXPECT_EQ(id_expected, id_actual);
}

TEST(divisor_sums_test, dirichlet_convolution_multiplicative) {
    auto pa = primes_table(n);
    vector<modx> phi(n); dirichlet_convolution_multiplicative(phi, f_id, f_mu, n, pa.data(), (int)pa.size());
    EXPECT_EQ(v_phi, phi);
}

TEST(divisor_sums_test, dirichlet_division_multiplicative) {
    auto pa = primes_table(n);
    vector<modx> phi1(n); dirichlet_division_multiplicative(phi1, f_id, f_1, n, pa.data(), (int)pa.size());
    EXPECT_EQ(v_phi, phi1);

    vector<modx> phi2(n); dirichlet_division_multiplicative(phi2, f_s, f_d, n, pa.data(), (int)pa.size());
    EXPECT_EQ(v_phi, phi2);
}

TEST(divisor_sums_test, dirichlet_inverse_multiplicative) {
    auto pa = primes_table(n);
    vector<modx> f_inv(n); dirichlet_inverse_multiplicative(f_inv, f_phi, n, pa.data(), (int)pa.size());
    EXPECT_EQ(to_modx(1009, { 0, 1, -1, -2, -1, -4, 2, -6, -1, -2, 4, -10, 2, -12, 6, 8, -1, -16, 2, - 18, 4 }), f_inv);
}

TEST(divisor_sums_test, calc_completely_multiplicative) {
    int n = 51;
    auto pf = prime_factor_table(n);
    auto id_expected = to_modx(1009, { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50 });
    auto id_actual = to_modx(1009, { 0, 1, 2, 3, 1, 5, 1, 7, 1, 1, 1, 11, 1, 13, 1, 1, 1, 17, 1, 19, 1, 1, 1, 23, 1, 1, 1, 1, 1, 29, 1, 31, 1, 1, 1, 1, 1, 37, 1, 1, 1, 41, 1, 43, 1, 1, 1, 47, 1, 1, 1 });
    calc_completely_multiplicative(id_actual, n, pf.data());
    EXPECT_EQ(id_expected, id_actual);
}

TEST(divisor_sums_test, dirichlet_convolution_completely_multiplicative) {
    auto pf = prime_factor_table(n);
    vector<modx> id(n); dirichlet_convolution_completely_multiplicative(id, f_mu, f_s, n, pf.data());
    EXPECT_EQ(v_id, id);
}

TEST(divisor_sums_test, dirichlet_division_completely_multiplicative) {
    auto pf = prime_factor_table(n);
    vector<modx> id(n); dirichlet_division_completely_multiplicative(id, f_phi, f_mu, n, pf.data());
    EXPECT_EQ(v_id, id);
}

TEST(divisor_sums_test, dirichlet_inverse_completely_multiplicative) {
    auto pf = prime_factor_table(n);
    auto f = [&](int n){ return v_mu[n] * n; };
    vector<modx> f_inv(n); dirichlet_inverse_completely_multiplicative(f_inv, f, n, pf.data());
    EXPECT_EQ(to_modx(1009, { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 }), f_inv);
}

TEST(divisor_sums_test, moebius_transform) {
    vector<int> actual(n); moebius_transform(actual, [](int n){ return n * (n + 2); }, n);
    EXPECT_EQ((vector<int>{0, 3, 5, 12, 16, 32, 28, 60, 56, 84, 80, 140, 104, 192, 156, 208, 208, 320, 228, 396, 304}), actual);
}

TEST(divisor_sums_test, moebius_transform_multiplicative) {
    auto pa = primes_table(n);
    vector<modx> actual(n); moebius_transform_multiplicative(actual, f_id, n, pa.data(), (int)pa.size());
    EXPECT_EQ(v_phi, actual);
}

TEST(divisor_sums_test, sieve_m_multiplicative) {
    auto pa = primes_table(n);
    auto t1 = to_func(to_modx(1009, { 0, 1, 2, 4, 6, 10, 12, 18, 22, 28, 32, 42, 46, 58, 64, 72, 80, 96, 102, 120, 128 }));
    auto p1 = [&](int n){ return modx(n * n * n, 1009); };
    vector<modx> actual1(n); sieve_m_multiplicative(actual1, t1, p1, n, pa.data(), (int)pa.size());
    EXPECT_EQ(to_modx(1009, { 0, 1, 1003, 978, 972, 851, 17, 689, 677, 629, 467, 155, 305, 138, 479, 477, 453, 601, 937, 150, 876 }), actual1);
    auto t2 = [&](int n){ return modx(n * n * (n + 1) * (n + 1) / 4, 1009); };
    auto p2 = [&](int n){ return modx(n * n, 1009); };
    vector<modx> actual2(n); sieve_m_multiplicative(actual2, t2, p2, n, pa.data(), (int)pa.size());
    EXPECT_EQ(to_modx(1009, { 0, 1, 5, 23, 55, 155, 227, 521, 777, 254, 654, 855, 422, 432, 599, 381, 411, 999, 925, 360, 533 }), actual2);
}

TEST(divisor_sums_test, make_sqrt_map) {
    int n = 200, M = 101;
    auto f = [&](int n) { return modx(n, M); };
    auto tbl = make_sqrt_map<int, modx>(f, n);
    for (int i = 1; i <= n; i++) {
        EXPECT_EQ(modx(n / i, M), tbl[n / i]) << "i:" << i;
        EXPECT_EQ(M, tbl[n / i].M()) << "i:" << i;
    }
    EXPECT_EQ(0, tbl[0].v);
    EXPECT_EQ(M, tbl[0].M());
}

TEST(divisor_sums_test, sieve_m) {
    auto t = [](int n){ return n * (n + 1) / 2; };
    vector<int> actual1(n); sieve_m(actual1, t, n);
    EXPECT_EQ((vector<int>{0, 1, 2, 4, 6, 10, 12, 18, 22, 28, 32, 42, 46, 58, 64, 72, 80, 96, 102, 120, 128}), actual1);

    typedef moduloX<int> modx;
    auto t2 = [](int n){ return modx(n * (n + 1) / 2, 1009); };
    auto p2 = [](int n){ return modx(n + 2, 1009); };
    vector<modx> actual2(n); sieve_m(actual2, t2, p2, n);
    EXPECT_EQ(to_modx(1009, { 0, 673, 449, 1, 973, 77, 264, 938, 540, 840, 205, 992, 170, 509, 61, 809, 482, 934, 112, 116, 490 }), actual2);
}

TEST(divisor_sums_test, sieve_sqfree_count) {
    int n = 31;
    auto pa = primes_table(isqrt(n) + 1);
    auto expected1 = vector<int>{0, 1, 2, 3, 3, 4, 5, 6, 6, 6, 7, 8, 8, 9, 10, 11, 11, 12, 12, 13, 13, 14, 15, 16, 16, 16, 17, 17, 17, 18, 19};
    vector<int> actual1(n); sieve_sqfree_count(actual1, n, pa.data(), (int)pa.size(), int(1));
    EXPECT_EQ(expected1, actual1);

    typedef moduloX<int> modx;
    auto expected2 = to_modx(1009, expected1);
    vector<modx> actual2(n); sieve_sqfree_count(actual2, n, pa.data(), (int)pa.size(), modx(1, 1009));
    EXPECT_EQ(expected2, actual2);
}

TEST(divisor_sums_test, sqfree_count) {
    int n = 30;
    auto v_S = to_modx(1009, { 0, 1, 2, 3, 3, 4, 5, 6, 6, 6, 7, 8, 8, 9, 10, 11, 11, 12, 12, 13, 13, 14, 15, 16, 16, 16, 17, 17, 17, 18, 19 });
    // preprocess `U = n^(3/5)` values
    int U = (int)pow(n, 0.6);
    sqrt_map<int, modx> ms(U, n);
    for (int k = 0; k < U; k++) {
        ms[k] = v_S[k];
    }
    // calc sqfree_count
    vector<modx> va;
    for (int k = 0; k <= n; k++) {
        ms.reset_max(k);
        va.push_back(sqfree_count(k, ms, modx(1, 1009)));
    }
    EXPECT_EQ(v_S, va);
}

TEST(divisor_sums_test, sum_phi_D_L) {
    auto id = field(1);
    auto vn = range<int64_t>(21);

    EXPECT_EQ((vector<field>{0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}), sum_phi_D_L(0, 0, vn, 0, id));
    EXPECT_EQ((vector<field>{0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}), sum_phi_D_L(0, 1, vn, 0, id));
    EXPECT_EQ((vector<field>{0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}), sum_phi_D_L(0, 2, vn, 0, id));

    EXPECT_EQ((vector<field>{0, 1, 2, 4, 6, 10, 12, 18, 22, 28, 32, 42, 46, 58, 64, 72, 80, 96, 102, 120, 128}), sum_phi_D_L(1, 0, vn, 0, id));
    EXPECT_EQ((vector<field>{0, 1, 3, 9, 17, 37, 49, 91, 123, 177, 217, 327, 375, 531, 615, 735, 863, 1135, 1243, 1585, 1745}), sum_phi_D_L(1, 1, vn, 0, id));
    EXPECT_EQ((vector<field>{0, 1, 5, 23, 55, 155, 227, 521, 777, 1263, 1663, 2873, 3449, 5477, 6653, 8453, 10501, 15125, 17069, 23567, 26767}), sum_phi_D_L(1, 2, vn, 0, id));

    EXPECT_EQ((vector<field>{0, 1, 3, 8, 15, 29, 42, 69, 95, 134, 172, 237, 287, 377, 452, 552, 652, 804, 915, 1104, 1252}), sum_phi_D_L(2, 0, vn, 0, id));
    EXPECT_EQ((vector<field>{0, 1, 5, 20, 48, 118, 196, 385, 593, 944, 1324, 2039, 2639, 3809, 4859, 6359, 7959, 10543, 12541, 16132, 19092}), sum_phi_D_L(2, 1, vn, 0, id));
    EXPECT_EQ((vector<field>{0, 1, 9, 54, 166, 516, 984, 2307, 3971, 7130, 10930, 18795, 25995, 41205, 55905, 78405, 104005, 147933, 183897, 252126, 311326}), sum_phi_D_L(2, 2, vn, 0, id));

    EXPECT_EQ(field(356214470), sum_phi_D_L(1, 0, 10000000, 0, id));
}

TEST(divisor_sums_test, sum_phi_D_L_modx) {
    typedef moduloX<int> modx;

    auto id = modx(1, 1009);
    auto vn = range<int64_t>(21);

    EXPECT_EQ(to_modx(1009, { 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }), sum_phi_D_L(0, 0, vn, 0, id));
    EXPECT_EQ(to_modx(1009, { 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }), sum_phi_D_L(0, 1, vn, 0, id));
    EXPECT_EQ(to_modx(1009, { 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }), sum_phi_D_L(0, 2, vn, 0, id));

    EXPECT_EQ(to_modx(1009, { 0, 1, 2, 4, 6, 10, 12, 18, 22, 28, 32, 42, 46, 58, 64, 72, 80, 96, 102, 120, 128 }), sum_phi_D_L(1, 0, vn, 0, id));
    EXPECT_EQ(to_modx(1009, { 0, 1, 3, 9, 17, 37, 49, 91, 123, 177, 217, 327, 375, 531, 615, 735, 863, 126, 234, 576, 736 }), sum_phi_D_L(1, 1, vn, 0, id));
    EXPECT_EQ(to_modx(1009, { 0, 1, 5, 23, 55, 155, 227, 521, 777, 254, 654, 855, 422, 432, 599, 381, 411, 999, 925, 360, 533 }), sum_phi_D_L(1, 2, vn, 0, id));

    EXPECT_EQ(to_modx(1009, { 0, 1, 3, 8, 15, 29, 42, 69, 95, 134, 172, 237, 287, 377, 452, 552, 652, 804, 915, 95, 243 }), sum_phi_D_L(2, 0, vn, 0, id));
    EXPECT_EQ(to_modx(1009, { 0, 1, 5, 20, 48, 118, 196, 385, 593, 944, 315, 21, 621, 782, 823, 305, 896, 453, 433, 997, 930 }), sum_phi_D_L(2, 1, vn, 0, id));
    EXPECT_EQ(to_modx(1009, { 0, 1, 9, 54, 166, 516, 984, 289, 944, 67, 840, 633, 770, 845, 410, 712, 78, 619, 259, 885, 554 }), sum_phi_D_L(2, 2, vn, 0, id));

    EXPECT_EQ(modx(984, 1009), sum_phi_D_L(1, 0, 10000000, 0, id));
}

TEST(divisor_sums_test, divisor_sigma) {
    int n = 30;
    auto pa = primes_table(n);
    vector<int> vds0(n);
    divisor_sigma<int>(vds0, 0, n, pa.data(), (int)pa.size());
    EXPECT_EQ((vector<int> {0, 1, 2, 2, 3, 2, 4, 2, 4, 3, 4, 2, 6, 2, 4, 4, 5, 2, 6, 2, 6, 4, 4, 2, 8, 3, 4, 4, 6, 2}), vds0);

    vector<int64_t> vds1(n);
    divisor_sigma<int64_t>(vds1, 1, n, pa.data(), (int)pa.size());
    EXPECT_EQ((vector<int64_t> {0, 1, 3, 4, 7, 6, 12, 8, 15, 13, 18, 12, 28, 14, 24, 24, 31, 18, 39, 20, 42, 32, 36, 24, 60, 31, 42, 40, 56, 30}), vds1);

    vector<int64_t> vds2(n);
    divisor_sigma<int64_t>(vds2, 2, n, pa.data(), (int)pa.size());
    EXPECT_EQ((vector<int64_t> {0, 1, 5, 10, 21, 26, 50, 50, 85, 91, 130, 122, 210, 170, 250, 260, 341, 290, 455, 362, 546, 500, 610, 530, 850, 651, 850, 820, 1050, 842 }), vds2);

    vector<modx> vds10(n);
    divisor_sigma(vds10, 10, n, pa.data(), (int)pa.size(), modx(1, 107));
    EXPECT_EQ(to_modx(107, { 0, 1, 62, 93, 38, 57, 95, 65, 72, 104, 3, 43, 3, 10, 71, 58, 6, 20, 28, 38, 26, 53, 98, 36, 62, 90, 85, 46, 9, 5 }), vds10);
}

TEST(divisor_sums_test, sum_multiplicative) {
    int M = 101;
    int n = 1000;
    int u = int(sqrt(n * log(n)));
    auto pa = primes_table(u);
    auto pa_all = primes_table(n + 1);
    auto zero = modx(0, M);
    // moebius
    auto g = [&](modx f_pe1, int p, int e) {
        return castOf(zero, (e > 1) ? 0 : -1);
    };
    // -primes_pi
    auto sg1 = [&](int64_t n) {
        int pi = int(std::upper_bound(pa_all.begin(), pa_all.end(), n) - pa_all.begin());
        return castOf(zero, -pi);
    };
    // mertens
    auto sg_tbl = sum_multiplicative<modx>(sg1, g, n, pa.data(), (int)pa.size(), identityOf(zero));
    vector<modx> v_mertens(n + 1); sieve_mertens(v_mertens, n + 1, pa_all.data(), (int)pa_all.size(), identityOf(zero));
    for (int i = 1; i <= n; i++) {
        int k = n / i;
        EXPECT_EQ(v_mertens[k], sg_tbl[k]) << "i:" << i;
        EXPECT_EQ(v_mertens[k], sum_multiplicative_34(sg1, g, k, pa.data(), (int)pa.size(), identityOf(zero))) << "i:" << i;
    }
}
