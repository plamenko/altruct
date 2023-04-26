#include "altruct/structure/math/modulo.h"
#include "altruct/structure/math/polynom.h"
#include "altruct/algorithm/math/polynom_mod.h"
#include "altruct/algorithm/search/binary_search.h"
#include "altruct/algorithm/random/xorshift.h"
#include "altruct/chrono/chrono.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

namespace {
constexpr bool kTestLarge = false; // slow

enum class Algorithm { Long, Karatsuba, FFT_Double_Split2, FFT_Double_Split3, FFT_CRT };

template<typename MOD>
polynom<MOD> do_polynom_mul(Algorithm a, const polynom<MOD>& p1, const polynom<MOD>& p2) {
    const int l1 = p1.deg(), l2 = p2.deg(), lr = l1 + l2;
    polynom<MOD> pr; pr.resize(lr + 1, p1.ZERO_COEFF);
    switch (a) {
    case Algorithm::Long:
        polynom<MOD>::_mul_long(pr.c.data(), lr, p1.c.data(), l1, p2.c.data(), l2);
        break;
    case Algorithm::Karatsuba:
        polynom<MOD>::_mul_karatsuba(pr.c.data(), lr, p1.c.data(), l1, p2.c.data(), l2);
        break;
    case Algorithm::FFT_Double_Split2:
        polynom_mul<MOD>::_mul_fft(pr.c.data(), lr, p1.c.data(), l1, p2.c.data(), l2);
        break;
    case Algorithm::FFT_Double_Split3:
        polynom_mul<MOD>::_mul_fft_big(pr.c.data(), lr, p1.c.data(), l1, p2.c.data(), l2);
        break;
    case Algorithm::FFT_CRT:
        polynom_mul<MOD>::_mul_fft_crt(pr.c.data(), lr, p1.c.data(), l1, p2.c.data(), l2);
        break;
    default:
        break;
    }
    return pr;
}

// a ((b x)^(l+1) - 1) / (b x - 1) == a + a b x + a b^2 x^2 + ... + a b^l x^l
template<typename MOD>
polynom<MOD> make_poly_0(int l, int a, int b, MOD zero) {
    polynom<MOD> p;
    p.resize(l + 1, zero);
    p[0] = castOf(zero, a);
    for (int i = 1; i <= l; i++) {
        p[i] = p[i - 1] * castOf(zero, b);
    }
    return p;
}
template<typename MOD>
MOD eval_poly_0(int l, int a, int b, MOD x) {
    MOD am = castOf(x, a);
    MOD bx = castOf(x, b) * x;
    MOD e1 = identityOf(x);
    if (bx == e1) return am * castOf(x, l + 1);
    return am * (powT(bx, l + 1) - e1) / (bx - e1);
}
template<typename MOD>
bool test_polynom_mul_0(MOD zero, Algorithm a, int l1, int l2, int a1 = 7, int b1 = 3, int a2 = 2, int b2 = 9) {
    auto p1 = make_poly_0(l1, a1, b1, zero);
    auto p2 = make_poly_0(l2, a2, b2, zero);
    auto pr = do_polynom_mul(a, p1, p2);
    for (int x : {0, 1, -1, 2, -2, 10, -10}) {
        MOD xm = castOf(zero, x);
        MOD v1 = eval_poly_0(l1, a1, b1, xm);
        MOD v2 = eval_poly_0(l2, a2, b2, xm);
        MOD v = pr.eval(xm);
        if (v1 * v2 != v) return false;
    }
    return true;
}

// (a + b x)^l
template<typename MOD>
polynom<MOD> make_poly_1(int l, int a, int b, MOD zero) {
    // binomial(l,i) a^(l-i) b^i x^i
    // l!/(l-i)!/i! a^(l-i) b^i x^i
    // binomial is computed in two passes so we avoid
    // doing a modular inverse for each coefficient
    polynom<MOD> p;
    p.resize(l + 1, zero);
    MOD am = castOf(zero, a);
    MOD bm = castOf(zero, b);
    p[0] = powT(am, l);
    p[l] = powT(bm, l);
    if (am == zero) return p;
    MOD ba = bm / am;
    MOD f = identityOf(zero);
    for (int i = 1; i <= l; i++) {
        p[i] = p[i - 1] * castOf(zero, l - i + 1) * ba; // l!/(l-i)! a^(l-i) b^i
        f *= castOf(zero, i);
    }
    MOD fi = f.inv();
    for (int i = l; i >= 1; i--) {
        p[i] *= fi; // 1/i!
        fi *= castOf(zero, i);
    }
    return p;
}
template<typename MOD>
MOD eval_poly_1(int l, int a, int b, MOD x) {
    return powT(castOf(x, a) + castOf(x, b) * x, l);
}
template<typename MOD>
bool test_polynom_mul_1(MOD zero, Algorithm a, int l1, int l2, int a1 = 7, int b1 = 3, int a2 = 2, int b2 = 9) {
    // we could just exponentiate the polynomials,
    // but since this is the very logic under test,
    // the polynomials are constructed manually.
    // poly p1 = powT(poly{ a1, b1 }, l1); // (a1 + b1 x)^l1
    // poly p2 = powT(poly{ a2, b2 }, l2); // (a2 + b2 x)^l2
    auto p1 = make_poly_1(l1, a1, b1, zero); // (a1 + b1 x)^l1
    auto p2 = make_poly_1(l2, a2, b2, zero); // (a2 + b2 x)^l2
    auto pr = do_polynom_mul(a, p1, p2);
    for (int x : {0, 1, -1, 2, -2, 10, -10}) {
        MOD xm = castOf(zero, x);
        MOD v1 = eval_poly_1(l1, a1, b1, xm);
        MOD v2 = eval_poly_1(l2, a2, b2, xm);
        MOD v = pr.eval(xm);
        if (v1 * v2 != v) return false;
    }
    return true;
}

template<typename MOD>
bool test_polynom_mul(MOD zero, Algorithm a, int l1, int l2, int a1 = 7, int b1 = 3, int a2 = 2, int b2 = 9) {
    return test_polynom_mul_0(zero, a, l1, l2, a1, b1, a2, b2) &&
        test_polynom_mul_0(zero, a, l1, l2, -1, 1, -1, 1) &&
        test_polynom_mul_1(zero, a, l1, l2, a1, b1, a2, b2);
}

template<typename MOD>
int find_max_size(MOD zero, Algorithm a, int max_iter = 1000) {
    using clk = altruct::chrono::rdtsc_clock<>;
    auto T0 = clk::now();
    altruct::random::xorshift_64star rng(12345);
    int max_l1 = (1 << 28);
    int min_l1 = 1;
    while (min_l1 < max_l1 && test_polynom_mul(zero, a, min_l1, min_l1)) {
        cerr << min_l1 << " passed... " << since(T0) << " sec" << endl;
        min_l1 = min_l1 * 2 + 1;
    }
    min_l1 = min(min_l1, max_l1);
    cerr << min_l1 << " failed... " << since(T0) << " sec" << endl;
    for (int iter = 0; iter < max_iter; iter++) {
        int a1 = rng.next() % 10000;
        int b1 = rng.next() % 10000;
        int a2 = rng.next() % 10000;
        int b2 = rng.next() % 10000;
        int l1 = altruct::search::binary_search_pred(1, min_l1, [&](int l1) {
            return !test_polynom_mul(zero, a, l1, l1, a1, b1, a2, b2);
        });
        if (l1 >= min_l1) continue;
        min_l1 = l1;
        iter = 0;
        cerr << l1 << " " << a1 << " " << b1 << " " << a2 << " " << b2 << " " << since(T0) << " sec" << endl;
    }
    return min_l1;
}

} // namespace

// 1000000007 = 10^9 + 7; commonly used prime smaller than 2^30
// 2147483629 = 2^31 - 19; second largest prime that fits int32_t

TEST(polynom_mod_test, polynom_mul__mod_int__long) {
    EXPECT_TRUE((test_polynom_mul(modulo<int, 1000000007, modulo_storage::CONSTANT>(0), Algorithm::Long, 10, 5)));
    EXPECT_TRUE((test_polynom_mul(modulo<int, 1000000007, modulo_storage::CONSTANT>(0), Algorithm::Long, 100, 30)));
}

TEST(polynom_mod_test, polynom_mul__mod_int__karatsuba) {
    EXPECT_TRUE((test_polynom_mul(modulo<int, 1000000007, modulo_storage::CONSTANT>(0), Algorithm::Karatsuba, 10, 5)));
    EXPECT_TRUE((test_polynom_mul(modulo<int, 1000000007, modulo_storage::CONSTANT>(0), Algorithm::Karatsuba, 100, 30)));
    EXPECT_TRUE((test_polynom_mul(modulo<int, 1000000007, modulo_storage::CONSTANT>(0), Algorithm::Karatsuba, 1000, 700)));
}

TEST(polynom_mod_test, polynom_mul__mod_int__fft_double_split2) {
    EXPECT_TRUE((test_polynom_mul(modulo<int, 1000000007, modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split2, 10, 5)));
    EXPECT_TRUE((test_polynom_mul(modulo<int, 1000000007, modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split2, 100, 30)));
    EXPECT_TRUE((test_polynom_mul(modulo<int, 1000000007, modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split2, 1000, 700)));
    EXPECT_TRUE((test_polynom_mul(modulo<int, 2147483629, modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split2, 1000, 700)));
    EXPECT_TRUE((test_polynom_mul(modulo<int, 2147483629, modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split2, 1000, 1000)));
    EXPECT_TRUE((test_polynom_mul(modulo<int, 2147483629, modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split2, 10000, 10000)));
    if (!kTestLarge) return;
    EXPECT_TRUE((test_polynom_mul(modulo<int, 2147483629, modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split2, 65535, 65535)));
    //find_max_size(modulo<int, 1000000007, modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split2);
    //find_max_size(modulo<int, 1073741789, modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split2);
    //find_max_size(modulo<int, 2147483629, modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split2);
}

TEST(polynom_mod_test, polynom_mul__mod_int__fft_double_split3) {
    EXPECT_TRUE((test_polynom_mul(modulo<int, 1000000007, modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split3, 4, 4)));
    EXPECT_TRUE((test_polynom_mul(modulo<int, 1000000007, modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split3, 10, 5)));
    EXPECT_TRUE((test_polynom_mul(modulo<int, 1000000007, modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split3, 100, 30)));
    EXPECT_TRUE((test_polynom_mul(modulo<int, 1000000007, modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split3, 1000, 700)));
    EXPECT_TRUE((test_polynom_mul(modulo<int, 2147483629, modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split3, 1000, 1000)));
    EXPECT_TRUE((test_polynom_mul(modulo<int, 2147483629, modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split3, 10000, 10000)));
    if (!kTestLarge) return;
    EXPECT_TRUE((test_polynom_mul(modulo<int, 2147483629, modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split3, 250000, 250000)));
    //find_max_size(modulo<int, 2147483629, modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split3);
}

TEST(polynom_mod_test, polynom_mul__mod_int__fft_crt) {
    EXPECT_TRUE((test_polynom_mul(modulo<int, 1000000007, modulo_storage::CONSTANT>(0), Algorithm::FFT_CRT, 4, 4)));
    EXPECT_TRUE((test_polynom_mul(modulo<int, 1000000007, modulo_storage::CONSTANT>(0), Algorithm::FFT_CRT, 10, 5)));
    EXPECT_TRUE((test_polynom_mul(modulo<int, 1000000007, modulo_storage::CONSTANT>(0), Algorithm::FFT_CRT, 100, 30)));
    EXPECT_TRUE((test_polynom_mul(modulo<int, 1000000007, modulo_storage::CONSTANT>(0), Algorithm::FFT_CRT, 1000, 700)));
    EXPECT_TRUE((test_polynom_mul(modulo<int, 2147483629, modulo_storage::CONSTANT>(0), Algorithm::FFT_CRT, 1000, 700)));
    EXPECT_TRUE((test_polynom_mul(modulo<int, 2147483629, modulo_storage::CONSTANT>(0), Algorithm::FFT_CRT, 1000, 1000)));
    EXPECT_TRUE((test_polynom_mul(modulo<int, 2147483629, modulo_storage::CONSTANT>(0), Algorithm::FFT_CRT, 10000, 10000)));
    if (!kTestLarge) return;
    EXPECT_TRUE((test_polynom_mul(modulo<int, 2147483629, modulo_storage::CONSTANT>(0), Algorithm::FFT_CRT, 250000, 250000)));
    //find_max_size(modulo<int, 2147483629, modulo_storage::CONSTANT>(0), Algorithm::FFT_CRT);
}

// 4294967291 = 2^32 - 5; largest prime that fits uint32_t

TEST(polynom_mod_test, polynom_mul__mod_uint32__long) {
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::Long, 10, 5)));
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::Long, 100, 30)));
}

TEST(polynom_mod_test, polynom_mul__mod_uint32__karatsuba) {
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::Karatsuba, 10, 5)));
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::Karatsuba, 100, 30)));
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::Karatsuba, 1000, 700)));
}

TEST(polynom_mod_test, polynom_mul__mod_uint32__fft_double_split2) {
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split2, 10, 5)));
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split2, 100, 30)));
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split2, 1000, 700)));
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split2, 1000, 1000)));
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split2, 10000, 10000)));
    if (!kTestLarge) return;
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split2, 65535, 65535)));
    //find_max_size(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split2);
}

TEST(polynom_mod_test, polynom_mul__mod_uint32__fft_double_split3) {
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split3, 4, 4)));
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split3, 10, 5)));
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split3, 100, 30)));
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split3, 1000, 700)));
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split3, 1000, 1000)));
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split3, 10000, 10000)));
    if (!kTestLarge) return;
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split3, 250000, 250000)));
    //find_max_size(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_Double_Split3);
}

TEST(polynom_mod_test, polynom_mul__mod_uint32__fft_crt) {
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_CRT, 4, 4)));
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_CRT, 10, 5)));
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_CRT, 100, 30)));
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_CRT, 1000, 700)));
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_CRT, 1000, 700)));
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_CRT, 1000, 1000)));
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_CRT, 10000, 10000)));
    if (!kTestLarge) return;
    EXPECT_TRUE((test_polynom_mul(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_CRT, 250000, 250000)));
    //find_max_size(modulo<uint32_t, UINT32_C(4294967291), modulo_storage::CONSTANT>(0), Algorithm::FFT_CRT);
}

TEST(polynom_mod_test, polynom_mul__modx_uint32__long) {
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::Long, 10, 5)));
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::Long, 100, 30)));
}

TEST(polynom_mod_test, polynom_mul__modx_uint32__karatsuba) {
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::Karatsuba, 10, 5)));
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::Karatsuba, 100, 30)));
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::Karatsuba, 1000, 700)));
}

TEST(polynom_mod_test, polynom_mul__modx_uint32__fft_double_split2) {
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_Double_Split2, 10, 5)));
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_Double_Split2, 100, 30)));
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_Double_Split2, 1000, 700)));
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_Double_Split2, 1000, 1000)));
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_Double_Split2, 10000, 10000)));
    if (!kTestLarge) return;
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_Double_Split2, 250000, 250000)));
    //find_max_size(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_Double_Split2);
}

TEST(polynom_mod_test, polynom_mul__modx_uint32__fft_double_split3) {
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_Double_Split3, 4, 4)));
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_Double_Split3, 10, 5)));
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_Double_Split3, 100, 30)));
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_Double_Split3, 1000, 700)));
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_Double_Split3, 1000, 1000)));
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_Double_Split3, 10000, 10000)));
    if (!kTestLarge) return;
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_Double_Split3, 250000, 250000)));
    //find_max_size(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_Double_Split3);
}

TEST(polynom_mod_test, polynom_mul__modx_uint32__fft_crt) {
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_CRT, 4, 4)));
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_CRT, 10, 5)));
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_CRT, 100, 30)));
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_CRT, 1000, 700)));
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_CRT, 1000, 700)));
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_CRT, 1000, 1000)));
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_CRT, 10000, 10000)));
    if (!kTestLarge) return;
    EXPECT_TRUE((test_polynom_mul(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_CRT, 250000, 250000)));
    //find_max_size(moduloX<uint32_t>(0, UINT32_C(4294967291)), Algorithm::FFT_CRT);
}
