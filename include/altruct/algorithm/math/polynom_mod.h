#pragma once

#include "modulos.h"
#include "altruct/algorithm/math/fft.h"
#include "altruct/structure/math/root_wrapper.h"
#include "altruct/structure/math/complex.h"
#include "altruct/structure/math/modulo.h"
#include "altruct/structure/math/polynom.h"
#include <array>

namespace altruct {
namespace math {

namespace {
// special modulo P is a prime of form q * 2^k for a large k; 2^31 < P < 2^32

template<typename MODP, typename MOD>
void convert_to_mod_P_hilo(MODP* hi, MODP* lo, const MOD* p, int l) {
    for (int i = 0; i <= l; i++) {
        hi[i].v = uint32_t(p[i].v) >> 16;
        lo[i].v = uint32_t(p[i].v) & 0xFFFF;
    }
}

template<typename MODP, typename MOD>
std::array<std::vector<MODP>, 3> polynom_mul_mod_P_hilo(const MOD* p1, int l1, const MOD* p2, int l2, int n, uint32_t primitive_root) {
    auto root = powT<MODP>(primitive_root, (MODP::M() - 1) / n);
    auto iroot = powT(root, n - 1);
    std::vector<MODP> hi1(n), lo1(n); convert_to_mod_P_hilo(hi1.data(), lo1.data(), p1, l1);
    std::vector<MODP> hi2(n), lo2(n); convert_to_mod_P_hilo(hi2.data(), lo2.data(), p2, l2);
    auto ni = MODP(n).inv();
    fft(hi1.data(), n, root);
    fft(lo1.data(), n, root);
    fft(hi2.data(), n, root);
    fft(lo2.data(), n, root);
    for (int i = 0; i < n; i++) {
        MODP ss = (hi1[i] + lo1[i]) * (hi2[i] + lo2[i]) * ni;
        hi1[i] *= hi2[i] * ni;
        lo1[i] *= lo2[i] * ni;
        hi2[i] = ss - hi1[i] - lo1[i];
    }
    fft(hi1.data(), n, iroot); // hi1 * hi2
    fft(hi2.data(), n, iroot); // hi1 * lo2 + lo1 * hi2
    fft(lo1.data(), n, iroot); // lo1 * lo2
    return { std::move(lo1), std::move(hi2), std::move(hi1) };
}
} // namespace

/**
 * polynom<modulo<int>> specialization
 */
template<uint64_t ID>
struct polynom_mul<modulo<int, ID, modulo_storage::CONSTANT>> {
    typedef modulo<int, ID, modulo_storage::CONSTANT> mod;
    typedef complex<double> cplx;

    static int next_pow2(int l) { int n = 1; while (n < l) n *= 2; return n; }
    static int64_t rnd(const cplx& z, int n) { return llround(z.a / n) % mod::M(); }

    // splits coefficients into three 10-bit blocks each to avoid overflow
    // works for `mod::M < 2^30` and `la, lb <= 2^30`;
    static void _mul_fft_big(mod* pr, int lr, const mod* pa, int la, const mod* pb, int lb) {
        int n = next_pow2(la + lb + 1);
        auto root = complex_root_wrapper<double>(n);
        root = powT(root, root.size / n);
        auto iroot = powT(root, n - 1);
        std::vector<cplx> a2(n), a1(n), a0(n), b2(n), b1(n), b0(n), tmp(n);
        for (int i = 0; i <= la; i++) a2[i] = pa[i].v >> 20, a1[i] = (pa[i].v >> 10) & 0x3FF, a0[i] = pa[i].v & 0x3FF;
        for (int i = 0; i <= lb; i++) b2[i] = pb[i].v >> 20, b1[i] = (pb[i].v >> 10) & 0x3FF, b0[i] = pb[i].v & 0x3FF;
        fft_rec(tmp.data(), a2.data(), n, root); std::swap(a2, tmp);
        fft_rec(tmp.data(), a1.data(), n, root); std::swap(a1, tmp);
        fft_rec(tmp.data(), a0.data(), n, root); std::swap(a0, tmp);
        fft_rec(tmp.data(), b2.data(), n, root); std::swap(b2, tmp);
        fft_rec(tmp.data(), b1.data(), n, root); std::swap(b1, tmp);
        fft_rec(tmp.data(), b0.data(), n, root); std::swap(b0, tmp);
        for (int i = 0; i < n; i++) {
            auto w22 = a2[i] * b2[i];
            auto w11 = a1[i] * b1[i];
            auto w00 = a0[i] * b0[i];
            auto w21 = (a2[i] + a1[i]) * (b2[i] + b1[i]);
            auto w10 = (a1[i] + a0[i]) * (b1[i] + b0[i]);
            auto w210 = (a2[i] + a1[i] + a0[i]) * (b2[i] + b1[i] + b0[i]);
            a2[i] = w22, a1[i] = w11, a0[i] = w00, b2[i] = w21, b1[i] = w10, b0[i] = w210;
        }
        fft_rec(tmp.data(), a2.data(), n, iroot); std::swap(a2, tmp);
        fft_rec(tmp.data(), a1.data(), n, iroot); std::swap(a1, tmp);
        fft_rec(tmp.data(), a0.data(), n, iroot); std::swap(a0, tmp);
        fft_rec(tmp.data(), b2.data(), n, iroot); std::swap(b2, tmp);
        fft_rec(tmp.data(), b1.data(), n, iroot); std::swap(b1, tmp);
        fft_rec(tmp.data(), b0.data(), n, iroot); std::swap(b0, tmp);
        for (int i = 0; i <= lr; i++) {
            auto z11 = ((rnd(a2[i], n) << 20) + (rnd(b2[i] - a2[i] - a1[i], n) << 10) + rnd(a1[i], n));
            auto z01 = ((rnd(a2[i], n) << 20) + (rnd(b0[i] - a2[i] - b1[i], n) << 10) + rnd(b1[i], n));
            auto z00 = rnd(a0[i], n);
            pr[i] = ((z11 % mod::M() << 20) + ((z01 - z11 - z00) % mod::M() << 10) + z00) % mod::M();
        }
    }

    // splits coefficients into two 16-bit blocks each to avoid overflow
    // works for `mod::M < 2^30` and `l1, l2 <= 2^18`; e.g.: `M = 10^9+7, l = 250.000`
    static void _mul_fft(mod* pr, int lr, const mod* p1, int l1, const mod* p2, int l2) {
        int n = next_pow2(l1 + l2 + 1);
        auto root = complex_root_wrapper<double>(n);
        root = powT(root, root.size / n);
        auto iroot = powT(root, n - 1);
        std::vector<cplx> hi1(n), lo1(n), hi2(n), lo2(n), tmp(n);
        for (int i = 0; i <= l1; i++) hi1[i] = p1[i].v >> 16, lo1[i] = p1[i].v & 0xFFFF;
        for (int i = 0; i <= l2; i++) hi2[i] = p2[i].v >> 16, lo2[i] = p2[i].v & 0xFFFF;
        fft_rec(tmp.data(), hi1.data(), n, root); std::swap(hi1, tmp);
        fft_rec(tmp.data(), lo1.data(), n, root); std::swap(lo1, tmp);
        fft_rec(tmp.data(), hi2.data(), n, root); std::swap(hi2, tmp);
        fft_rec(tmp.data(), lo2.data(), n, root); std::swap(lo2, tmp);
        for (int i = 0; i < n; i++) {
            auto h = hi1[i] * hi2[i];
            auto l = lo1[i] * lo2[i];
            auto m = lo1[i] * hi2[i] + lo2[i] * hi1[i];
            lo1[i] = l, hi1[i] = m, hi2[i] = h;
        }
        fft_rec(tmp.data(), hi1.data(), n, iroot); std::swap(hi1, tmp);
        fft_rec(tmp.data(), lo1.data(), n, iroot); std::swap(lo1, tmp);
        fft_rec(tmp.data(), hi2.data(), n, iroot); std::swap(hi2, tmp);
        for (int i = 0; i <= lr; i++) {
            auto h = rnd(hi2[i], n), m = rnd(hi1[i], n), l = rnd(lo1[i], n);
            pr[i] = ((l << 0) + (m << 16) + (h << 32) % mod::M()) % mod::M();
        }
    }

    // splits coefficients into two 16-bit blocks each to avoid overflow; then
    // performs two separate convolutions modulo P1 and P2 and combines the result with CRT
    // works for: `mod::M < 2^32` and `l1 + l2 + 1 <= 2^28`
    // ~ 14 n log2 n + 18 n mad-operations (modulo mul + modulo add)
    static void _mul_fft_crt(mod* pr, int lr, const mod* p1, int l1, const mod* p2, int l2) {
        // n must divide 2^28, the largest power of two that divides both phi(P1) and phi(P2)
        int n = next_pow2(l1 + l2 + 1);
        const uint32_t P1 = UINT32_C(3221225473), root1 = 5u; // P1 = 3 * 2^30 + 1
        const uint32_t P2 = UINT32_C(3489660929), root2 = 3u; // P2 = 13 * 2^28 + 1
        using modP1 = modulo<uint32_t, P1, modulo_storage::CONSTANT>;
        using modP2 = modulo<uint32_t, P2, modulo_storage::CONSTANT>;
        auto hmlP1 = polynom_mul_mod_P_hilo<modP1>(p1, l1, p2, l2, n, root1);
        auto hmlP2 = polynom_mul_mod_P_hilo<modP2>(p1, l1, p2, l2, n, root2);
        modP1 P2i = -12; modP2 P1i = 13; const uint64_t PP = uint64_t(P1) * P2;
        auto crt = [&](modP1 v1, modP2 v2) {
            uint64_t r1 = uint64_t(P2) * (P2i * v1).v;
            uint64_t r2 = uint64_t(P1) * (P1i * v2).v;
            // r1 + r2 < 2 * PP, but can still overflow uint64_t
            return mod((r1 < PP - r2) ? r1 + r2 : (r1 + r2 - PP), p1->M());
        };
        mod w = powT(mod(2, p1->M()), 16); // 2^16
        for (int i = 0; i <= lr; i++) {
            mod lo = crt(hmlP1[0][i], hmlP2[0][i]);
            mod mi = crt(hmlP1[1][i], hmlP2[1][i]);
            mod hi = crt(hmlP1[2][i], hmlP2[2][i]);
            pr[i] = (hi * w + mi) * w + lo; // hi * 2^32 + mi * 2^16 + lo
        }
    }

    static void _mul_long(mod* pr, int lr, const mod* p1, int l1, const mod* p2, int l2) {
        for (int i = lr; i >= 0; i--) {
            int64_t r = 0;
            int jmax = std::min(i, l1);
            int jmin = std::max(0, i - l2);
            for (int j = jmax; j >= jmin; j--) {
                r += (int64_t(p1[j].v) * p2[i - j].v) % mod::M();
            }
            pr[i] = r % mod::M();
        }
    }

    static double cost_karatsuba(int l1, int l2) { return 0.25 * l1 * pow(l2, 0.5849625); }
    static double cost_fft(int l1, int l2) { int n = next_pow2(l1 + l2 + 1); return 0.5 * n * log2(n); }

    static void impl(mod* pr, int lr, const mod* p1, int l1, const mod* p2, int l2) {
        if (l2 < 16) {
            _mul_long(pr, lr, p1, l1, p2, l2);
        } else if (l2 < 300 || cost_karatsuba(l1, l2) < cost_fft(l1, l2)) {
            polynom<mod>::_mul_karatsuba(pr, lr, p1, l1, p2, l2);
        } else if (l1 <= 250000) {
            _mul_fft(pr, lr, p1, l1, p2, l2);
        } else {
            _mul_fft_big(pr, lr, p1, l1, p2, l2);
        }
    }
};

} // math
} // altruct
