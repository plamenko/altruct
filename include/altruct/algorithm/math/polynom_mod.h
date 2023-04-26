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
    std::vector<MODP> hi1(n), lo1(n), hi2(n), lo2(n);
    convert_to_mod_P_hilo(hi1.data(), lo1.data(), p1, l1);
    convert_to_mod_P_hilo(hi2.data(), lo2.data(), p2, l2);
    auto ni = MODP(n).inv();
    fft(hi1.data(), n, root);
    fft(lo1.data(), n, root);
    fft(hi2.data(), n, root);
    fft(lo2.data(), n, root);
    for (int i = 0; i < n; i++) {
        MODP hi = hi1[i] * hi2[i] * ni;
        MODP lo = lo1[i] * lo2[i] * ni;
        MODP mi = (lo1[i] * hi2[i] + lo2[i] * hi1[i]) * ni;
        lo1[i] = lo, hi1[i] = mi, hi2[i] = hi;
    }
    fft(hi2.data(), n, iroot); // hi1 * hi2
    fft(hi1.data(), n, iroot); // hi1 * lo2 + lo1 * hi2
    fft(lo1.data(), n, iroot); // lo1 * lo2
    return { std::move(lo1), std::move(hi1), std::move(hi2) };
}
} // namespace

/**
 * polynom<modulo<integral_type, ...>> specialization
 */
template<typename I, uint64_t ID, int STORAGE_TYPE>
struct polynom_mul<modulo<I, ID, STORAGE_TYPE>, std::enable_if_t<std::is_integral<I>::value>> {
    typedef modulo<I, ID, STORAGE_TYPE> mod;
    typedef complex<double> cplx;

    static int next_pow2(int l) { int n = 1; while (n < l) n *= 2; return n; }
    static uint64_t rnd(const cplx& z, int n, uint32_t M) { return uint64_t(llround(z.a / n)) % M; }
    static uint64_t shl_sub(uint64_t v) { return (v << 10) - v; }

    static void convert_to_cplx_210(cplx* c2, cplx* c1, cplx* c0, const mod* p, int l) {
        for (int i = 0; i <= l; i++) {
            c2[i] = uint32_t(p[i].v) >> 20;
            c1[i] = (uint32_t(p[i].v) >> 10) & 0x3FF;
            c0[i] = uint32_t(p[i].v) & 0x3FF;
        }
    }

    // splits coefficients into three 10-bit blocks each to avoid overflow
    // works for `mod::M < 2^30` and `la, lb <= 2^30`;
    static void _mul_fft_big(mod* pr, int lr, const mod* pa, int la, const mod* pb, int lb) {
        I M = pa->M();
        int n = next_pow2(la + lb + 1);
        auto root = complex_root_wrapper<double>(n);
        root = powT(root, root.size / n);
        auto iroot = powT(root, n - 1);
        std::vector<cplx> a2(n), a1(n), a0(n), b2(n), b1(n), b0(n);
        convert_to_cplx_210(a2.data(), a1.data(), a0.data(), pa, la);
        convert_to_cplx_210(b2.data(), b1.data(), b0.data(), pb, lb);
        fft(a2.data(), n, root);
        fft(a1.data(), n, root);
        fft(a0.data(), n, root);
        fft(b2.data(), n, root);
        fft(b1.data(), n, root);
        fft(b0.data(), n, root);
        for (int i = 0; i < n; i++) {
            cplx w22 = a2[i] * b2[i];
            cplx w11 = a1[i] * b1[i];
            cplx w00 = a0[i] * b0[i];
            cplx w21 = (a2[i] + a1[i]) * (b2[i] + b1[i]);
            cplx w10 = (a1[i] + a0[i]) * (b1[i] + b0[i]);
            cplx w210 = (a2[i] + a1[i] + a0[i]) * (b2[i] + b1[i] + b0[i]);
            a2[i] = w22, a1[i] = w11, a0[i] = w00, b2[i] = w21, b1[i] = w10, b0[i] = w210;
        }
        fft(a2.data(), n, iroot);
        fft(a1.data(), n, iroot);
        fft(a0.data(), n, iroot);
        fft(b2.data(), n, iroot);
        fft(b1.data(), n, iroot);
        fft(b0.data(), n, iroot);
        mod w = powT(mod(2, M), 20); // 2^20
        for (int i = 0; i <= lr; i++) {
            // r = 2^40 * (w22)
            //   + 2^30 * (w21 - w22 - w11)
            //   + 2^20 * (w210 - w21 - w10 + 2 * w11)
            //   + 2^10 * (w10 - w11 - w00)
            //   + 2^00 * (w00)
            uint64_t z22 = shl_sub(rnd(a2[i], n, M));                    // (w22 << 10) - w22
            uint64_t z11 = shl_sub(rnd(a1[i], n, M)) + rnd(b1[i], n, M); // (w11 << 10) - w11 + w10
            uint64_t z00 = shl_sub(rnd(a0[i], n, M));                    // (w00 << 10) - w00
            uint64_t z21 = shl_sub(rnd(b2[i], n, M)) + rnd(b0[i], n, M); // (w21 << 10) - w21 + w210
            uint64_t z10 = shl_sub(z11);                                 // (z11 << 10) - z11
            // r = (z22 << 30) + (z21 << 20) - (z10 << 10) - z00;
            pr[i] = mod((z22 << 10) + z21, M) * w - mod((z10 << 10) + z00, M);
        }
    }

    static void convert_to_cplx_hilo(cplx* hi, cplx* lo, const mod* p, int l) {
        for (int i = 0; i <= l; i++) {
            hi[i].a = uint32_t(p[i].v) >> 16;
            lo[i].a = uint32_t(p[i].v) & 0xFFFF;
        }
    }
    
    // splits coefficients into two 16-bit blocks each to avoid overflow
    // works for `mod::M < 2^32` and `l1+l2 < 2^17`; e.g.: `M = 2^32 - 5,  l1+l2 < 132.072`
    // works for `mod::M < 2^31` and `l1+l2 < 2^18`; e.g.: `M = 2^31 - 19, l1+l2 < 262.144`
    static void _mul_fft(mod* pr, int lr, const mod* p1, int l1, const mod* p2, int l2) {
        I M = p1->M();
        int n = next_pow2(l1 + l2 + 1);
        auto root = complex_root_wrapper<double>(n);
        root = powT(root, root.size / n);
        auto iroot = powT(root, n - 1);
        std::vector<cplx> hi1(n), lo1(n), hi2(n), lo2(n);
        convert_to_cplx_hilo(hi1.data(), lo1.data(), p1, l1);
        convert_to_cplx_hilo(hi2.data(), lo2.data(), p2, l2);
        fft(hi1.data(), n, root);
        fft(lo1.data(), n, root);
        fft(hi2.data(), n, root);
        fft(lo2.data(), n, root);
        for (int i = 0; i < n; i++) {
            cplx hi = hi1[i] * hi2[i];
            cplx lo = lo1[i] * lo2[i];
            cplx mi = lo1[i] * hi2[i] + lo2[i] * hi1[i];
            lo1[i] = lo, hi1[i] = mi, hi2[i] = hi;
        }
        fft(hi1.data(), n, iroot);
        fft(lo1.data(), n, iroot);
        fft(hi2.data(), n, iroot);
        for (int i = 0; i <= lr; i++) {
            uint64_t hi = rnd(hi2[i], n, M);
            uint64_t mi = rnd(hi1[i], n, M);
            uint64_t lo = rnd(lo1[i], n, M);
            //pr[i] = mod((((hi << 16) + mi) << 16) + lo, M); // can overflow by 1 bit
            pr[i] = mod(hi << 32, M) + mod((mi << 16) + lo, M);
        }
    }

    // splits coefficients into two 16-bit blocks each to avoid overflow; then
    // performs two separate convolutions modulo P1 and P2 and combines the result with CRT
    // works for: `mod::M < 2^32` and `l1+l2 < 2^28`
    // ~ 14 n log2 n + 18 n mad-operations (modulo mul + modulo add)
    static void _mul_fft_crt(mod* pr, int lr, const mod* p1, int l1, const mod* p2, int l2) {
        // n must divide 2^28, the largest power of two that divides both phi(P1) and phi(P2)
        I M = p1->M();
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
            return mod((r1 < PP - r2) ? r1 + r2 : (r1 + r2 - PP), M);
        };
        mod w = powT(mod(2, M), 16); // 2^16
        for (int i = 0; i <= lr; i++) {
            mod lo = crt(hmlP1[0][i], hmlP2[0][i]);
            mod mi = crt(hmlP1[1][i], hmlP2[1][i]);
            mod hi = crt(hmlP1[2][i], hmlP2[2][i]);
            pr[i] = (hi * w + mi) * w + lo; // hi * 2^32 + mi * 2^16 + lo
        }
    }

    static double cost_karatsuba(int l1, int l2) { return 0.25 * l1 * pow(l2, 0.5849625); }
    static double cost_fft(int l1, int l2) { int n = next_pow2(l1 + l2 + 1); return 0.75 * n * log2(n); }

    //static int LONG_THRESHOLD;
    static void impl(mod* pr, int lr, const mod* p1, int l1, const mod* p2, int l2) {
        if (l2 < 48) { // LONG_THRESHOLD
            polynom<mod>::_mul_long(pr, lr, p1, l1, p2, l2);
        } else if (l2 < 275 || l1 < 900 || cost_karatsuba(l1, l2) < cost_fft(l1, l2)) {
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
