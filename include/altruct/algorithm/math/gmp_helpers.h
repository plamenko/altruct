#pragma once

#include <cstdint>
#include <iostream>
#include <gmpxx.h>
#include "altruct/algorithm/math/base.h"
#include "altruct/structure/math/modulo.h"

typedef mpz_class mpz;
typedef mpq_class mpq;
typedef mpf_class mpf;

// casts
template<> unsigned long altruct::math::castT<unsigned long, mpz>::of(const mpz& x) {
    unsigned long res = (unsigned long)x.get_ui();
    if (mpz_sgn(x.get_mpz_t()) < 0) res = -res;
    return res;
};
template<> long altruct::math::castT<long, mpz>::of(const mpz& x) {
    return (long)castOf<unsigned long>(x);
};
template<> unsigned int altruct::math::castT<unsigned int, mpz>::of(const mpz& x) {
    return (unsigned int)castOf<unsigned long>(x);
};
template<> int altruct::math::castT<int, mpz>::of(const mpz& x) {
    return (int)castOf<unsigned int>(x);
};
template<> unsigned long long altruct::math::castT<unsigned long long, mpz>::of(const mpz& x) {
    // allocating on the stack is kind of free, assuming 100 words is enough; possible buffer overrun
    unsigned long long res[100]; res[0] = 0;
    mpz_export(res, 0, -1, sizeof res[0], 0, 0, x.get_mpz_t());
    if (mpz_sgn(x.get_mpz_t()) < 0) res[0] = -res[0];
    return res[0];
};
template<> long long altruct::math::castT<long long, mpz>::of(const mpz& x) {
    return (long long)castOf<unsigned long long>(x);
};
template<> mpz altruct::math::castT<mpz, unsigned long long>::of(const unsigned long long& x) {
    mpz res(0);
    mpz_import(res.get_mpz_t(), 1, -1, sizeof x, 0, 0, &x);
    return res;
};
template<> mpz altruct::math::castT<mpz, long long>::of(const long long& x) {
    return (x < 0) ? -castOf<mpz>((unsigned long long) - x) : castOf<mpz>((unsigned long long)x);
};

//inline mpz int64_to_mpz(int64_t x) { return (mpz(int32_t(x >> 32)) << 32) | mpz(uint32_t(x & 0xFFFFFFFF)); }
//inline int64_t mpz_to_int64(mpz x) { mpz xh = x >> 32, xl = x & 0xFFFFFFFF; return (int64_t(xh.get_si()) << 32) | uint64_t(xl.get_ui()); }
inline mpz int64_to_mpz(int64_t x) { return altruct::math::castOf<mpz>(x); }
inline int64_t mpz_to_int64(mpz x) { return altruct::math::castOf<int64_t>(x); }
inline mpz z_mulmod(mpz x, mpz y, mpz m) { return x * y % m; }
inline int64_t z_mulmod(int64_t x, int64_t y, int64_t m) { return ((x | y) >> 31 == 0) ? (x * y) % m : mpz_to_int64(z_mulmod(int64_to_mpz(x), int64_to_mpz(y), int64_to_mpz(m))); }
inline mpz z_powmod(mpz x, mpz y, mpz m) { mpz r; mpz_powm(r.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t(), m.get_mpz_t()); return r; }
inline int64_t z_powmod(int64_t x, int64_t y, int64_t m) { return mpz_to_int64(z_powmod(int64_to_mpz(x), int64_to_mpz(y), int64_to_mpz(m))); }
inline int is_prime(mpz x, int iter) { return mpz_probab_prime_p(x.get_mpz_t(), iter); }
inline int is_prime(int64_t x, int iter) { return is_prime(int64_to_mpz(x), iter); }
inline int64_t next_prime(int64_t x) { mpz p; mpz_nextprime(p.get_mpz_t(), int64_to_mpz(x).get_mpz_t()); return mpz_to_int64(p); }
inline mpz z_gcd(mpz x, mpz y) { mpz r; mpz_gcd(r.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t()); return r; }
inline mpz z_lcm(mpz x, mpz y) { return x * y / z_gcd(x, y); }
inline mpz z_abs(mpz x) { mpz r; mpz_abs(r.get_mpz_t(), x.get_mpz_t()); return r; }
inline mpz z_cb(mpz x) { return x*x*x; }
inline mpz z_sq(mpz x) { return x*x; }
inline mpz z_sqrt(mpz x) { mpz r; mpz_sqrt(r.get_mpz_t(), x.get_mpz_t()); return r; }
inline mpz z_sqrtc(mpz x) { mpz r = z_sqrt(x); if (r*r < x) r++; return r; }
inline mpz z_pow(mpz x, int n) { mpz r; mpz_pow_ui(r.get_mpz_t(), x.get_mpz_t(), n); return r; }
inline mpz z_inverse(mpz x, mpz m) { mpz r; mpz_invert(r.get_mpz_t(), x.get_mpz_t(), m.get_mpz_t()); return r; }
inline int z_testbit(mpz z, int i) { return mpz_tstbit(z.get_mpz_t(), i); }
inline int f_int(mpf x) { return mpf_get_si(x.get_mpf_t()); }
inline mpf f_sqrt(int n) { mpf r; mpf_sqrt_ui(r.get_mpf_t(), n); return r; }
inline mpf f_pow(int n) { mpf r; mpf_sqrt_ui(r.get_mpf_t(), n); return r; }
inline mpf f_abs(mpf x) { mpf r; mpf_abs(r.get_mpf_t(), x.get_mpf_t()); return r; }
inline mpf f_floor(mpf x) { mpf r; mpf_floor(r.get_mpf_t(), x.get_mpf_t()); return r; }
template<typename T> inline mpf& f_set(mpf &x, const T &val, int prec) { x.set_prec(prec); x = val; return x; }
inline mpf& f_set(mpf &x, const mpf &val, int prec) { x.set_prec(prec ? prec : val.get_prec()); x = val; return x; }
inline mpz z_div_floor(mpz a, mpz b) { mpz q; mpz_fdiv_q(q.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t()); return q; }
inline mpz z_div_ceil(mpz a, mpz b) { mpz q; mpz_cdiv_q(q.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t()); return q; }

inline mpz operator % (const mpz& lhs, const int64_t& rhs) { return lhs % int64_to_mpz(rhs); }
inline mpz operator % (const int64_t& lhs, const mpz& rhs) { return int64_to_mpz(lhs) % rhs; }
template<> inline int64_t altruct::math::modulo_mul(const int64_t& x, const int64_t& y, const int64_t& M) { return z_mulmod(x, y, M); }

template<> inline mpz altruct::math::sqrtT(mpz x, mpz) { return z_sqrt(x); }
template<> inline void altruct::math::modulo_normalize(mpz* v, const mpz& M) { altruct::math::modulo_normalize_int<mpz>(v, M); }
template<> inline mpz altruct::math::modulo_sub(const mpz& x, const mpz&y, const mpz& M) { return modulo_sub_int(x, y, M); }
template<> inline mpz altruct::math::modulo_inv(const mpz& x, const mpz& M) { return modulo_inv_int(x, M); }
template<> inline mpz altruct::math::modulo_div(const mpz& x, const mpz& y, const mpz& M) { return modulo_div_int(x, y, M); }
template<> inline mpz altruct::math::modulo_power(const mpz& x, const mpz& y, const mpz& M) { return z_powmod(x, y, M); }
template<> inline mpz altruct::math::modulo_power(const mpz& x, const int64_t& y, const mpz& M) { return z_powmod(x, int64_to_mpz(y), M); }

std::ostream& operator << (std::ostream& os, const mpz& rhs) { return os << rhs.get_str(); }
std::ostream& operator << (std::ostream& os, const mpq& rhs) { return os << rhs.get_str(); }
std::ostream& operator << (std::ostream& os, const mpf& rhs) { mp_exp_t exp = 0; return os << "0." << rhs.get_str(exp) << "e" << exp; }
