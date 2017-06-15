#pragma once

#include "altruct/algorithm/math/primes.h"
#include "altruct/algorithm/math/ranges.h"

#include <algorithm>
#include <vector>

namespace altruct {
namespace math {

class prime_holder {
private:
    typedef std::pair<int, int> fact_pair;

    int sz;                // upper bound (exclusive)
    int m;                 // number of primes up to sz
    std::vector<int> vp;   // primes
    std::vector<char> vq;  // prime flags
    std::vector<int> vpf;  // biggest prime factor
    std::vector<int> vpi;  // prime pi
    std::vector<int> vphi; // euler phi (totient)
    std::vector<int> vmu;  // moebius mu
    std::vector<int> vmer; // mertens

    void ensure_pq();
    std::vector<int>& ensure(std::vector<int> &v, void(*f)(int*, int, const int*, int));

public:
    prime_holder(int sz) : sz(sz) {}

    int size() { return sz; }
    int primes() { ensure_pq(); return m; }

    std::vector<int>& p() { ensure_pq(); return vp; }
    std::vector<char>& q() { ensure_pq(); return vq; }
    std::vector<int>& pf() { return ensure(vpf, altruct::math::factor); }
    std::vector<int>& pi() { return ensure(vpi, altruct::math::prime_pi); }
    std::vector<int>& phi() { return ensure(vphi, altruct::math::euler_phi); }
    std::vector<int>& mu() { return ensure(vmu, altruct::math::moebius_mu); }
    std::vector<int>& mertens() { if (vmer.empty()) { vmer = mu(); altruct::math::accumulate(vmer.begin(), vmer.end()); } return vmer; }

    int p(int i) { return p().at(i); }
    int q(int i) { return q().at(i); }
    int pf(int i) { return pf().at(i); }
    int pi(int i) { return pi().at(i); }
    int phi(int i) { return phi().at(i); }
    int mu(int i) { return mu().at(i); }
    int mertens(int i) { return mertens().at(i); }

    std::vector<fact_pair> factor_integer(int n);
    std::vector<fact_pair> factor_integer(std::vector<int> vn);
    template<typename I = int>
    std::vector<I> divisors(int n, I maxd = 0);
    template<typename I = int>
    std::vector<I> divisors(const std::vector<int> &vn, I maxd = 0);
    template<typename I = int>
    std::vector<I> divisors(const std::vector<fact_pair> &vf, I maxd = 0);
};

inline void prime_holder::ensure_pq() {
    if (vq.empty()) {
        m = int(sz / (log(sz) - 1.1)) + 5; // upper bound on pi(sz)
        if (sz < 40) m = sz / 2 + 2; // more accurate for small sz
        vp.resize(m);
        vq.resize(sz);
        m = altruct::math::primes(vp.data(), vq.data(), sz);
        vp.resize(m);
    }
}

inline std::vector<int>& prime_holder::ensure(std::vector<int> &v, void(*f)(int*, int, const int*, int)) {
    if (v.empty()) {
        v.resize(sz);
        f(v.data(), sz, p().data(), primes());
    }
    return v;
}

inline std::vector<prime_holder::fact_pair> prime_holder::factor_integer(std::vector<int> vn) {
    std::vector<prime_holder::fact_pair> vf;
    altruct::math::factor_integer(vf, vn, pf().data());
    std::sort(vf.begin(), vf.end());
    return vf;
}

inline std::vector<prime_holder::fact_pair> prime_holder::factor_integer(int n) {
    std::vector<prime_holder::fact_pair> vf;
    altruct::math::factor_integer(vf, n, pf().data());
    std::sort(vf.begin(), vf.end());
    return vf;
}

template<typename I>
std::vector<I> prime_holder::divisors(int n, I maxd) {
    return divisors(factor_integer(n), maxd);
}

template<typename I>
std::vector<I> prime_holder::divisors(const std::vector<int> &vn, I maxd) {
    return divisors(factor_integer(vn), maxd);
}

template<typename I>
std::vector<I> prime_holder::divisors(const std::vector<prime_holder::fact_pair> &vf, I maxd) {
    std::vector<I> vd;
    altruct::math::divisors(vd, vf, maxd);
    sort(vd.begin(), vd.end());
    return vd;
}

} // math
} // altruct
