#include "altruct/algorithm/math/primes.h"

#include <cmath>
#include <climits>
#include <stdint.h>
#include <algorithm>

namespace altruct {
namespace math {

int primes(int *p, char *q, int n) {
    // if q not specified, use p for both
    // last n chars of p are used for q
    if (!q) q = (char *) (p + n) - n;
    // initialize q
    if (n > 0) q[0] = 0;
    if (n > 1) q[1] = 0;
    for (int i = 2; i < n; i++) {
        q[i] = 1;
    }
    // perform sieving
    int m = 0;
    for (int i = 2; i < n; i++) {
        if (!q[i]) continue;
        if (p) p[m] = i;
        m++;
        if (i > n / i) continue;
        for (int j = i * i; j < n; j += i)
            q[j] = 0;
    }
    if (p) p[m] = 0;
    return m;
}

void prime_pi(int *pi, int n, const int *p, int m) {
    for (int i = 0, l = 0; i < n; i++) {
        if (i == p[l]) l++;
        pi[i] = l;
    }
}

void euler_phi(int *phi, int n, const int *p, int m) {
    for (int i = 0; i < n; i++)
        phi[i] = i;
    for (int i = 0; i < m; i++)
        for (int j = p[i]; j < n; j += p[i])
            phi[j] = phi[j] / p[i] * (p[i] - 1);
}

void moebius_mu(int *mu, int n, const int* p, int m) {
    mu[0] = 0;
    for (int i = 1; i < n; i++)
        mu[i] = 1;
    for (int i = 2; i * i < n; i++) {
        if (mu[i] != 1) continue;
        int i2 = i * i;
        for (int j = 0; j < n; j += i2)
            mu[j] = 0;
        for (int j = 0; j < n; j += i)
            mu[j] *= -i;
    }
    for (int i = 2; i < n; i++) {
             if (mu[i] == +i) mu[i] = +1;
        else if (mu[i] == -i) mu[i] = -1;
        else if (mu[i] < 0) mu[i] = +1; // correction for a big prime factor
        else if (mu[i] > 0) mu[i] = -1; // correction for a big prime factor
    }
}

void prime_nu(int *nu, int n, const int* p, int m) {
    for (int i = 0; i < n; i++)
        nu[i] = 0;
    for (int i = 0; i < m; i++)
        for (int j = p[i]; j < n; j += p[i])
            nu[j]++;
}

void segmented_q(char* q, int64_t b, int64_t e, const int *p, int m) {
    char* _q = q - b;
    if (b == 0) _q[b++] = 0;
    if (b == 1) _q[b++] = 0;
    for (int64_t a = b; a < e; a++)
        _q[a] = 1;
    for (int i = 0; i < m; i++) {
        b = std::max(b, isq(p[i])); if (b >= e) break;
        for (int64_t a = multiple<int64_t>(p[i], b); a < e; a += p[i]) {
            _q[a] = 0;
        }
    }
}

void segmented_phi(int64_t *phi, int64_t *tmp, int64_t b, int64_t e, const int *p, int m) {
    int64_t *_phi = phi - b, *_tmp = tmp - b;
    if (b == 0) _phi[b++] = 0;
    for (int64_t q = b; q < e; q++)
        _phi[q] = 1, _tmp[q] = q;
    for (int i = 0; i < m; i++) {
        for (int64_t q = multiple<int64_t>(p[i], b); q < e; q += p[i]) {
            _phi[q] *= p[i] - 1, _tmp[q] /= p[i];
            while (_tmp[q] % p[i] == 0)
                _phi[q] *= p[i], _tmp[q] /= p[i];
        }
    }
    // correction for a large prime factor (p > sqrt(e))
    for (int64_t q = b; q < e; q++)
        if (_tmp[q] > 1) _phi[q] *= _tmp[q] - 1;
}

void segmented_mu(int64_t *mu, int64_t b, int64_t e, const int *p, int m) {
    int64_t *_mu = mu - b;
    if (b == 0) _mu[b++] = 0;
    for (int64_t q = b; q < e; q++)
        _mu[q] = 1;
    for (int i = 0; i < m; i++) {
        for (int64_t q = multiple<int64_t>(p[i], b); q < e; q += p[i]) {
            _mu[q] *= -p[i];
        }
        int64_t p2 = isq(p[i]);
        for (int64_t q = multiple<int64_t>(p2, b); q < e; q += p2) {
            _mu[q] = 0;
        }
    }
    // correction for a large prime factor (p > sqrt(e))
    for (int64_t q = b; q < e; q++) {
        if (_mu[q] < 0 && _mu[q] != -q) _mu[q] = +q;
        if (_mu[q] > 0 && _mu[q] != +q) _mu[q] = -q;
    }
    // normalize to +/-1
    for (int64_t q = b; q < e; q++) {
        if (_mu[q] < 0) _mu[q] = -1;
        if (_mu[q] > 0) _mu[q] = +1;
    }
}

void divisor_sigma0(int *ds0, int n) {
    for (int i = 1; i < n; i++)
        ds0[i] = 0;
    for (int i = 1; i < n; i++)
        for (int j = i; j < n; j += i)
            ds0[j] += 1;
}

void divisor_sigma1(int64_t *ds1, int n) {
    for (int i = 1; i < n; i++)
        ds1[i] = 0;
    for (int i = 1; i < n; i++)
        for (int j = i; j < n; j += i)
            ds1[j] += i;
}

void factor(int *bpf, int n, const int *p, int m) {
    bpf[0] = 0, bpf[1] = 1;
    for (int i = 0; i < m; i++)
        for (int j = p[i]; j < n; j += p[i])
            bpf[j] = p[i];
}

void factor_integer(std::vector<std::pair<int, int>> &vf, int n, const int *pf) {
    while (n > 1) {
        int p = pf[n], e = 0;
        while (n % p == 0) {
            n /= p, e++;
        }
        vf.push_back({ p, e });
    }
}

void factor_integer(std::vector<std::pair<int, int>> &vf, std::vector<int> vn, const int *pf) {
    for (auto &n : vn) {
        while (n > 1) {
            int p = pf[n], e = 0;
            for (auto &m : vn) {
                while (m % p == 0) {
                    m /= p, e++;
                }
            }
            vf.push_back({ p, e });
        }
    }
}

} // math
} // altruct
