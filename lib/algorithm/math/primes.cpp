#include "algorithm/math/primes.h"

#include <cmath>
#include <stdint.h>

namespace altruct {
namespace math {

typedef long long ll;

int primes(int *p, int *q, int n) {
	// if q not specified, use p for both
	if (!q) q = p;
	// initialize q
	q[0] = q[1] = 0;
	for (int i = 2; i < n; i++)
		q[i] = 1;
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
	return m;
}

void prime_pi(int *pi, int n, int *p, int m) {
	for (int i = 0, l = 0; i < n; i++) {
		if (i == p[l]) l++;
		pi[i] = l;
	}
}

void euler_phi(int *phi, int n, int *p, int m) {
	for (int i = 0; i < n; i++)
		phi[i] = i;
	for (int i = 0; i < m; i++)
		for (int j = p[i]; j < n; j += p[i])
			phi[j] = phi[j] / p[i] * (p[i] - 1);
}

int moebius_mu(int *mu, int n, int *p, int m) {
	int i, j;
	mu[0] = 0;
	for (i = 1; i < n; i++)
		mu[i] = 1;
	for (i = 0; i < m; i++)
	for (j = 0; j < n; j += p[i])
		mu[j] = -mu[j];
	for (i = 0; i < m; i++) {
		int p2 = p[i] * p[i]; if (p2 >= n) break;
		for (j = 0; j < n; j += p2)
			mu[j] = 0;
	}
	return 0;
}

void segmented_phi(ll *phi, ll *tmp, ll b, ll e, int *p, int m) {
	ll *_phi = phi - b, *_tmp = tmp - b;
	for (ll q = b; q < e; q++)
		_phi[q] = 1, _tmp[q] = q;
	for (int i = 0; i < m; i++) {
		for (ll q = div_ceil(b, p[i]) * p[i]; q < e; q += p[i]) {
			_phi[q] *= p[i] - 1, _tmp[q] /= p[i];
			while (_tmp[q] % p[i] == 0)
				_phi[q] *= p[i], _tmp[q] /= p[i];
		}
	}
	// correction for a large prime factor (p > sqrt(e))
	for (ll q = b; q < e; q++)
		if (_tmp[q] > 1) _phi[q] *= _tmp[q] - 1;
}

void segmented_mu(ll *mu, ll b, ll e, int *p, int m) {
	ll *_mu = mu - b;
	for (ll q = b; q < e; q++)
		_mu[q] = 1;
	for (int i = 0; i < m; i++) {
		for (ll q = div_ceil(b, p[i]) * p[i]; q < e; q += p[i]) {
			_mu[q] *= -p[i];
		}
		ll p2 = isq(p[i]);
		for (ll q = div_ceil(b, p2) * p2; q < e; q += p2) {
			_mu[q] = 0;
		}
	}
	// correction for a large prime factor (p > sqrt(e))
	for (ll q = b; q < e; q++) {
		if (_mu[q] < 0 && _mu[q] != -q) _mu[q] = +q;
		if (_mu[q] > 0 && _mu[q] != +q) _mu[q] = -q;
	}
	// normalize to +/-1
	for (ll q = b; q < e; q++) {
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

void divisor_sigma1(ll *ds1, int n) {
	for (int i = 1; i < n; i++)
		ds1[i] = 0;
	for (int i = 1; i < n; i++)
		for (int j = i; j < n; j += i)
			ds1[j] += i;
}

void factor(int *bpf, int n, int *p, int m) {
	bpf[0] = 0, bpf[1] = 1;
	for (int i = 0; i < m; i++)
		for (int j = p[i]; j < n; j += p[i])
			bpf[j] = p[i];
}

void factor_integer(std::vector<std::pair<int, int>> &vf, int n, int *pf) {
	while (n > 1) {
		int p = pf[n];
		int e = 0; while (n % p == 0) n /= p, e++;
		vf.push_back({ p, e });
	}
}

} // math
} // altruct
