#pragma once

#include "structure/math/polynom.h"
#include "algorithm/math/recurrence.h"

#include <vector>

namespace altruct {
namespace math {
		
/**
 * Searches for argument `x` within the monotonic interval `[b, e]` so that `p(x) = y`.
 *
 * Note, this method only works for floating point arguments.
 *
 * @param p - polynom (or arbitrary functor) to evaluate
 * @param b, e - range [b, e] to be searched. p(x) must be monotonic within that range
 * @param y - value for which argument has to be found
 * @param epsy - absolute tolerance for p(x)
 * @param epsx - absolute tolerance for x
 */
template<typename P, typename F>
F monotonic_search(const P &p, const F &b, const F& e, const F& y, const F& epsy = 0, const F& epsx = 0) {
	F lo = b, hi = e, mid = e;
	F val_b = p(b), val_e = p(e);
	if (-epsy <= val_b && val_b <= epsy) return b;
	if (-epsy <= val_e && val_e <= epsy) return e;
	int s = (val_b > val_e) ? -1 : +1;
	while (hi - lo > epsx) {
		mid = (lo + hi) / 2;
		if (mid == lo || mid == hi) return mid;
		F val = p(mid);
		if (-epsy <= val && val <= epsy) return mid;
		if ((s > 0) ? (val < y) : val > y) {
			lo = mid;
		}
		else {
			hi = mid;
		}
	}
	return mid;
}

/**
 * Finds polynom zeros.
 *
 * Note, this method only works for floating point arguments.
 *
 * @param p - polynom
 * @param inf - zeros are searched within [-inf, +inf] range only
 * @param epsy - absolute tolerance for p(x)
 * @param epsx - absolute tolerance for x
 */
template<typename T, typename F>
std::vector<F> find_zeros(const polynom<T>& p, const F& inf, const F& epsy = 0, const F& epsx = 0) {
	int l = p.deg();
	// derivations
	std::vector<polynom<T>> pd(l);
	pd[0] = p;
	for (int i = 1; i < l; i++)
		pd[i] = pd[i - 1].derivative();
	// derivations' zeros
	std::vector<F> z(l + 1);
	for (int i = l - 1; i >= 0; i--) {
		z[0] = -inf; z[l - i] = inf;
		for (int j = l - i; j >= 1; j--) {
			z[j] = monotonic_search<polynom<T>, F>(pd[i], z[j - 1], z[j], 0, epsy, epsx);
		}
	}
	// keep only actual zeros
	int n = 0;
	for (int i = 0; i <= l; i++) {
		F y = p(z[i]);
		if (-epsy <= y && y <= epsy) {
			z[n++] = z[i];
		}
	}
	z.resize(n);
	return z;
}

/**
 * Discrete integral of the polynom `p`
 *
 * s(n) = sum{p(k), {k, 1, n}}
 */
template<typename T>
polynom<T> polynom_sum(const polynom<T>& p) {
	polynom<T> s;
	std::vector<T> b = bernoulli_b<T>(p.deg());
	for (int m = p.deg(); m >= 0; m--) {
		T c = p[m] / (m + 1);
		if (c == 0) continue;
		for (int k = 0; k <= m; k++) {
			s[m + 1 - k] += c * b[k];
			c *= (m + 1 - k);
			c /= (k + 1);
		}
	}
	return s;
}

} // math
} // altruct
