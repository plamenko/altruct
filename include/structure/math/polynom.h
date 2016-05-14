#pragma once

#include "algorithm/math/fft.h"

#include <type_traits>
#include <vector>

namespace altruct {
namespace math {

template<typename T>
class polynom {
public:

	// zero coefficient
	// TODO: perhaps we should keep this as c.back() instead of a static member.
	// Right now we have issues with polynom<moduloX>.
	static T ZERO_COEFF;
	
	// FFT root and order; if set, FFT will be used for multiplication
	static T FFT_ROOT;
	static int FFT_ORDER;

	// p(x) = sum{c[i] * x^i}
	std::vector<T> c;

	polynom() {}
	polynom(polynom&& rhs) : c(std::move(rhs.c)) {}
	polynom(const polynom& rhs) : c(rhs.c) {}
	polynom(std::vector<T>&& rhs) : c(std::move(rhs)) {}
	polynom(const std::vector<T>& rhs) : c(rhs) {}
	template<typename It> polynom(It begin, It end) : c(begin, end) {}
	polynom(const T& c0) { c.push_back(c0); }
	// construct from int, but only if T is not integral to avoid constructor clashing
	template <typename = std::enable_if_t<!std::is_integral<T>::value>>
	polynom(int c0) { c.push_back(c0); } // to allow constructing from 0 and 1
	polynom(std::initializer_list<T> list) : c(list) {}

	polynom& swap(polynom &rhs) { c.swap(rhs.c); return *this; }
	polynom& shrink_to_fit() { c.resize(deg() + 1, ZERO_COEFF); return *this; }
	polynom& reserve(int sz) { if (sz > size()) c.resize(sz, ZERO_COEFF); return *this; }

	int size() const { return (int)c.size(); }
	const T& at(int index) const { return (0 <= index && index < size()) ? c[index] : ZERO_COEFF; }
	const T& operator [] (int index) const { return at(index); }
	T& operator [] (int index) { reserve(index + 1); return c[index]; }
	int deg() const { for (int i = size() - 1; i > 0; i--) if (!(c[i] == ZERO_COEFF)) return i; return 0; }
	const T& leading_coeff() const { return at(deg()); }
	bool is_power() const { for (int i = deg() - 1; i >= 0; i--) if (!(c[i] == ZERO_COEFF)) return false; return leading_coeff() == identityT<T>::of(ZERO_COEFF); }

	// compares p1 and p2; O(l1 + l2)
	static int cmp(const polynom &p1, const polynom &p2) {
		int l1 = p1.deg(), l2 = p2.deg(); int l = max(l1, l2);
		for (int i = l; i >= 0; i--) {
			if (p1[i] < p2[i]) return -1;
			if (p2[i] < p1[i]) return +1;
		}
		return 0;
	}
	
	// pr = -p1; O(l1)
	// it is allowed for `p1`, and `pr` to be the same instance
	static void neg(polynom &pr, const polynom &p1) {
		int lr = p1.deg();
		pr.c.resize(lr + 1, ZERO_COEFF);
		for (int i = 0; i <= lr; i++) {
			pr[i] = -p1[i];
		}
	}

	// pr = p1 + p2; O(l1 + l2)
	// it is allowed for `p1`, `p2` and `pr` to be the same instance
	static void add(polynom &pr, const polynom &p1, const polynom &p2) {
		int l1 = p1.deg(), l2 = p2.deg(); int lr = max(l1, l2);
		pr.c.resize(lr + 1, ZERO_COEFF);
		for (int i = 0; i <= lr; i++) {
			pr[i] = p1[i] + p2[i];
		}
	}
	
	// pr = p1 - p2; O(l1 + l2)
	// it is allowed for `p1`, `p2` and `pr` to be the same instance
	static void sub(polynom &pr, const polynom &p1, const polynom &p2) {
		int l1 = p1.deg(), l2 = p2.deg(); int lr = max(l1, l2);
		pr.c.resize(lr + 1, ZERO_COEFF);
		for (int i = 0; i <= lr; i++) {
			pr[i] = p1[i] - p2[i];
		}
	}
	
	// pr = p1 * p2; O(lr log lr)
	// it is allowed for `p1`, `p2` and `pr` to be the same instance
	static void mul_fft(polynom &pr, const polynom &p1, const polynom &p2, int lr = -1) {
		int l1 = p1.deg(), l2 = p2.deg(); if (lr < 0) lr = l1 + l2;
		int sizeFFT = 1; while (sizeFFT <= lr) sizeFFT *= 2;
		std::vector<T> temp1(p1.c); temp1.resize(sizeFFT);
		std::vector<T> temp2(p2.c); temp2.resize(sizeFFT);
		pr.c.assign(sizeFFT, ZERO_COEFF);
		fft_cyclic_convolution(&pr.c[0], &temp1[0], &temp2[0], sizeFFT, FFT_ROOT, FFT_ORDER);
		pr.c.resize(lr + 1, ZERO_COEFF);
	}
	
	// pr = p1 * p2; O(lr ^ 1.59); or more accurate: O(l1 * l2 ^ 0.59) for l1 >= l2
	// it is allowed for `p1`, `p2` and `pr` to be the same instance
	static void mul_karatsuba(polynom &pr, const polynom &p1, const polynom &p2, int lr = -1) {
		int l1 = p1.deg(), l2 = p2.deg(); if (lr < 0) lr = l1 + l2;
		if (l1 < l2) return mul_karatsuba(pr, p2, p1, lr); // ensure l1 >= l2
		l1 = min(l1, lr), l2 = min(l2, lr); int k = l1 / 2 + 1;
		if (l2 == 0) {
			pr.c.resize(lr + 1, ZERO_COEFF);
			for (int i = lr; i > l1; i--) pr[i] = ZERO_COEFF;
			for (int i = l1; i >= 0; i--) pr[i] = p1[i] * p2[0];
		} else if (l2 < k) {
			polynom lo(p1.c.begin(), p1.c.begin() + k);
			polynom hi(p1.c.begin() + k, p1.c.begin() + l1 + 1);
			mul(lo, lo, p2, l2 + k - 1);
			mul(hi, hi, p2, min(lr, l2 + l1) - k);
			pr.c.assign(lr + 1, ZERO_COEFF);
			for (int i = lo.deg(); i >= 0; i--) pr[i] = lo[i];
			for (int i = hi.deg(); i >= 0; i--) pr[k + i] += hi[i];
		} else {
			polynom lo1(p1.c.begin(), p1.c.begin() + k);
			polynom hi1(p1.c.begin() + k, p1.c.begin() + l1 + 1);
			polynom lo2(p2.c.begin(), p2.c.begin() + k);
			polynom hi2(p2.c.begin() + k, p2.c.begin() + l2 + 1);
			polynom m0; mul(m0, lo1, lo2);
			lo1 += hi1; lo2 += hi2;
			polynom m1; mul(m1, lo1, lo2, min(lr - k, l1));
			polynom m2; mul(m2, hi1, hi2, min(lr - k, l2));
			pr.c.assign(lr + 1, ZERO_COEFF);
			for (int i = m0.deg(); i >= 0; i--) pr[i] = m0[i];
			for (int i = min(lr - k, m0.deg()); i >= 0; i--) pr[k + i] -= m0[i];
			for (int i = min(lr - k, m1.deg()); i >= 0; i--) pr[k + i] += m1[i];
			for (int i = min(lr - k, m2.deg()); i >= 0; i--) pr[k + i] -= m2[i];
			for (int i = min(lr - k - k, m2.deg()); i >= 0; i--) pr[k + k + i] += m2[i];
		}
	}

	// pr = p1 * p2; O(l1 * l2)
	// it is allowed for `p1`, `p2` and `pr` to be the same instance
	static void mul_long(polynom &pr, const polynom &p1, const polynom &p2, int lr = -1) {
		int l1 = p1.deg(), l2 = p2.deg(); if (lr < 0) lr = l1 + l2;
		pr.c.resize(lr + 1, ZERO_COEFF);
		for (int i = lr; i >= 0; i--) {
			T r = ZERO_COEFF;
			for (int j = min(i, l1); j >= max(0, i - l2); j--)
				r += p1[j] * p2[i - j];
			pr[i] = r;
		}
	}

	// it is allowed for `p1`, `p2` and `pr` to be the same instance
	static void mul(polynom &pr, const polynom &p1, const polynom &p2, int lr = -1) {
		int l1 = p1.deg(), l2 = p2.deg(); if (lr < 0) lr = l1 + l2;
		auto cost = min<long long>(l1, lr) * min(l2, lr);
		if (FFT_ORDER > lr && cost > 10000) {
			mul_fft(pr, p1, p2, lr);
		} else if (cost > 300) {
			mul_karatsuba(pr, p1, p2, lr);
		} else {
			mul_long(pr, p1, p2, lr);
		}
	}

	// pr = p1 % p2 | p1 / p2; O((l1 - lm) * lm)
	// it is allowed for `p1` and `pr` to be the same instance
	// `pr` and `pm` must not be the same instance
	static void quot_rem(polynom &pr, const polynom &p1, const polynom &pm) {
		int l1 = p1.deg(), lm = pm.deg(); int lr = l1 - lm;
		pr.c = p1.c; pr.c.resize(l1 + 1, ZERO_COEFF);
		if (lr < 0 || pm.is_power()) return;
		for (int i = l1; i >= lm; i--) {
			T s = pr[i] /= pm[lm]; if (s == ZERO_COEFF) continue;
			for (int j = 1; j <= lm; j++) {
				pr[i - j] = pr[i - j] - s * pm[lm - j];
			}
		}
	}

	// pr = p1 / p2; O((l1 - lm) * lm)
	// it is allowed for `p1` and `pr` to be the same instance
	// `pr` and `pm` must not be the same instance
	static void div(polynom &pr, const polynom &p1, const polynom &pm) {
		int l1 = p1.deg(), lm = pm.deg(); int lr = l1 - lm;
		if (lr < 0) { pr.c.clear(); return; }
		quot_rem(pr, p1, pm);
		for (int i = 0; i <= lr; i++) {
			pr[i] = pr[i + lm];
		}
		pr.c.resize(lr + 1, ZERO_COEFF);
	}

	// pr = p1 % pm; O((l1 - lm) * lm)
	// it is allowed for `p1` and `pr` to be the same instance
	// `pr` and `pm` must not be the same instance
	static void mod(polynom &pr, const polynom &p1, const polynom &pm) {
		int l1 = p1.deg(), lm = pm.deg(); int lr = lm - 1;
		quot_rem(pr, p1, pm);
		if (lr < l1) pr.c.resize(lr + 1, ZERO_COEFF);
	}

	// pr = p1 * s; O(l1)
	// it is allowed for `p1` and `pr` to be the same instance
	static void mul(polynom &pr, const polynom &p1, const T &s) {
		int lr = p1.deg();
		pr.c.resize(lr + 1, ZERO_COEFF);
		for (int i = 0; i <= lr; i++) {
			pr[i] = p1[i] * s;
		}
	}

	// pr = p1 / s; O(l1)
	// it is allowed for `p1` and `pr` to be the same instance
	static void div(polynom &pr, const polynom &p1, const T &s) {
		int lr = p1.deg();
		pr.c.resize(lr + 1, ZERO_COEFF);
		for (int i = 0; i <= lr; i++) {
			pr[i] = p1[i] / s;
		}
	}

	bool operator == (const polynom &rhs) const { return cmp(*this, rhs) == 0; }
	bool operator != (const polynom &rhs) const { return cmp(*this, rhs) != 0; }
	bool operator <  (const polynom &rhs) const { return cmp(*this, rhs) <  0; }
	bool operator >  (const polynom &rhs) const { return cmp(*this, rhs) >  0; }
	bool operator <= (const polynom &rhs) const { return cmp(*this, rhs) <= 0; }
	bool operator >= (const polynom &rhs) const { return cmp(*this, rhs) >= 0; }
	
	polynom  operator +  (const polynom &rhs) const { polynom t(*this); t += rhs; return t; }
	polynom  operator -  (const polynom &rhs) const { polynom t(*this); t -= rhs; return t; }
	polynom  operator -  ()                   const { polynom t(*this); neg(t, t); return t; }
	polynom  operator *  (const polynom &rhs) const { polynom t(*this); t *= rhs; return t; }
	polynom  operator /  (const polynom &rhs) const { polynom t(*this); t /= rhs; return t; }
	polynom  operator %  (const polynom &rhs) const { polynom t(*this); t %= rhs; return t; }

	polynom  operator *  (const T &val) const { polynom t(*this); t *= val; return t; }
	polynom  operator /  (const T &val) const { polynom t(*this); t /= val; return t; }

	polynom& operator += (const polynom &rhs) { add(*this, *this, rhs); return *this; }
	polynom& operator -= (const polynom &rhs) { sub(*this, *this, rhs); return *this; }
	polynom& operator *= (const polynom &rhs) { mul(*this, *this, rhs); return *this; }
	polynom& operator /= (const polynom &rhs) { div(*this, *this, rhs); return *this; }
	polynom& operator %= (const polynom &rhs) { mod(*this, *this, rhs); return *this; }
	
	polynom& operator *= (const T &val) { mul(*this, *this, val); return *this; }
	polynom& operator /= (const T &val) { div(*this, *this, val); return *this; }

	template<typename A>
	A operator () (const A& x) const { return eval<A>(x); }
	
	template<typename A>
	A eval(const A& x) const {
		A r = zeroT<A>::of(x);
		for (int i = deg(); i >= 0; i--) {
			r = r * x + A(c[i]);
		}
		return r;
	}

	polynom derivative() const {
		polynom r;
		for (int i = deg(); i > 0; i--) {
			r[i - 1] = c[i] * i;
		}
		return r;
	}
	
	polynom integral(const T& c0 = ZERO_COEFF) const {
		polynom r;
		r[0] = c0;
		for (int i = deg(); i >= 0; i--) {
			r[i + 1] = c[i] / (i + 1);
		}
		return r;
	}
};

template<typename T> T polynom<T>::ZERO_COEFF = 0;
template<typename T> T polynom<T>::FFT_ROOT = 0;
template<typename T> int polynom<T>::FFT_ORDER = 0;

template<typename T>
struct identityT<polynom<T>> {
	static polynom<T> of(const polynom<T>& p) {
		return polynom<T>(identityT<T>::of(polynom<T>::ZERO_COEFF));
	}
};

template<typename T>
struct zeroT<polynom<T>> {
	static polynom<T> of(const polynom<T>& p) {
		return polynom<T>(polynom<T>::ZERO_COEFF);
	}
};

} // math
} // altruct
