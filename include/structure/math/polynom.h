#pragma once

#include "algorithm/math/base.h"
#include <limits>
#include <type_traits>
#include <vector>
#include <algorithm>

namespace altruct {
namespace math {

template<typename T> struct polynom_mul;

/**
 * Polynomial with coefficients in T.
 */
template<typename T>
class polynom {
public:

	template<typename It>
	static T make_zero(It begin, It end) { return (begin != end) ? zeroT<T>::of(*begin) : T(0); }
	static T make_zero(const T& c) { return zeroT<T>::of(c); }

	T ZERO_COEFF;
	
	// p(x) = sum{c[i] * x^i}
	std::vector<T> c;

	polynom(const T& c0 = T(0)) : ZERO_COEFF(make_zero(c0)) { c.push_back(c0); }
	// construct from int, but only if T is not integral to avoid constructor clashing
	template <typename I = T, typename = std::enable_if_t<!std::is_integral<I>::value>>
	polynom(int c0) : ZERO_COEFF(make_zero(c0)) { c.push_back(c0); } // to allow constructing from 0 and 1
	polynom(std::vector<T>&& rhs) : ZERO_COEFF(make_zero(rhs.begin(), rhs.end())), c(std::move(rhs)) {}
	polynom(const std::vector<T>& rhs) : ZERO_COEFF(make_zero(rhs.begin(), rhs.end())), c(rhs) {}
	template<typename It> polynom(It begin, It end) : ZERO_COEFF(make_zero(begin, end)), c(begin, end) {}
	polynom(std::initializer_list<T> list) : ZERO_COEFF(make_zero(list.begin(), list.end())), c(list) {}

	polynom& swap(polynom &rhs) { std::swap(ZERO_COEFF, rhs.ZERO_COEFF); c.swap(rhs.c); return *this; }
	polynom& shrink_to_fit() { c.resize(deg() + 1, ZERO_COEFF); return *this; }
	polynom& reserve(int sz) { if (sz > size()) c.resize(sz, ZERO_COEFF); return *this; }
	polynom& resize(int sz) { if (sz != size()) c.resize(sz, ZERO_COEFF); return *this; }
	polynom& resize(int sz, const T& _ZERO_COEFF) { ZERO_COEFF = _ZERO_COEFF; return resize(sz); }

	int size() const { return (int)c.size(); }
	const T& at(int index) const { return (0 <= index && index < size()) ? c[index] : ZERO_COEFF; }
	const T& operator [] (int index) const { return at(index); }
	T& operator [] (int index) { reserve(index + 1); return c[index]; }
	int deg() const { for (int i = size() - 1; i > 0; i--) if (!(c[i] == ZERO_COEFF)) return i; return 0; }
	const T& leading_coeff() const { return at(deg()); }
	bool is_power() const { for (int i = deg() - 1; i >= 0; i--) if (!(c[i] == ZERO_COEFF)) return false; return leading_coeff() == identityT<T>::of(ZERO_COEFF); }

	// compares p1 and p2; O(l1 + l2)
	static int cmp(const polynom &p1, const polynom &p2) {
		int l1 = p1.deg(), l2 = p2.deg(); int l = std::max(l1, l2);
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
		pr.resize(lr + 1, p1.ZERO_COEFF);
		for (int i = 0; i <= lr; i++) {
			pr[i] = -p1[i];
		}
	}

	// pr = p1 + p2; O(l1 + l2)
	// it is allowed for `p1`, `p2` and `pr` to be the same instance
	static void add(polynom &pr, const polynom &p1, const polynom &p2) {
		int l1 = p1.deg(), l2 = p2.deg(); int lr = std::max(l1, l2);
		pr.resize(lr + 1, p1.ZERO_COEFF);
		for (int i = 0; i <= lr; i++) {
			pr[i] = p1[i] + p2[i];
		}
	}
	
	// pr = p1 - p2; O(l1 + l2)
	// it is allowed for `p1`, `p2` and `pr` to be the same instance
	static void sub(polynom &pr, const polynom &p1, const polynom &p2) {
		int l1 = p1.deg(), l2 = p2.deg(); int lr = std::max(l1, l2);
		pr.resize(lr + 1, p1.ZERO_COEFF);
		for (int i = 0; i <= lr; i++) {
			pr[i] = p1[i] - p2[i];
		}
	}

	// pr[lm + 1 : lr] = 0; O(l)
	static void _zero(T* pr, int lm, int lr, const T& ZERO_COEFF) {
		for (int i = lm + 1; i <= lr; i++) pr[i] = ZERO_COEFF;
	}

	// pr += p2; O(l2)
	// it is allowed for `pr` and `p2` to be the same instance
	static void _add_to(T* pr, const T* p2, int l2) {
		for (int i = 0; i <= l2; i++) pr[i] += p2[i];
	}

	// pr -= p2; O(l2)
	// it is allowed for `pr` and `p2` to be the same instance
	static void _sub_from(T* pr, const T* p2, int l2) {
		for (int i = 0; i <= l2; i++) pr[i] -= p2[i];
	}

	// pr = p1 * p2; O(l1 * l2)
	// it is allowed for `p1`, `p2` and `pr` to be the same instance
	static void _mul_long(T* pr, int lr, const T* p1, int l1, const T* p2, int l2) {
		auto ZERO_COEFF = make_zero(*p1);
		for (int i = lr; i >= 0; i--) {
			T r = ZERO_COEFF;
			int jmax = std::min(i, l1);
			int jmin = std::max(0, i - l2);
			for (int j = jmax; j >= jmin; j--) {
				r += p1[j] * p2[i - j];
			}
			pr[i] = r;
		}
	}

	// pr = p1 * p2; O(lr ^ 1.59); or more accurate: O(l1 * l2 ^ 0.59)
	// `0 <= l2 <= l1 <= lr <= l1 + l2` must hold
	// it is allowed for `p1`, `p2` and `pr` to be the same instance
	static void _mul_karatsuba(T* pr, int lr, const T* p1, int l1, const T* p2, int l2) {
		auto ZERO_COEFF = make_zero(*p1);
		int k = l1 / 2 + 1; // k > l1 - k >= 0
		if (l2 == 0) {
			for (int i = lr; i >= 0; i--) pr[i] = p1[i] * p2[0];
		} else if (l2 < k) {
			std::vector<T> MM(lr - k + 1, ZERO_COEFF);
			_mul(MM.data(), lr - k, p1 + k, l1 - k, p2, l2);
			_mul(pr, std::min(lr, l2 + k - 1), p1, k - 1, p2, l2);
			_zero(pr, l2 + k - 1, lr, ZERO_COEFF);
			_add_to(pr + k, MM.data(), lr - k);
		} else {
			std::vector<T> S1(p1, p1 + k);
			_add_to(S1.data(), p1 + k, l1 - k);
			std::vector<T> S2(p2, p2 + k);
			_add_to(S2.data(), p2 + k, l2 - k);
			int mm_l = std::min(lr - k, k - 1 + k - 1);
			std::vector<T> MM(mm_l + 1);
			_mul(MM.data(), mm_l, S1.data(), k - 1, S2.data(), k - 1);
			int hh_l = std::min(mm_l, l1 - k + l2 - k);
			std::vector<T> HH(hh_l + 1, ZERO_COEFF);
			_mul(HH.data(), hh_l, p1 + k, l1 - k, p2 + k, l2 - k);
			_mul(pr, k - 1 + k - 1, p1, k - 1, p2, k - 1);
			_zero(pr, k - 1 + k - 1, lr, ZERO_COEFF);
			_sub_from(MM.data(), pr, std::min(mm_l, k - 1 + k - 1));
			_sub_from(MM.data(), HH.data(), hh_l);
			_add_to(pr + k, MM.data(), mm_l);
			_add_to(pr + k + k, HH.data(), lr - k - k);
		}
	}

	// ensures `l2 <= l1 <= lr <= l1 + l2` and delegates to `polynom_mul<T>::impl`
	// it is allowed for `p1`, `p2` and `pr` to be the same instance
	static void _mul(T* pr, int lr, const T* p1, int l1, const T* p2, int l2) {
		if (l2 > l1) return _mul(pr, lr, p2, l2, p1, l1);   // ensure `l2 <= l1`
		l1 = std::min(l1, lr); l2 = std::min(l2, lr);       // ensure `l2 <= l1 <= lr`
		_zero(pr, l1 + l2, lr, make_zero(*p1));
		lr = std::min(lr, l1 + l2);                         // ensure `lr <= l1 + l2`
		polynom_mul<T>::impl(pr, lr, p1, l1, p2, l2);
	}

	// pr = p1 * p2;
	// it is allowed for `p1`, `p2` and `pr` to be the same instance
	// @param lr - the required degree of the resulting polynomial;
	//             if -1, the result will be of degree l1 + l2
	static void mul(polynom &pr, const polynom &p1, const polynom &p2, int lr = -1) {
		int l1 = p1.deg(), l2 = p2.deg(); if (lr < 0) lr = l1 + l2;
		pr.resize(lr + 1, p1.ZERO_COEFF);
		if (p1.size() == 0 || p2.size() == 0) {
			_zero(pr.c.data(), -1, lr, p1.ZERO_COEFF);
		} else {
			_mul(pr.c.data(), lr, p1.c.data(), l1, p2.c.data(), l2);
		}
	}

	// pr = p1 % p2 | p1 / p2; O((l1 - lm) * lm)
	// it is allowed for `p1` and `pr` to be the same instance
	// `pr` and `pm` must not be the same instance
	static void quot_rem(polynom &pr, const polynom &p1, const polynom &pm) {
		int l1 = p1.deg(), lm = pm.deg(); int lr = l1 - lm;
		pr = p1;
		if (lr < 0 || pm.is_power()) return;
		for (int i = l1; i >= lm; i--) {
			T s = pr[i] /= pm[lm]; if (s == p1.ZERO_COEFF) continue;
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
		pr.resize(lr + 1);
	}

	// pr = p1 % pm; O((l1 - lm) * lm)
	// it is allowed for `p1` and `pr` to be the same instance
	// `pr` and `pm` must not be the same instance
	static void mod(polynom &pr, const polynom &p1, const polynom &pm) {
		int l1 = p1.deg(), lm = pm.deg(); int lr = lm - 1;
		quot_rem(pr, p1, pm);
		if (lr < l1) pr.resize(lr + 1);
	}

	// pr = p1 * s; O(l1)
	// it is allowed for `p1` and `pr` to be the same instance
	static void mul(polynom &pr, const polynom &p1, const T &s) {
		int lr = p1.deg();
		pr.resize(lr + 1, p1.ZERO_COEFF);
		for (int i = 0; i <= lr; i++) {
			pr[i] = p1[i] * s;
		}
	}

	// pr = p1 / s; O(l1)
	// it is allowed for `p1` and `pr` to be the same instance
	static void div(polynom &pr, const polynom &p1, const T &s) {
		int lr = p1.deg();
		pr.resize(lr + 1, p1.ZERO_COEFF);
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
		polynom r(ZERO_COEFF);
		for (int i = deg(); i > 0; i--) {
			r[i - 1] = c[i] * i;
		}
		return r;
	}
	
	polynom integral() const { return integral(ZERO_COEFF); }
	polynom integral(const T& c0) const {
		polynom r(c0);
		for (int i = deg(); i >= 0; i--) {
			r[i + 1] = c[i] / (i + 1);
		}
		return r;
	}
};

/**
 * `polynom<T>` multiplication implementation.
 *
 * Specialize this template for a custom or tweaked implementation.
 *
 * If you need to call multiplication recursively, don't call
 * `impl` directly, but call `polynom<T>::_mul` instead as it
 * ensures the invariants before delegating to this `impl`.
 *
 * You may also call one of the already provided implementations:
 * `polynom<T>::_mul_long` or `polynom<T>::_mul_karatsuba`,
 * or the utility methods: `_add_to`, `_sub_from` and `_zero`.
 */
template<typename T>
struct polynom_mul {
	// @param pX - The polynomials to perform multiplication on: `pr = p1 * p2`.
	//             It is allowed for `p1`, `p2` and `pr` to be the same instance.
	// @param lX - The lengths of the polynomials: `0 <= l2 <= l1 <= lr <= l1 + l2`.
	//             Coefficients of `pr` in the range [0, lr] inclusive must be set.
	//             Truncate or pad with 0 if necessary.
	static void impl(T* pr, int lr, const T* p1, int l1, const T* p2, int l2) {
		if (l2 < 15 || int64_t(l1) * l2 < 300) {  // just l2 < 16 ???
			polynom<T>::_mul_long(pr, lr, p1, l1, p2, l2);
		} else {
			polynom<T>::_mul_karatsuba(pr, lr, p1, l1, p2, l2);
		}
	}
};

template<typename T>
struct identityT<polynom<T>> {
	static polynom<T> of(const polynom<T>& p) {
		return polynom<T>(identityT<T>::of(p.ZERO_COEFF));
	}
};

template<typename T>
struct zeroT<polynom<T>> {
	static polynom<T> of(const polynom<T>& p) {
		return polynom<T>(p.ZERO_COEFF);
	}
};

} // math
} // altruct
