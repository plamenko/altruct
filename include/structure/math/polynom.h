#pragma once

#include "algorithm/math/fft.h"

#include <type_traits>
#include <vector>

namespace altruct {
namespace math {

template<typename T>
class polynom {
public:

	// zero coefficient used for const reference return
	static T ZERO_COEFF;
	
	// FFT root and order; if set, FFT will be used for multiplication
	static T FFT_ROOT;
	static int FFT_ORDER;

	// p(x) = sum{c[i] * x^i}
	std::vector<T> c;

	polynom() {}
	polynom(const polynom& rhs) : c(rhs.c) {}
	polynom(const std::vector<T>& c) : c(c) {}
	template<typename It> polynom(It begin, It end) : c(begin, end) {}
	// construct from T, but only if T is not integral to avoid clashing with `polynom(int c0)`
	template <typename = std::enable_if_t<!std::is_integral<T>::value>>
	polynom(const T& c0) { c.push_back(c0); }
	polynom(int c0) { c.push_back(c0); } // to allow constructing from 0 and 1
	polynom(std::initializer_list<T> list) : c(list) {}
	
	polynom& swap(polynom &rhs) { c.swap(rhs.c); return *this; }
	polynom& shrink_to_fit() { c.resize(deg() + 1); return *this; }
	polynom& reserve(int sz) { if (sz > size()) c.resize(sz); return *this; }

	int size() const { return (int)c.size(); }
	const T& at(int index) const { return (0 <= index && index < size()) ? c[index] : ZERO_COEFF; }
	const T& operator [] (int index) const { return at(index); }
	T& operator [] (int index) { reserve(index + 1); return c[index]; }
	int deg() const { for (int i = size() - 1; i > 0; i--) if (c[i] != 0) return i; return 0; }
	const T& leading_coeff() const { return at(deg()); }
	bool is_power() const { for (int i = deg() - 1; i >= 0; i--) if (c[i] != 0) return false; return leading_coeff() == 1; }

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
		pr.c.resize(lr + 1);
		for (int i = 0; i <= lr; i++)
			pr[i] = -p1[i];
	}

	// pr = p1 + p2; O(l1 + l2)
	// it is allowed for `p1`, `p2` and `pr` to be the same instance
	static void add(polynom &pr, const polynom &p1, const polynom &p2) {
		int l1 = p1.deg(), l2 = p2.deg(); int lr = max(l1, l2);
		pr.c.resize(lr + 1); 
		for (int i = 0; i <= lr; i++)
			pr[i] = p1[i] + p2[i];
	}
	
	// pr = p1 - p2; O(l1 + l2)
	// it is allowed for `p1`, `p2` and `pr` to be the same instance
	static void sub(polynom &pr, const polynom &p1, const polynom &p2) {
		int l1 = p1.deg(), l2 = p2.deg(); int lr = max(l1, l2);
		pr.c.resize(lr + 1); 
		for (int i = 0; i <= lr; i++)
			pr[i] = p1[i] - p2[i];
	}
	
	// pr = p1 * p2; O(lr log lr)
	// it is allowed for `p1`, `p2` and `pr` to be the same instance
	static void mul_fft(polynom &pr, const polynom &p1, const polynom &p2) {
		int sizeR = p1.deg() + p2.deg() + 1;
		int sizeFFT = 1; while (sizeFFT < sizeR) sizeFFT *= 2;
		std::vector<T> temp1(p1.c); temp1.resize(sizeFFT);
		std::vector<T> temp2(p2.c); temp2.resize(sizeFFT);
		pr.c.assign(sizeFFT, 0);
		fft_convolve(&pr.c[0], &temp1[0], &temp2[0], sizeFFT, FFT_ROOT, FFT_ORDER);
		pr.c.resize(sizeR);
	}
	
	// pr = p1 * p2; O(l1 * l2)
	// it is allowed for `p1`, `p2` and `pr` to be the same instance
	static void mul_long(polynom &pr, const polynom &p1, const polynom &p2) {
		int l1 = p1.deg(), l2 = p2.deg(); int lr = l1 + l2;
		pr.c.resize(lr + 1);
		for (int i = lr; i >= 0; i--) {
			T r = 0;
			for (int j = min(i, l1); j >= max(0, i - l2); j--)
				r += p1[j] * p2[i - j];
			pr[i] = r;
		}
	}

	// it is allowed for `p1`, `p2` and `pr` to be the same instance
	static void mul(polynom &pr, const polynom &p1, const polynom &p2) {
		int sizeR = p1.deg() + p2.deg() + 1;
		long long cost = (long long) p1.deg() * p2.deg();
		if (FFT_ORDER >= sizeR && cost > 10000) {
			mul_fft(pr, p1, p2);
		} else {
			mul_long(pr, p1, p2);
		}
	}

	// pr = p1 % p2 | p1 / p2; O((l1 - lm) * lm)
	// it is allowed for `p1` and `pr` to be the same instance
	// `pr` and `pm` must not be the same instance
	static void quot_rem(polynom &pr, const polynom &p1, const polynom &pm) {
		int l1 = p1.deg(), lm = pm.deg(); int lr = l1 - lm;
		pr.c = p1.c; pr.c.resize(l1 + 1);
		if (lr < 0 || pm.is_power()) return;
		for (int i = l1; i >= lm; i--) {
			T s = pr[i] /= pm[lm]; if (s == 0) continue;
			for (int j = 1; j <= lm; j++)
				pr[i - j] = pr[i - j] - s * pm[lm - j];
		}
	}

	// pr = p1 / p2; O((l1 - lm) * lm)
	// it is allowed for `p1` and `pr` to be the same instance
	// `pr` and `pm` must not be the same instance
	static void div(polynom &pr, const polynom &p1, const polynom &pm) {
		int l1 = p1.deg(), lm = pm.deg(); int lr = l1 - lm;
		if (lr < 0) { pr.c.clear(); return; }
		quot_rem(pr, p1, pm);
		for (int i = 0; i <= lr; i++)
			pr[i] = pr[i + lm];
		pr.c.resize(lr + 1);
	}

	// pr = p1 % pm; O((l1 - lm) * lm)
	// it is allowed for `p1` and `pr` to be the same instance
	// `pr` and `pm` must not be the same instance
	static void mod(polynom &pr, const polynom &p1, const polynom &pm) {
		int l1 = p1.deg(), lm = pm.deg(); int lr = lm - 1;
		quot_rem(pr, p1, pm);
		if (lr < l1) pr.c.resize(lr + 1);
	}

	// pr = p1 * s; O(l1)
	// it is allowed for `p1` and `pr` to be the same instance
	static void mul(polynom &pr, const polynom &p1, const T &s) {
		int lr = p1.deg();
		pr.c.resize(lr + 1);
		for (int i = 0; i <= lr; i++)
			pr[i] = p1[i] * s;
	}

	// pr = p1 / s; O(l1)
	// it is allowed for `p1` and `pr` to be the same instance
	static void div(polynom &pr, const polynom &p1, const T &s) {
		int lr = p1.deg();
		pr.c.resize(lr + 1);
		for (int i = 0; i <= lr; i++)
			pr[i] = p1[i] / s;
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
		A r = 0;
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
};

template<typename T> T polynom<T>::ZERO_COEFF = 0;
template<typename T> T polynom<T>::FFT_ROOT = 0;
template<typename T> int polynom<T>::FFT_ORDER = 0;

} // math
} // altruct
