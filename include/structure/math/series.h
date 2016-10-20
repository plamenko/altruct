#pragma once

#include "structure/math/polynom.h"

namespace altruct {
namespace math {

/**
 * Power series
 */
template<typename T, int N>
class series {
public:
	// s(x) = p(x) + O(x^N)
	polynom<T> p;

	series() : p() { p.reserve(N); }
	series(series&& rhs) : p(std::move(rhs.p)) { p.reserve(N); }
	series(const series& rhs) : p(rhs.p) { p.reserve(N); }
	series(polynom<T>&& rhs) : p(std::move(rhs)) { p.reserve(N); }
	series(const polynom<T>& rhs) : p(rhs) { p.reserve(N); }
	series(std::vector<T>&& c) : p(std::move(c)) { p.reserve(N); }
	series(const std::vector<T>& c) : p(c) { p.reserve(N); }
	template<typename It> series(It begin, It end) : p(begin, end) { p.reserve(N); }
	series(const T& c0) : p(c0) { p.reserve(N); }
	// construct from int, but only if T is not integral to avoid constructor clashing
	template <typename I = T, typename = std::enable_if_t<!std::is_integral<I>::value>>
	series(int c0) : p(c0) { p.reserve(N); } // to allow constructing from 0 and 1
	series(std::initializer_list<T> list) : p(list) { p.reserve(N); }
	series& operator = (const series& rhs) { p = rhs.p; return *this; }

	series& swap(series& rhs) { p.swap(rhs.p); return *this; }

	int size() const { return N; }
	const T& at(int index) const { return p.at(index); }
	const T& operator [] (int index) const { return p[index]; }
	T& operator [] (int index) { return p[index]; }

	bool operator == (const series& rhs) const { return p == rhs.p; }
	bool operator != (const series& rhs) const { return p != rhs.p; }
	bool operator <  (const series& rhs) const { return p <  rhs.p; }
	bool operator >  (const series& rhs) const { return p >  rhs.p; }
	bool operator <= (const series& rhs) const { return p <= rhs.p; }
	bool operator >= (const series& rhs) const { return p >= rhs.p; }
	
	series  operator +  (const series &rhs) const { series t(*this); t += rhs; return t; }
	series  operator -  (const series &rhs) const { series t(*this); t -= rhs; return t; }
	series  operator -  ()                  const { return series(-p); }
	series  operator *  (const series& rhs) const { series t(*this); t *= rhs; return t; }
	series  operator /  (const series& rhs) const { series t(*this); t /= rhs; return t; }
	
	series  operator *  (const T &val) const { series t(*this); t *= val; return t; }
	series  operator /  (const T &val) const { series t(*this); t /= val; return t; }
	
	series& operator += (const series &rhs) { p += rhs.p; return *this; }
	series& operator -= (const series &rhs) { p -= rhs.p; return *this; }
	series& operator *= (const series& rhs) { polynom<T>::mul(p, p, rhs.p, N - 1); return *this; }
	series& operator /= (const series& rhs) { return *this *= rhs.inverse(); }

	series& operator *= (const T &val) { p *= val; return *this; }
	series& operator /= (const T &val) { p /= val; return *this; }

	series derivative() const { return series(p.derivative()); }
	series integral(const T& c0 = polynom<T>::ZERO_COEFF) const { return series(p.integral(c0)); }

	// t(x) so that s(x) * t(x) == 1 + O(x^N); O(M(N) log N)
	series inverse() const {
		// ensure that p[0] is 1 before inverting
		if (p[0] == zeroT<T>::of(p[0])) return series();
		if (p[0] != identityT<T>::of(p[0])) return (*this / p[0]).inverse() / p[0];
		polynom<T> r{1}, t;
		for (int l = 1; l < N * 2; l *= 2) {
			int m = std::min(N - 1, l), k = l / 2 + 1;
			t.c.assign(p.c.begin(), p.c.begin() + m + 1);
			polynom<T>::mul(t, t, r, l + 1);
			t.c.erase(t.c.begin(), t.c.begin() + k);
			polynom<T>::mul(t, t, r, l - k);
			for (int i = m; i >= k; i--) {
				r[i] = -t[i - k];
			}
		}
		return series(r);
	}

	// series expansion of natural logarithm of s(x)
	series ln(const T& c0 = polynom<T>::ZERO_COEFF) const {
		series s = (derivative() / *this).integral(c0);
		s.p.c.resize(N);
		return s;
	}

	// series expansion of exp(a*x) = Sum[a^n * x^n / n!, n]
	static series exp(const T& a) {
		series s;
		s[0] = id_coeff();
		for (int i = 1; i < N; i++) {
			s[i] = s[i - 1] * a;
		}
		return s.make_exponential();
	}

	// converts the ordinary generating function to the exponential
	// by dividing coefficient of each x^n by n!
	series make_exponential() const {
		series s = *this;
		T fact = id_coeff();
		for (int i = 1; i < N; i++) {
			fact *= i;
		}
		T ifact = id_coeff() / fact;
		for (int i = N - 1; i > 0; i--) {
			s[i] *= ifact;
			ifact *= i;
		}
		return s;
	}

	// converts the exponential generating function to the ordinary
	// by multiplying coefficient of each x^n by n!
	series make_ordinary() const {
		series s = *this;
		T fact = id_coeff();
		for (int i = 1; i < N; i++) {
			fact *= i;
			s[i] *= fact;
		}
		return s;
	}

	// identity coefficient
	static T id_coeff() {
		return identityT<T>::of(polynom<T>::ZERO_COEFF);
	}
};

template<typename T, int N>
struct identityT<series<T, N>> {
	static series<T, N> of(const series<T, N>& s) {
		return series<T, N>(identityT<polynom<T>>::of(s.p));
	}
};

template<typename T, int N>
struct zeroT<series<T, N>> {
	static series<T, N> of(const series<T, N>& s) {
		return series<T, N>(zeroT<polynom<T>>::of(s.p));
	}
};

} // math
} // altruct
