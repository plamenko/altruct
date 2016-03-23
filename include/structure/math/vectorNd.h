#pragma once

#include "algorithm/math/base.h"

#include <array>

namespace altruct {
namespace math {

/**
 * A vector of type T and fixed size N.
 */
template<typename T, int N>
class vectorNd {
public:
	std::array<T, N> a;

	vectorNd() { a.fill(0); }
	vectorNd(const std::array<T, N>& rhs) : a(rhs) { }
	vectorNd(const vectorNd& rhs) : a(rhs.a) { }
	// construct from T, but only if T is not integral to avoid clashing with `vectorNd(int v0)`
	template <typename = std::enable_if_t<!std::is_integral<T>::value>>
	vectorNd(const T& a0) { a.fill(a0); }
	vectorNd(int a0) { a.fill(a0); } // to allow constructing from 0 and 1
	vectorNd(std::initializer_list<T> list) { assign(list.begin()); }

	int size() const { return N; }
	const T& operator [] (int index) const { return a[index]; }
	T& operator [] (int index) { return a[index]; }

	void assign(const T* rhs) {
		for (int i = 0; i < N; i++) {
			a[i] = rhs[i];
		}
	}

	int cmp(const vectorNd& rhs) const {
		for (int i = 0; i < N; i++) {
			if (a[i] < rhs.a[i]) return -1;
			if (a[i] > rhs.a[i]) return +1;
		}
		return 0;
	}

	bool operator == (const vectorNd& rhs) const { return cmp(rhs) == 0; }
	bool operator != (const vectorNd& rhs) const { return cmp(rhs) != 0; }
	bool operator <  (const vectorNd& rhs) const { return cmp(rhs) < 0; }
	bool operator >  (const vectorNd& rhs) const { return cmp(rhs) > 0; }
	bool operator <= (const vectorNd& rhs) const { return cmp(rhs) <= 0; }
	bool operator >= (const vectorNd& rhs) const { return cmp(rhs) >= 0; }

	vectorNd  operator +  (const vectorNd& rhs) const { vectorNd t(*this); t += rhs; return t; }
	vectorNd  operator -  (const vectorNd& rhs) const { vectorNd t(*this); t -= rhs; return t; }
	vectorNd  operator -  ()                    const { vectorNd t; for (int i = 0; i < N; i++) t.a[i] = -a[i]; return t; }
	vectorNd  operator *  (const vectorNd& rhs) const { vectorNd t(*this); t *= rhs; return t; }
	vectorNd  operator /  (const vectorNd& rhs) const { vectorNd t(*this); t /= rhs; return t; }
	vectorNd  operator %  (const vectorNd& rhs) const { vectorNd t(*this); t %= rhs; return t; }

	vectorNd  operator +  (const T& rhs) const { vectorNd t(*this); t += rhs; return t; }
	vectorNd  operator -  (const T& rhs) const { vectorNd t(*this); t -= rhs; return t; }
	vectorNd  operator *  (const T& rhs) const { vectorNd t(*this); t *= rhs; return t; }
	vectorNd  operator /  (const T& rhs) const { vectorNd t(*this); t /= rhs; return t; }
	vectorNd  operator %  (const T& rhs) const { vectorNd t(*this); t %= rhs; return t; }

	vectorNd& operator += (const vectorNd& rhs) { for (int i = 0; i < N; i++) a[i] += rhs.a[i]; return *this; }
	vectorNd& operator -= (const vectorNd& rhs) { for (int i = 0; i < N; i++) a[i] -= rhs.a[i]; return *this; }
	vectorNd& operator *= (const vectorNd& rhs) { for (int i = 0; i < N; i++) a[i] *= rhs.a[i]; return *this; }
	vectorNd& operator /= (const vectorNd& rhs) { for (int i = 0; i < N; i++) a[i] /= rhs.a[i]; return *this; }
	vectorNd& operator %= (const vectorNd& rhs) { for (int i = 0; i < N; i++) a[i] %= rhs.a[i]; return *this; }

	vectorNd& operator += (const T& rhs) { for (int i = 0; i < N; i++) a[i] += rhs; return *this; }
	vectorNd& operator -= (const T& rhs) { for (int i = 0; i < N; i++) a[i] -= rhs; return *this; }
	vectorNd& operator *= (const T& rhs) { for (int i = 0; i < N; i++) a[i] *= rhs; return *this; }
	vectorNd& operator /= (const T& rhs) { for (int i = 0; i < N; i++) a[i] /= rhs; return *this; }
	vectorNd& operator %= (const T& rhs) { for (int i = 0; i < N; i++) a[i] %= rhs; return *this; }

	T abs2() const { T r = 0; for (int i = 0; i < N; i++) r += a[i] * a[i]; return r; }
};

} // math
} // altruct
