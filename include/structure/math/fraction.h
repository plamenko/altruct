#pragma once

#include "algorithm/math/base.h"

namespace altruct {
namespace math {

template <typename T>
class fraction {
public:
	T p, q;
	fraction(const T &p = 0, const T &q = 1) : p(p), q(q) { reduce(); }
	
	void reduce() {
		T g = gcd(p, q);
		if (g != 1) p /= g, q /= g;
		if (q < 0) p = -p, q = -q;
	}
	
	bool operator == (const fraction &f) const { return (p * f.q == f.p * q); }
	bool operator != (const fraction &f) const { return (p * f.q != f.p * q); }
	bool operator <  (const fraction &f) const { return (p * f.q <  f.p * q); }
	bool operator <= (const fraction &f) const { return (p * f.q <= f.p * q); }
	bool operator >  (const fraction &f) const { return (p * f.q >  f.p * q); }
	bool operator >= (const fraction &f) const { return (p * f.q >= f.p * q); }
	
	fraction& operator += (const fraction &f) { p = p * f.q + f.p * q; q *= f.q; reduce(); return *this; }
	fraction& operator -= (const fraction &f) { p = p * f.q - f.p * q; q *= f.q; reduce(); return *this; }
	fraction& operator *= (const fraction &f) { p = p * f.p;           q *= f.q; reduce(); return *this; }
	fraction& operator /= (const fraction &f) { T f_p = f.p; p *= f.q; q *= f_p; reduce(); return *this; }
	fraction& operator %= (const fraction &f) { p = 0; q = 1;                              return *this; }

	fraction  operator +  (const fraction &f) const { return fraction(p * f.q + f.p * q, q * f.q); }
	fraction  operator -  (const fraction &f) const { return fraction(p * f.q - f.p * q, q * f.q); }
	fraction  operator -  (                 ) const { return fraction(-p, q); }
	fraction  operator *  (const fraction &f) const { return fraction(p * f.p,           q * f.q); }
	fraction  operator /  (const fraction &f) const { return fraction(p * f.q,           q * f.p); }
	fraction  operator %  (const fraction &f) const { return fraction(0); }
};

} // math
} // altruct
