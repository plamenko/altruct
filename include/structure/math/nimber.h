#pragma once

#include <algorithm>
#include <type_traits>

#include "algorithm/math/base.h"

namespace altruct {
namespace math {

/**
 * Nimber
 */
template<typename I, int BITS = sizeof(I) * 8>
class nimber {
	typedef typename std::make_unsigned<I>::type U;
public:
	U v;

	nimber(const U& v = 0) : v(v) {}
	nimber(const nimber& rhs) : v(rhs.v) {}

	static U pow2(int n) {
		return (n < BITS) ? U(1) << n : 0;
	}

	nimber high(int k) const {
		return nimber(v >> pow2(k));
	}

	nimber low(int k) const {
		return nimber(v & (pow2(pow2(k)) - 1));
	}

	static nimber join(nimber hi, nimber lo, int k) {
		return nimber((hi.v << pow2(k)) ^ lo.v);
	}
	
	// Multiplication by pow2(pow2(k) - 1);
	// 
	// This is faster than general multiplication because of the simple
	// structure of the multiplicand. Its high_part is the same structure
	// at the next level down, and its low_part is zero.
	//
	// Time T(n) = O(n) + 3 T(n/2) = O(n^(lb 3)).
	static nimber lift(nimber a, int k) {
		//return a * pow2(pow2(k) - 1);
		if (--k < 0) return a;
		nimber ah = a.high(k);
		nimber al = a.low(k);
		if (ah == 0) return join(lift(al, k), 0, k);
		return join(lift(ah + al, k), lift(lift(ah, k), k), k);
	}

	// Square root.
	// 
	// Note that all nimbers have unique square roots.
	//
	// Time T(n) = O(n^(lb 3)).
	nimber sqrt() const {
		int k = level() - 1; if (k < 0) return *this;
		nimber ah = high(k);
		nimber al = low(k);
		nimber as = lift(ah, k) + al;
		return join(ah.sqrt(), as.sqrt(), k);
	}

	// Multiplicative inverse.
	//
	// If we treat nimbers as column vectors of their high and low parts, then
	// the product (ah,al) (bh,bl) can be found as the matrix multiplication:
	//
	//         ( ah + al    ah )  ( bh )
	//         ( lift(ah)   al )  ( bl )
	//
	// Multiplication forms a group, so the matrix is non-singular. The inverse
	// of (ah,al) is found by inverting the matrix, then multiplying by (0,1).
	// When (ah,al) == 0 we return zero for lack of a better value.
	// 
	// Time T(n) = O(n^(lb 3) ln n) + T(n/2) = O(n^(lb 3) ln n).
	nimber inverse() const {
		int k = level() - 1; if (k < 0) return *this;
		nimber ah = high(k);
		nimber al = low(k);
		nimber as = ah + al;
		nimber inv_det = (as * al + lift(ah * ah, k)).inverse();
		return join(ah * inv_det, as * inv_det, k);
	}

	// Karatsuba multiplication.
	//
	// (ah,al)(bh,bl)
	//     = (ah,0)(bh,0) + (ah,0)(0,bl) + (0,al)(bh,0) + (0,al)(0,bl)
	//     = (ah bh, lift(ah bh)) + (ah bl,0) + (al bh,0) + (0,al bl)
	//     = (ah bh + ah bl + al bh, lift(ah bh) + al bl)
	//     = ((ah + al)(bh + bl) + al bl, lift(ah bh) + al bl)
	//
	// Time T(n) = O(n^(lb 3)) + 3 T(n/2) = O(n^(lb 3) lb n).
	// Note, 1.7x slower than schoolbok multiplication for small nimbers!
	static nimber mul2(nimber x, nimber y) {
		if (x < y) std::swap(x, y);
		if (y.v == 0) return nimber(0);
		if (y.v == 1) return x;
		int k = x.level() - 1;
		if (k > y.level() - 1) {
			nimber xh = x.high(k);
			nimber xl = x.low(k);
			return join(mul2(xh, y), mul2(xl, y), k);
		} else { // equal levels
			nimber xh = x.high(k);
			nimber xl = x.low(k);
			nimber yh = y.high(k);
			nimber yl = y.low(k);
			nimber xl_yl = mul2(xl, yl);
			nimber xh_yh = mul2(xh, yh);
			nimber xs_ys = mul2(xh + xl, yh + yl);
			return join(xs_ys + xl_yl, lift(xh_yh, k) + xl_yl, k);
		}
	}

	// Multiplication of powers of 2. 2^m * 2^n.
	static U mul_pow2_impl(int m, int n) {
		U r = pow2(m ^ n);
		for (int i = m & n; i > 0; i &= i - 1) {
			int j = i & -i;
			r = mul(r, pow2(j) ^ pow2(j - 1));
		}
		return r;
	}

	// Memoized multiplication of powers of 2. 2^m * 2^n.
	static U mul_pow2(int m, int n) {
		static U tbl[64][64] = { { 0 } };
		if (m >= 64 || n >= 64) return mul_pow2_impl(m, n);
		if (m < n) std::swap(m, n);
		return tbl[m][n] ? tbl[m][n] : tbl[m][n] = mul_pow2_impl(m, n);
	}

	// Schoolbook multiplication.
	static U mul(U a, U b) {
		if (a <= 1 || b <= 1) return a * b;
		U r = 0;
		for (int i = 0; i < BITS && pow2(i) <= a; i++) {
			if ((a & pow2(i)) == 0) continue;
			for (int j = 0; j < BITS && pow2(j) <= b; j++) {
				if (b & pow2(j)) r ^= mul_pow2(i, j);
			}
		}
		return r;
	}

	bool operator == (const nimber& rhs) const { return (v == rhs.v); }
	bool operator != (const nimber& rhs) const { return (v != rhs.v); }
	bool operator <  (const nimber& rhs) const { return (v <  rhs.v); }
	bool operator >  (const nimber& rhs) const { return (v >  rhs.v); }
	bool operator <= (const nimber& rhs) const { return (v <= rhs.v); }
	bool operator >= (const nimber& rhs) const { return (v >= rhs.v); }

	nimber  operator +  (const nimber& rhs) const { nimber t(*this); t += rhs; return t; }
	nimber  operator -  (const nimber& rhs) const { nimber t(*this); t -= rhs; return t; }
	nimber  operator -  ()                  const { nimber t(*this);           return t; }
	nimber  operator *  (const nimber& rhs) const { nimber t(*this); t *= rhs; return t; }
	nimber  operator /  (const nimber& rhs) const { nimber t(*this); t /= rhs; return t; }
	nimber  operator %  (const nimber& rhs) const { nimber t(*this); t %= rhs; return t; }

	nimber& operator += (const nimber& rhs) { v ^= rhs.v;        return *this; }
	nimber& operator -= (const nimber& rhs) { v ^= rhs.v;        return *this; }
	nimber& operator *= (const nimber& rhs) { v = mul(v, rhs.v); return *this; }
	nimber& operator /= (const nimber& rhs) { return *this *= rhs.inverse(); }
	nimber& operator %= (const nimber& rhs) { v = 0;             return *this; }

	int level() const { int k = 0; for (U w = 2; w && w <= v; w *= w) k++; return k; }
};

template<typename I>
struct identityT<nimber<I>> {
	static nimber<I> of(const nimber<I>& n) {
		return nimber<I>(identityT<I>::of(n.v));
	}
};

template<typename I>
struct zeroT<nimber<I>> {
	static nimber<I> of(const nimber<I>& n) {
		return nimber<I>(zeroT<I>::of(n.v));
	}
};

} // math
} // altruct
