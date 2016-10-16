#pragma once

#include "algorithm/math/base.h"
#include "algorithm/math/bits.h"

#include <stdint.h>
#include <algorithm>
#include <string>
#include <typeinfo>
#include <type_traits>

namespace altruct {
namespace math {

/** TODO **

mul
	mulw; long multiplication by a single machine word

div
	half machine-word size step; (O(n) m-w divisions, O(n^2) m-w multiplications)
	div_full

to string
	add buffer length to be safe
	use proper format specifier for printf: <inttypes.h>
	use divide-and-conquer to achieve O(n log^2 n) instead O(n^2) time
*/

template<typename T>
class double_int {
public:
	typedef T half_type;
	static const int type_bits = half_type::type_bits * 2;
	static inline int sign(bool is_negative) { return is_negative ? -1 : 0; }

	// DATA
	T hi, lo;

	// CONSTRUCTORS
	double_int(int val = 0) : hi(sign(val < 0)), lo(val) {}
	double_int(const T& val) : hi(sign(val.is_negative())), lo(val) {}
	double_int(const T& rhs_hi, const T& rhs_lo) : hi(rhs_hi), lo(rhs_lo) {}
	double_int(const double_int& rhs) : hi(rhs.hi), lo(rhs.lo) {}
	double_int(const char *buff) : double_int(from_string10(buff)) {}

	// GETTERS
	bool is_negative() const { return hi.is_negative(); }

	// EQUALITY COMPARISON
	bool operator == (const double_int& rhs) const { return (hi == rhs.hi && lo == rhs.lo); }
	bool operator != (const double_int& rhs) const { return (hi != rhs.hi || lo != rhs.lo); }
	
	// SIGNED COMPARISON
	bool operator <  (const double_int& rhs) const { return hi < rhs.hi || (hi == rhs.hi && lo.unsigned_lt(rhs.lo)); }
	bool operator >  (const double_int& rhs) const { return hi > rhs.hi || (hi == rhs.hi && lo.unsigned_gt(rhs.lo)); }
	bool operator <= (const double_int& rhs) const { return hi < rhs.hi || (hi == rhs.hi && lo.unsigned_lte(rhs.lo)); }
	bool operator >= (const double_int& rhs) const { return hi > rhs.hi || (hi == rhs.hi && lo.unsigned_gte(rhs.lo)); }
	
	// UNSIGNED COMPARISON
	bool unsigned_lt(const double_int& rhs) const { return hi.unsigned_lt(rhs.hi) || (hi == rhs.hi && lo.unsigned_lt(rhs.lo)); }
	bool unsigned_gt(const double_int& rhs) const { return hi.unsigned_gt(rhs.hi) || (hi == rhs.hi && lo.unsigned_gt(rhs.lo)); }
	bool unsigned_lte(const double_int& rhs) const { return hi.unsigned_lt(rhs.hi) || (hi == rhs.hi && lo.unsigned_lte(rhs.lo)); }
	bool unsigned_gte(const double_int& rhs) const { return hi.unsigned_gt(rhs.hi) || (hi == rhs.hi && lo.unsigned_gte(rhs.lo)); }
	
	// ADDITION / SUBTRACTION WITH CARRY
	double_int& assign_adc(const double_int& rhs, int& carry) {
		lo.assign_adc(rhs.lo, carry);
		hi.assign_adc(rhs.hi, carry);
		return *this;
	}
	double_int& assign_sbb(const double_int& rhs, int& borrow) {
		lo.assign_sbb(rhs.lo, borrow);
		hi.assign_sbb(rhs.hi, borrow);
		return *this;
	}
	
	// INCREMENT / DECREMENT
	double_int& operator ++ () { int carry = 1; assign_adc(0, carry); return *this; }
	double_int& operator -- () { int borrow = 1; assign_sbb(0, borrow); return *this; }
	double_int operator ++ (int) { double_int r(*this); return ++r; }
	double_int operator -- (int) { double_int r(*this); return --r; }

	// ADDITION / SUBTRACTION
	double_int& operator += (const double_int& rhs) { int carry = 0; return assign_adc(rhs, carry); }
	double_int& operator -= (const double_int& rhs) { int borrow = 0; return assign_sbb(rhs, borrow); }
	double_int operator + (const double_int& rhs) const { double_int r(*this); return r += rhs; }
	double_int operator - (const double_int& rhs) const { double_int r(*this); return r -= rhs; }
	double_int operator + () const { return *this; }
	double_int operator - () const { return double_int(0) -= *this; }
	double_int& negate() { return *this = -*this; }
	
	// MULTIPLICATION
	static double_int<double_int> unsigned_mul_full(const double_int& lhs, const double_int& rhs) {
		int lhs_cy = 0; auto lhs_su = lhs.lo; lhs_su.assign_adc(lhs.hi, lhs_cy);
		int rhs_cy = 0; auto rhs_su = rhs.lo; rhs_su.assign_adc(rhs.hi, rhs_cy);
		auto m0 = T::unsigned_mul_full(lhs.lo, rhs.lo);
		auto m2 = T::unsigned_mul_full(lhs.hi, rhs.hi);
		auto m1 = T::unsigned_mul_full(lhs_su, rhs_su);
		int borrow0 = 0; m1.assign_sbb(m0, borrow0);
		int borrow2 = 0; m1.assign_sbb(m2, borrow2);
		double_int<double_int> r(m2, m0);
		int carry = 0;
		r.lo.hi.assign_adc(m1.lo, carry);
		r.hi.lo.assign_adc(m1.hi, carry);
		r.hi.hi += (lhs_cy & rhs_cy) + carry - borrow0 - borrow2;
		if (rhs_cy) r.hi += double_int(0, lhs_su);
		if (lhs_cy) r.hi += double_int(0, rhs_su);
		return r;
	}

	static double_int unsigned_mul(const double_int& lhs, const double_int& rhs) {
		auto r = T::unsigned_mul_full(lhs.lo, rhs.lo);
		r.hi += T::unsigned_mul(lhs.lo, rhs.hi);
		r.hi += T::unsigned_mul(lhs.hi, rhs.lo);
		return r;
	}

	double_int operator * (const double_int& rhs) const { return unsigned_mul(*this, rhs); }
	double_int& operator *= (const double_int& rhs) { return *this = *this * rhs; }

	// DIVISION
	static double_int unsigned_div(const double_int& a0, const double_int& b0, double_int *r = 0) {
		if (a0.unsigned_lt(b0)) return 0;
		double_int q = 0, a = a0, b = b0;

		// TODO: do half machine-word divisions, not half double-int divisions !!

		//int r_shift = 0;
		//while (a.unsigned_gte(b)) {
		//  // shifting both a and b by the same amount doesn't affect the quotient
		//	int a_lzc = a.leading_zeros_count(); a <<= a_lzc;
		//	int b_lzc = b.leading_zeros_count(); b <<= a_lzc;
		//	r_shift += a_lzc;
		//	if (b_lzc == a_lzc) { a -= b; q += 1; continue; }
		//  // shifting b left so that (b << k).hi >= 2^(T::type_bits/2-1)
		//	int k = max(0, b_lzc - a_lzc - T::type_bits / 2);
		//	double_int c = b << k;
		//	double_int t = T::unsigned_div(a.get_hi(), c.get_hi() + 1);
		//	t <<= k; // correcting the quotient because shifting b left shifted the quotient right
		//	a -= b * t;
		//	q += t;
		//}
		//if (r) *r = a >> r_shift;
		return q;
	}

	static double_int signed_div(const double_int& a0, const double_int& b0, double_int *r = 0) {
		double_int a = a0; if (a0.is_negative()) a.negate();
		double_int b = b0; if (b0.is_negative()) b.negate();
		double_int q = unsigned_div(a, b, r);
		if (r && a0.is_negative()) r->negate();
		return (a0.is_negative() != b0.is_negative()) ? -q : q;
	}
	
	double_int operator / (const double_int& rhs) const { return signed_div(*this, rhs); }
	double_int operator % (const double_int& rhs) const { double_int r; signed_div(*this, rhs, &r); return r; }
	double_int& operator /= (const double_int& rhs) { return *this = *this / rhs; }
	double_int& operator %= (const double_int& rhs) { return *this = *this % rhs; }

	// BITWISE
	double_int& operator &= (const double_int& rhs) { hi &= rhs.hi; lo &= rhs.lo; return *this; }
	double_int& operator |= (const double_int& rhs) { hi |= rhs.hi; lo |= rhs.lo; return *this; }
	double_int& operator ^= (const double_int& rhs) { hi ^= rhs.hi; lo ^= rhs.lo; return *this; }
	double_int operator & (const double_int& rhs) const { double_int r(*this); return r &= rhs; }
	double_int operator | (const double_int& rhs) const { double_int r(*this); return r |= rhs; }
	double_int operator ^ (const double_int& rhs) const { double_int r(*this); return r ^= rhs; }
	double_int operator ~ () const { return double_int(~hi, ~lo); }

	// SHIFTS
	double_int& operator <<= (int cnt) {
		if (cnt >= 2 * T::type_bits) {
			hi = 0;
			lo = 0;
		} else if (cnt > T::type_bits) {
			hi = lo;
			hi <<= cnt - T::type_bits;
			lo = 0;
		} else if (cnt == T::type_bits) {
			hi = lo;
			lo = 0;
		} else if (cnt > 0) {
			hi <<= cnt;
			hi |= double_int(lo).assign_unsigned_shr(T::type_bits - cnt);
			lo <<= cnt;
		}
		return *this;
	} 
	
	double_int& operator >>= (int cnt) {
		return assign_extended_shr(cnt, sign(is_negative()));
	}

	double_int& assign_unsigned_shr(int cnt) {
		return assign_extended_shr(cnt, 0);
	}
	
	double_int& assign_extended_shr(int cnt, int ext) {
		if (cnt >= 2 * T::type_bits) {
			lo = ext;
			hi = ext;
		} else if (cnt > T::type_bits) {
			lo = hi;
			lo.assign_extended_shr(cnt - T::type_bits, ext);
			hi = ext;
		} else if (cnt == T::type_bits) {
			lo = hi;
			hi = ext;
		} else if (cnt > 0) {
			lo.assign_unsigned_shr(cnt);
			lo |= hi << (T::type_bits - cnt);
			hi.assign_extended_shr(cnt, ext);
		}
		return *this;
	}
	
	double_int operator << (int cnt) const { double_int r(*this); return r <<= cnt; } 
	double_int operator >> (int cnt) const { double_int r(*this); return r >>= cnt; } 

	int leading_zeros_count() const {
		int hi_lzc = hi.leading_zeros_count();
		return (hi_lzc < T::type_bits) ? hi_lzc : T::type_bits + lo.leading_zeros_count();
	}

	// INPUT / OUTPUT

	uint64_t to_uint64() const {
		if (T::type_bits >= 64) return lo.to_uint64();
		return (hi.to_uint64() << T::type_bits) | lo.to_uint64();
	}
	
	//static double_int pow10(int e) {
	//	double_int r = 1;
	//	while (e-- > 0) {
	//		r += r;  if (r.is_negative()) return r;
	//		double_int r2 = r;
	//		r += r;  if (r.is_negative()) return r;
	//		r += r;  if (r.is_negative()) return r;
	//		r += r2; if (r.is_negative()) return r;
	//	}
	//	return r;
	//}

	//// WARNING: does not check for buffer overflow!
	//// LIMITATION: outputs "-OF" for minimum value
	//char *to_string10(char *buff) const {
	//	char *wr = buff;
	//	double_int t = *this;
	//	if (t == 0) {
	//		*(wr++) = '0';
	//		*(wr++) = 0;
	//		return buff;
	//	}
	//	if (t.is_negative()) {
	//		*(wr++) = '-';
	//		t = -t;
	//	}
	//	if (t.is_negative()) {
	//		*(wr++) = 'O';
	//		*(wr++) = 'F';
	//		*(wr++) = 0;
	//		return buff;
	//	}
	//	int digits = 0;
	//	double_int p10;
	//	do {
	//		digits++;
	//		p10 = pow10(digits);
	//	} while (!p10.is_negative() && p10 <= t);
	//	while (digits-- > 0) {
	//		p10 = pow10(digits);
	//		int d = 0;
	//		while (t >= p10) {
	//			t -= p10;
	//			d++;
	//		}
	//		*(wr++) = '0' + d;
	//	}
	//	*(wr++) = 0;
	//	return buff;
	//}

	//static double_int from_string10(const char *buff) {
	//	double_int r = 0;
	//	bool is_neg = false;
	//	if (*buff == '-') {
	//		buff++;
	//		is_neg = true;
	//	}
	//	while (*buff) {
	//		int d = *(buff++) -'0';
	//		r += r;
	//		double_int r2 = r;
	//		r += r;
	//		r += r;
	//		r += r2;
	//		r += d;
	//	}
	//	return is_neg ? -r : r;
	//}

	std::string to_string16() const {
		return hi.to_string16() + lo.to_string16();
	}
};


/**
 * A wrapper around a signed primitive integral type
 * that adds some unsigned facilities.
 *
 * Accepts:
 *   int8_t, uint8_t, uint16_t
 *   int16_t, uint16_t, uint32_t
 *   int32_t, uint32_t, uint64_t
 *   int64_t, uint64_t, <uint128_t>
 *
 * @param sT - signed variant of the wrapped type
 * @param uT - unsigned variant of the wrapped type
 * @param INTR - intrinsics for `adc`, `sbb` and `umul`.
 */
template<typename sT, typename uT, typename INTR>
class prim_int {
public:
	static const int type_bits = sizeof(sT) * 8;

	// DATA
	sT v;

	// CONSTRUCTORS
	prim_int() : v(0) {}
	prim_int(const sT& val) : v(val) {}
	prim_int(const uT& val) : v(sT(val)) {}
	// construct from int, but only if sT is not int itself to avoid constructor clashing
	template <typename sI = sT, typename = std::enable_if_t<!std::is_same<sI, int>::value>>
	prim_int(int val) : v(sT(val)) {}

	// GETTERS
	bool is_negative() const { return v < 0; }

	// FORWARDING
	bool operator == (const prim_int& rhs) const { return (v == rhs.v); }
	bool operator != (const prim_int& rhs) const { return (v != rhs.v); }
	
	bool operator <  (const prim_int& rhs) const { return (v <  rhs.v); }
	bool operator >  (const prim_int& rhs) const { return (v >  rhs.v); }
	bool operator <= (const prim_int& rhs) const { return (v <= rhs.v); }
	bool operator >= (const prim_int& rhs) const { return (v >= rhs.v); }

	prim_int operator ++ (int) { return prim_int(v++); } // prim_int++
	prim_int operator -- (int) { return prim_int(v--); } // prim_int--

	prim_int operator - () const { return prim_int(-v); }
	prim_int operator + () const { return prim_int(+v); }

	prim_int operator + (const prim_int& rhs) const { return prim_int(v + rhs.v); }
	prim_int operator - (const prim_int& rhs) const { return prim_int(v - rhs.v); }
	prim_int operator * (const prim_int& rhs) const { return prim_int(v * rhs.v); }
	prim_int operator / (const prim_int& rhs) const { return prim_int(v / rhs.v); }
	prim_int operator % (const prim_int& rhs) const { return prim_int(v % rhs.v); }

	prim_int operator ~ () const { return prim_int(~v); }

	prim_int operator & (const prim_int& rhs) const { return prim_int(v & rhs.v); }
	prim_int operator | (const prim_int& rhs) const { return prim_int(v | rhs.v); }
	prim_int operator ^ (const prim_int& rhs) const { return prim_int(v ^ rhs.v); }

	prim_int operator << (int cnt) const { return prim_int(v << cnt); }
	prim_int operator >> (int cnt) const { return prim_int(v >> cnt); }

	prim_int& operator ++ () { ++v; return *this; } // ++prim_int
	prim_int& operator -- () { --v; return *this; } // --prim_int

	prim_int& operator += (const prim_int& rhs) { v += rhs.v; return *this; }
	prim_int& operator -= (const prim_int& rhs) { v -= rhs.v; return *this; }
	prim_int& operator *= (const prim_int& rhs) { v *= rhs.v; return *this; }
	prim_int& operator /= (const prim_int& rhs) { v /= rhs.v; return *this; }
	prim_int& operator %= (const prim_int& rhs) { v %= rhs.v; return *this; }

	prim_int& negate() { return prim_int(-v); }

	prim_int& operator &= (const prim_int& rhs) { v &= rhs.v; return *this; }
	prim_int& operator |= (const prim_int& rhs) { v |= rhs.v; return *this; }
	prim_int& operator ^= (const prim_int& rhs) { v ^= rhs.v; return *this; }

	prim_int& operator <<= (int cnt) { v <<= cnt; return *this; }
	prim_int& operator >>= (int cnt) { v >>= cnt; return *this; }

	// UNSIGNED COMPARISON
	bool unsigned_lt(const prim_int& rhs) const { return (uT(v) < uT(rhs.v)); }
	bool unsigned_gt(const prim_int& rhs) const { return (uT(v) > uT(rhs.v)); }
	bool unsigned_lte(const prim_int& rhs) const { return (uT(v) <= uT(rhs.v)); }
	bool unsigned_gte(const prim_int& rhs) const { return (uT(v) >= uT(rhs.v)); }

	// UNSIGNED SHIFT RIGHT
	prim_int& assign_unsigned_shr(int cnt) { v = uT(v) >> cnt; return *this; }
	prim_int& assign_extended_shr(int cnt, int ext) { v = (uT(v) >> cnt) | (uT(ext) << (type_bits - cnt)); return *this; } // undefined for cnt = 0
	int leading_zeros_count() const { return lzc(uT(v)); }

	// ADDITION / SUBTRACTION WITH CARRY
	prim_int& assign_adc(const prim_int& rhs, int& carry) {
		v = INTR().adc(uT(v), uT(rhs.v), carry);
		return *this;
	}
	prim_int& assign_sbb(const prim_int& rhs, int& borrow) {
		v = INTR().sbb(uT(v), uT(rhs.v), borrow);
		return *this;
	}

	// UNSIGNED MULTIPLICATION
	static double_int<prim_int> unsigned_mul_full(const prim_int& lhs, const prim_int& rhs) {
		double_int<prim_int> r;
		r.lo = INTR().umul(uT(lhs.v), uT(rhs.v), (uT*)&(r.hi.v));
		return r;
	}
	static prim_int unsigned_mul(const prim_int& lhs, const prim_int& rhs) {
		return uT(lhs.v) * uT(rhs.v);
	}

	// UNSIGNED DIVISION
	static prim_int unsigned_div_full(const prim_int& lhs_hi, const prim_int& lhs_lo, const prim_int& rhs, prim_int *r = 0) {
		// TODO
	}
	static prim_int unsigned_div(const prim_int& lhs, const prim_int& rhs, prim_int *r = 0) {
		if (r) *r = uT(lhs.v) % uT(rhs.v);
		return uT(lhs.v) / uT(rhs.v);
	}
	
	// OUTPUT
	uint64_t to_uint64() const { return uT(v); }
	
	std::string to_string16() const {
		static const char hex_digits[16] = {
			'0', '1', '2', '3', '4', '5', '6', '7',
			'8', '9', 'a', 'b', 'c', 'd', 'e', 'f' };
		int num_digits = type_bits / 4;
		std::string s(num_digits, '0');
		auto u = to_uint64();
		for (int i = num_digits - 1; i >= 0; i--) {
			s[i] = hex_digits[u & 0xF];
			u >>= 4;
		}
		return s;
	}
};

template<typename T>
struct identityT<double_int<T>> {
	static double_int<T> of(const double_int<T>& x) {
		return double_int<T>(1);
	}
};

template<typename T>
struct zeroT<double_int<T>> {
	static double_int<T> of(const double_int<T>& x) {
		return double_int<T>(0);
	}
};

} // math
} // altruct
