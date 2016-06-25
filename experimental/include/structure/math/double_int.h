#pragma once

#include "algorithm/math/base.h"

#include <stdint.h>
#include <algorithm>
#include <type_traits>

namespace altruct {
namespace math {

/** TODO **

add/sub ?
	add_lo, add_hi
	add_l0 -> add
	add_hi -> adc (add with carry)
	no other instructions may be done there in order to preserve the carry flag

mul
	mul_full(a, b) -> (hi, lo)
	mull_mw; long multiplication by a single machine word

div
	half machine-word size step; (O(n) m-w divisions, O(n^2) m-w multiplications)

lzc ?
	built-in instruction for primitive_int
	
*/

template<typename T>
class double_int {
public:
	// STATIC
	static inline int type_bits() { return T::type_bits() * 2; }
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
	
	//INCREMENT
	
	double_int& operator ++ () { // ++double_int
		++lo;
		if (lo.unsigned_lt(1)) ++hi; // carry
		return *this;
	}

	double_int operator ++ (int) { // double_int++
		double_int t = *this; ++*this; return t;
	}
	
	// ADDITION

	double_int& operator += (const double_int& rhs) {
		T t = lo; hi += rhs.hi; lo += rhs.lo;
		if (lo.unsigned_lt(t)) ++hi; // carry
		return *this;
	}

	double_int& operator += (const T& rhs) {
		T t = lo; lo += rhs;
		if (lo.unsigned_lt(t)) ++hi; // carry
		return *this;
	}

	double_int& operator += (int rhs) {
		T t = lo; lo += rhs;
		if (lo.unsigned_lt(t)) ++hi; // carry
		return *this;
	}

	double_int operator + (const double_int& rhs) const { double_int r(*this); return r += rhs; }

	double_int operator + () const { return *this; }

	// DECREMENT
	
	double_int& operator -- () { // --double_int
		T t = lo; --lo;
		if (lo.unsigned_gt(t)) --hi; // borrow
		return *this;
	}
	
	double_int operator -- (int) { // double_int--
		double_int t = *this; --*this; return t;
	}

	// SUBTRACTION
	
	double_int& operator -= (const double_int& rhs) {
		T t = lo; hi -= rhs.hi; lo -= rhs.lo;
		if (lo.unsigned_gt(t)) --hi; // borrow
		return *this;
	}
	
	double_int operator - (const double_int& rhs) const { double_int r(*this); return r -= rhs; }

	double_int operator - () const { return double_int(0) - *this; }
	
	// MULTIPLICATION
	
	static std::pair<double_int, double_int> unsigned_mul_full(const double_int& lhs, const double_int& rhs) const {
		// 2^(3*L/2) (hm2) + 2^(2*L/2) (lm2 +  hm2 + hm0 - hm1) + 2 ^ (1 * L / 2) (hm0 + lm0 + lm2 - lm1) + 2 ^ (0 * L / 2) (lm0)
		auto m2 = T::unsigned_mul_full(lhs.hi, rhs.hi);
		auto m1 = T::unsigned_mul_full(lhs.hi - lhs.lo, rhs.hi - rhs.lo); // TODO: addition / subtraction can overflow
		auto m0 = T::unsigned_mul_full(lhs.lo, rhs.lo);
		//int64_t x = to_int64();
		//int64_t y = rhs.to_int64();
		//double_int r;
		//T a(0, lo.get_hi()), b(0, lo.get_lo());
		//T c(0, rhs.lo.get_hi()), d(0, rhs.lo.get_lo());
		//r.lo = b.unsigned_mul(d);
		//T ad = a.unsigned_mul(d);
		//r += T(0, ad.get_lo()) << (T::type_bits() / 2);
		//r.hi += T(0, ad.get_hi());
		//T bc = b.unsigned_mul(c);
		//r += T(0, bc.get_lo()) << (T::type_bits() / 2);
		//r.hi += T(0, bc.get_hi());
		//r.hi += a.unsigned_mul(c);
		//r.hi += hi.unsigned_mul(rhs.lo);
		//r.hi += lo.unsigned_mul(rhs.hi);
		return r;
	}

	static std::pair<double_int, double_int> unsigned_mul(const double_int& lhs, const double_int& rhs) const {
	
	}
	
	double_int operator * (const double_int& rhs) const {
		return unsigned_mul(rhs);
	}

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
		//  // shifting b left so that (b << k).hi >= 2^(T::type_bits()/2-1)
		//	int k = max(0, b_lzc - a_lzc - T::type_bits() / 2);
		//	double_int c = b << k;
		//	double_int t = T::unsigned_div(a.get_hi(), c.get_hi() + 1);
		//	t <<= k; // correcting the quotient because shifting b left shifted the quotient right
		//	a -= b * t;
		//	q += t;
		//}
		//if (r) *r = a >> r_shift;
		return q;
	}
	
	double_int operator / (const double_int& rhs) const {
		int s = 1;
		double_int a = *this; if (a.is_negative()) a = -a, s =-s;
		double_int b = rhs;   if (b.is_negative()) b = -b, s =-s;
		double_int q = unsigned_div(a, b);
		return (s < 0) ? -q : q;
	}

	double_int& operator /= (const double_int& rhs) { return *this = *this / rhs; }

	double_int operator % (const double_int& rhs) const {
		int s = 1;
		double_int r;
		double_int a = *this; if (a.is_negative()) a = -a, s =-s;
		double_int b = rhs;   if (b.is_negative()) b = -b;
		unsigned_div(a, b, &r);
		return (s < 0) ? -r : r;
	}

	double_int& operator %= (const double_int& rhs) { return *this = *this % rhs; }

	// BITWISE
	
	double_int operator ~ () const { return double_int(~hi, ~lo); }
	
	double_int& operator &= (const double_int& rhs) { hi &= rhs.hi; lo &= rhs.lo; return *this; }
	double_int& operator |= (const double_int& rhs) { hi |= rhs.hi; lo |= rhs.lo; return *this; }
	double_int& operator ^= (const double_int& rhs) { hi ^= rhs.hi; lo ^= rhs.lo; return *this; }

	double_int operator & (const double_int& rhs) const { double_int r(*this); return r &= rhs; }
	double_int operator | (const double_int& rhs) const { double_int r(*this); return r |= rhs; }
	double_int operator ^ (const double_int& rhs) const { double_int r(*this); return r ^= rhs; }

	// SHIFTS
	
	double_int& operator <<= (int cnt) {
		if (cnt >= 2 * T::type_bits()) {
			hi = 0;
			lo = 0;
		} else if (cnt > T::type_bits()) {
			hi = lo;
			hi <<= cnt - T::type_bits();
			lo = 0;
		} else if (cnt == T::type_bits()) {
			hi = lo;
			lo = 0;
		} else if (cnt > 0) {
			hi <<= cnt;
			hi |= double_int(lo).assign_unsigned_shr(T::type_bits() - cnt);
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
		if (cnt >= 2 * T::type_bits()) {
			lo = ext;
			hi = ext;
		} else if (cnt > T::type_bits()) {
			lo = hi;
			lo.assign_extended_shr(cnt - T::type_bits(), ext);
			hi = ext;
		} else if (cnt == T::type_bits()) {
			lo = hi;
			hi = ext;
		} else if (cnt > 0) {
			lo.assign_unsigned_shr(cnt);
			lo |= hi << (T::type_bits() - cnt);
			hi.assign_extended_shr(cnt, ext);
		}
		return *this;
	}
	
	double_int operator << (int cnt) const { double_int r(*this); return r <<= cnt; } 
	double_int operator >> (int cnt) const { double_int r(*this); return r >>= cnt; } 

	int leading_zeros_count() const {
		int hi_lzc = hi.leading_zeros_count();
		return (hi_lzc < T::type_bits()) ? hi_lzc : T::type_bits() + lo.leading_zeros_count();
	}

	// INPUT / OUTPUT

	int64_t to_int64() const {
		if (T::type_bits() >= 64) return lo.to_int64();
		int64_t lo_mask = (1LL << T::type_bits()) - 1;
		return (hi.to_int64() << T::type_bits()) | (lo.to_int64() & lo_mask);
	}
	
	static double_int pow10(int e) {
		double_int r = 1;
		while (e-- > 0) {
			r += r;  if (r.is_negative()) return r;
			double_int r2 = r;
			r += r;  if (r.is_negative()) return r;
			r += r;  if (r.is_negative()) return r;
			r += r2; if (r.is_negative()) return r;
		}
		return r;
	}

	// WARNING: does not check for buffer overflow!
	// LIMITATION: outputs "-OF" for minimum value
	char *to_string10(char *buff) const {
		char *wr = buff;
		double_int t = *this;
		if (t == 0) {
			*(wr++) = '0';
			*(wr++) = 0;
			return buff;
		}
		if (t.is_negative()) {
			*(wr++) = '-';
			t = -t;
		}
		if (t.is_negative()) {
			*(wr++) = 'O';
			*(wr++) = 'F';
			*(wr++) = 0;
			return buff;
		}
		int digits = 0;
		double_int p10;
		do {
			digits++;
			p10 = pow10(digits);
		} while (!p10.is_negative() && p10 <= t);
		while (digits-- > 0) {
			p10 = pow10(digits);
			int d = 0;
			while (t >= p10) {
				t -= p10;
				d++;
			}
			*(wr++) = '0' + d;
		}
		*(wr++) = 0;
		return buff;
	}

	static double_int from_string10(const char *buff) {
		double_int r = 0;
		bool is_neg = false;
		if (*buff == '-') {
			buff++;
			is_neg = true;
		}
		while (*buff) {
			int d = *(buff++) -'0';
			r += r;
			double_int r2 = r;
			r += r;
			r += r;
			r += r2;
			r += d;
		}
		return is_neg ? -r : r;
	}
};


// Wrapper around signed primitive integral type
// that adds some unsigned facilities.
// Accepts:
//   int8_t, uint8_t
//   int16_t, uint16_t
//   int32_t, uint32_t
//   int64_t, uint64_t
template<typename sT, typename uT>
class prim_int {
public:
	// STATIC
	static inline int type_bits() { return sizeof(sT) * 8; }

	// DATA
	sT v;

	// CONSTRUCTORS
	prim_int(const sT& val = 0) : v(val) {}
	prim_int(const sT& hi, const sT& lo) : v(lo) {}

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

	prim_int& operator &= (const prim_int& rhs) { v &= rhs.v; return *this; }
	prim_int& operator |= (const prim_int& rhs) { v |= rhs.v; return *this; }
	prim_int& operator ^= (const prim_int& rhs) { v ^= rhs.v; return *this; }

	prim_int& operator <<= (int cnt) { v <<= cnt; return *this; }
	prim_int& operator >>= (int cnt) { v >>= cnt; return *this; }

	// UNSIGNED MULTIPLICATION
	static void _add(uT& hi, uT& lo, uT hi2, uT lo2) {
		lo += lo2; if (lo < lo2) hi++; hi += hi2;
	}
	static std::pair<prim_int, prim_int> unsigned_mul_full(const prim_int& lhs, const prim_int& rhs) const {
		const int L2 = type_bits() / 2;
		uT hi1 = lhs.get_hi(), lo1 = lhs.get_lo(), hi2 = rhs.get_hi(), lo2 = rhs.get_lo();
		uT hi(hi1 * hi2), lo(lo1 * lo2), m1(hi1 * lo2), m2(lo1 * hi2);
		_add(hi, lo, m1 >> L2, m1 << L2); _add(hi, lo, m2 >> L2, m2 << L2);
		return{ hi, lo };
	}
	static prim_int unsigned_mul(const prim_int& lhs, const prim_int& rhs) const {
		return (uT)lhs.v * (uT)rhs.v;
	}

	// UNSIGNED DIVISION
	static prim_int unsigned_div(const prim_int& lhs, const prim_int& rhs, prim_int *r = 0) {
		if (r) *r = (uT) lhs.v % (uT) rhs.v;
		return (uT) lhs.v / (uT) rhs.v;
	}
	
	// UNSIGNED COMPARISON
	bool unsigned_lt(const prim_int& rhs) const { return ((uT)v < (uT)rhs.v); }
	bool unsigned_gt(const prim_int& rhs) const { return ((uT)v >(uT)rhs.v); }
	bool unsigned_lte(const prim_int& rhs) const { return ((uT)v <= (uT)rhs.v); }
	bool unsigned_gte(const prim_int& rhs) const { return ((uT)v >= (uT)rhs.v); }

	// UNSIGNED SHIFT RIGHT
	static uT _ushr(uT v, int cnt) { return v >> cnt; }
	prim_int& assign_unsigned_shr(int cnt) { v = _ushr(v, cnt); return *this; }
	prim_int& assign_extended_shr(int cnt, int ext) { v = _ushr(v, cnt) | (ext << (type_bits() - cnt)); return *this; }

	int leading_zeros_count() const {
		uT u = (uT) v;
		switch (type_bits()) {
			case 64:
				if (u >> 3*16) return 1*16 - len_lookup(int(u >> 3*16));
				if (u >> 2*16) return 2*16 - len_lookup(int(u >> 2*16));
				if (u >> 1*16) return 3*16 - len_lookup(int(u >> 1*16));
				               return 4*16 - len_lookup(int(u));
			case 32:
				if (u >> 16) return 16 - len_lookup(int(u >> 16));
				             return 32 - len_lookup(int(u));
			case 16:
				return 16 - len_lookup(int(u));
			case 8:
				return 8 - len_lookup(int(u));
			default:
				return 0;
		}
	}

	static int len_lookup(int index) {
		static bool len_lookup_initialized = false;
		static char len_lookup_table[1 << 16] = {};
		if (!len_lookup_initialized) {
			for (int i = 0; i < (1 << 16); i++) {
				int cnt = 0;
				for (int j = i; j > 0; j >>= 1) {
					cnt++;
				}
				len_lookup_table[i] = cnt;
			}
			len_lookup_initialized = true;
		}
		return len_lookup_table[index];
	}
	
	// OUTPUT
	int64_t to_int64() const { return int64_t(v); }
	char *to_string10(char *buff) const { sprintf(buff, "%lld", int64_t(v)); return buff;}
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
		return double_int<T>(1);
	}
};

} // math
} // altruct
