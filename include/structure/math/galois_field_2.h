#pragma once

#include "algorithm/math/base.h"

namespace altruct {
namespace math {

/**
 * Galois field of two elements: GF(2), F2, Z/2Z, Z2
 *
 * The two elements are 0 and 1, being the additive and multiplicative identities respectively.
 * The field's addition operation corresponds to the logical XOR operation.
 * The field's multiplication operation corresponds to the logical AND operation.
 * Note: to allow for SWAR, operations are performed on a word instead of a single bit.
 *
 * @param I - the underlying type
 */
template<typename I = uint32_t>
class galois_field_2 {
public:
	I v;

	galois_field_2(const I& v = 0) : v(v) {}

	bool operator == (const galois_field_2 &rhs) const { return (v == rhs.v); }
	bool operator != (const galois_field_2 &rhs) const { return (v != rhs.v); }
	bool operator <  (const galois_field_2 &rhs) const { return (v <  rhs.v); }
	bool operator >  (const galois_field_2 &rhs) const { return (v >  rhs.v); }
	bool operator <= (const galois_field_2 &rhs) const { return (v <= rhs.v); }
	bool operator >= (const galois_field_2 &rhs) const { return (v >= rhs.v); }

	galois_field_2  operator +  (const galois_field_2 &rhs) const { galois_field_2 r(*this); return r += rhs; }
	galois_field_2  operator -  (const galois_field_2 &rhs) const { galois_field_2 r(*this); return r -= rhs; }
	galois_field_2  operator -  ()                          const { galois_field_2 r(*this); return r; }
	galois_field_2  operator *  (const galois_field_2 &rhs) const { galois_field_2 r(*this); return r *= rhs; }
	galois_field_2  operator /  (const galois_field_2 &rhs) const { galois_field_2 r(*this); return r /= rhs; }

	galois_field_2& operator += (const galois_field_2 &rhs) { v ^= rhs.v; return *this; }
	galois_field_2& operator -= (const galois_field_2 &rhs) { v ^= rhs.v; return *this; }
	galois_field_2& operator *= (const galois_field_2 &rhs) { v &= rhs.v; return *this; }
	galois_field_2& operator /= (const galois_field_2 &rhs) { v |= ~rhs.v; return *this; }
};

template<typename I>
struct identityT<galois_field_2<I>> {
	static galois_field_2<I> of(const galois_field_2<I>& x) {
		return galois_field_2<I>(~I(0));
	}
};

template<typename I>
struct zeroT<galois_field_2<I>> {
	static galois_field_2<I> of(const galois_field_2<I>& x) {
		return galois_field_2<I>(I(0));
	}
};

} // math
} // altruct
