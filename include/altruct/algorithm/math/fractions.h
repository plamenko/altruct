#pragma once

#include "altruct/algorithm/math/base.h"
#include "altruct/structure/math/fraction.h"

#include <vector>
#include <algorithm>

namespace altruct {
namespace math {

/**
 * The next element in the Farey sequence of order n.
 *
 * if f_prev is the left neighbour or -inf, f_next is the right neighbour of f
 * if f_prev is the right neighbour or +inf, f_next is the left neighbour of f
 */
template<typename I>
fraction<I> farey_neighbour(const I& n, const fraction<I>& f_prev, const fraction<I>& f) {
	I p = f_prev.p, q = f_prev.q;
	if (f_prev.q == 0) {
		gcd_ex(f.q, f.p, &p, &q);
		(f_prev.p < 0) ? p = -p : q = -q;
	}
	I k = (n + q) / f.q;
	return fraction<I>(k * f.p - p, k * f.q - q);
}

} // math
} // altruct
