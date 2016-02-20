#pragma once

#include "algorithm/math/base.h"
#include "structure/math/fraction.h"

#include <vector>
#include <algorithm>

namespace altruct {
namespace math {

/**
 * The next element in the Farey sequence of order n.
 *
 * if f_prev is the left neighbour or -inf, f_next is the right neighbour
 * if f_prev is the right neighbour or +inf, f_next is the left neighbour
 */
template<typename T>
fraction<T> farey_neighbour(const T& n, const fraction<T>& f_prev, const fraction<T>& f) {
	T p = f_prev.p, q = f_prev.q;
	if (f_prev.q == 0) {
		gcd_ex(f.q, f.p, &p, &q);
		(f_prev.p < 0) ? p = -p : q = -q;
	}
	T k = (n + q) / f.q;
	return fraction<T>(k * f.p - p, k * f.q - q);
}

} // math
} // altruct
