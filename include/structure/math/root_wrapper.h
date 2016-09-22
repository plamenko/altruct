#pragma once

#include "algorithm/math/base.h"
#include "structure/math/complex.h"

#include <cmath>
#include <vector>

namespace altruct {
namespace math {

/**
 * A root wrapper that wraps array of all `k`
 * k-th roots of unity and can be used instead of T.
 * Note: the size of the array must be a power of two.
 *
 * `roots[i] = w^i`, where `w` is the principal k-th root.
 *
 * This is useful to avoid precision issues as all roots
 * can be computed upfront with sufficient precission,
 * instead of deriving them from the principal k-th root
 * (by exponentiation) which would accumulate the error.
 */
template<typename T>
class root_wrapper {
public:
	int index, size;
	const T* roots;
	root_wrapper(const T* roots, int size, int index) : roots(roots), size(size), index(index) {}
	root_wrapper& operator *= (const root_wrapper& rhs) { index += rhs.index; index &= size - 1; return *this; }
	root_wrapper operator * (const root_wrapper& rhs) const { return root_wrapper(*this) *= rhs; }
	operator T() const { return roots[index]; }
};
template<typename T>
struct identityT<root_wrapper<T>> {
	static root_wrapper<T> of(const root_wrapper<T>& w) {
		return root_wrapper<T>(w.roots, w.size, 0);
	}
};

/**
 * Returns a root_wrapper<cplx> of principal k-th root
 * of unity for some power of 2 `k` no smaller than `l`.
 */
template<typename F>
root_wrapper<complex<F>> complex_root_wrapper(int l) {
	typedef complex<F> cplx;
	static const auto _2_PI = 2 * acos(F(-1));
	static std::vector<cplx> roots{F(1)};
	int size = (int)roots.size();
	if (size < l) {
		while (size < l) size *= 2;
		roots.resize(size);
		for (int i = 0; i < size; i++) {
			auto A = _2_PI * i / size;
			roots[i] = cplx(cos(A), sin(A));
		}
	}
	return root_wrapper<cplx>(roots.data(), size, 1);
}

} // math
} // altruct
