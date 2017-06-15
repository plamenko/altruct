#pragma once

#include "altruct/algorithm/math/base.h"
#include "altruct/algorithm/math/sums.h"

namespace altruct {
namespace math {
namespace sequences {

/**
 * Some common sequences.
 */

template<typename R, typename I>
R castR(I n) { return castOf<R>(n); }

// Delta function: delta(n) = [n == 0]
template<typename R, typename I>
static R delta(I n) { auto r = castR<R>(n); return (n == 0) ? identityOf(r) : zeroOf(r); }

// Dirichlet multiplicative identity: dirichlet_id(n) = [n == 1]
template<typename R, typename I>
static R dirichlet_id(I n) { auto r = castR<R>(n); return (n == 1) ? identityOf(r) : zeroOf(r); }

// Constant function: one(n) = 1
template<typename R, typename I>
static R zero(I n) { auto r = castR<R>(n); return zeroOf(r); }

// Constant function: one(n) = 1
template<typename R, typename I>
static R one(I n) { auto r = castR<R>(n); return identityOf(r); }

// Identity function: identity(n) = n
template<typename R, typename I>
static R identity(I n) { auto r = castR<R>(n); return r; }

// Suare numbers: square(n) = n^2
template<typename R, typename I>
static R square(I n) { auto r = castR<R>(n); return r * r; }

// Cube numbers: cube(n) = n^3
template<typename R, typename I>
static R cube(I n) { auto r = castR<R>(n); return r * r * r; }

// Triangular numbers: triangular(n) = Sum[k, {k,1,n}] = n(n+1)/2
template<typename R, typename I>
static R triangular(I n) { auto r = castR<R>(n); return r * (r + 1) / 2; }

// Tetrahedral (Triangular pyramidal) numbers: tetrahedral(n) = n*(n+1)*(n+2)/6
template<typename R, typename I>
static R tetrahedral(I n) { auto r = castR<R>(n); return r * (r + 1) * (r + 2) / 6; }

// Square pyramidal numbers: pyramidal(n) = Sum[k^2, {k,1,n}]
template<typename R, typename I>
static R pyramidal(I n) { auto r = castR<R>(n); return r * (r + 1) * (r * 2 + 1) / 6; }

// Octahedral numbers: octahedral(n) = n*(2*n^2+1)/3
template<typename R, typename I>
static R octahedral(I n) { auto r = castR<R>(n); return r * (r * r * 2 + 1) / 3; }

// Dodecahedral numbers: dodecahedral(n) = n*(3*n-1)*(3*n-2)/2
template<typename R, typename I>
static R dodecahedral(I n) { auto r = castR<R>(n); return r * (r * 3 - 1) * (r * 3 - 2) / 2; }

// Icosahedral numbers: icosahedral(n) = n*(5*n^2-5*n+2)/2
template<typename R, typename I>
static R icosahedral(I n) { auto r = castR<R>(n); return r * (r * (r - 1) * 5 + 2) / 2; }

// Partial sums of DivisorSigma0 function. O(sqrt n)
template<typename R, typename I>
static R sum_sigma0(I n) { auto r = castR<R>(n); return sum_sqrt2m(one<R, I>, identity<R, I>, identity<R, I>, n, zeroOf(r)); }

// Partial sums of DivisorSigma1 function. O(sqrt n)
template<typename R, typename I>
static R sum_sigma1(I n) { auto r = castR<R>(n); return sum_sqrt2m(identity<R, I>, triangular<R, I>, identity<R, I>, n, zeroOf(r)); }

// Partial sums of DivisorSigma2 function. O(sqrt n)
template<typename R, typename I>
static R sum_sigma2(I n) { auto r = castR<R>(n); return sum_sqrt2m(square<R, I>, pyramidal<R, I>, identity<R, I>, n, zeroOf(r)); }

} // sequences
} // math
} // altruct
