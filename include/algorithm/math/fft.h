#pragma once

#include "base.h"
#include <iterator>
#include <vector>

namespace altruct {
namespace math {

/**
 * Inplace Fast Fourier Transform of a sequence
 *
 * @param data - data to transform, array of length `size`
 * @param size - number of elements, must be a power of two
 * @param root - a principal n-th root of unity in the ring T
 */
template<typename T, typename R>
void fft(T *data, int size, R root) {
	T w, e1 = identityT<R>::of(root);
	int m, h, i, j, k;
	for (m = size; h = m / 2, m > 1; m /= 2, root *= root) {
		for (i = 0, w = e1; i < h; i++, w *= root) {
			for (j = i; j < size; j += m) {
				k = j + h;
				T t = data[j] - data[k];
				data[j] += data[k];
				data[k] = t * w;
			}
		}
	}
	// reorder: swap(a[j], a[bitrev(j)])
	for (int i = 0, j = 1; j < size - 1; j++) {
		for (int k = size / 2; (i ^= k) < k; k /= 2);
		if (j < i) std::swap(data[i], data[j]);
	}
}

/**
 * Fast Fourier Transform of a sequence
 *
 * This implementation has better numerical stability than `fft`.
 *
 * @param dest - destination array of length `size` for result
 * @param src  - source data to transform, array of length `size`
 * @param size - number of elements, must be a power of two
 * @param root - a principal n-th root of unity in the ring T
 */
template<typename T, typename R>
void fft_rec(T *dest, T *src, int size, const R& root, int off = 1) {
	if (size <= 1) { *dest = *src; return; }
	int h = size / 2;
	R root2 = root * root;
	R rooti = identityT<R>::of(root);
	fft_rec(dest, src, h, root2, off * 2);
	fft_rec(dest + h, src + off, h, root2, off * 2);
	for (int i = 0; i < h; i++, rooti *= root) {
		T z = dest[i + h] * rooti;
		dest[i + h] = dest[i] - z;
		dest[i] += z;
	}
}

/**
 * FFT Cyclic Convolution of two sequences
 *
 * Result is stored in `dataR`. All `dataR`, `data1` and `data2` are modified.
 * dataR[k] = Sum[data1[i] * data2[(k - i) % size], {i, 0, size - 1}]
 *
 * Note: in case `z1 + z2 >= size - 1`, where `z1` and `z2` are number of
 * trailing zero elements in `data1` and `data2` respectively, cyclic
 * convolution is equal to normal convolution. That means that normal
 * convolution can be calulated by applying cyclic convolution after the
 * input is zero padded with a sufficient amount of zeros.
 *
 * @param dataR - result, array of length `size`
 * @param data1 - data1, array of length `size`
 * @param data2 - data2, array of length `size`
 * @param size - number of elements, must be a power of two
 * @param root_base - a principal k-th root of unity in the ring T
 * @param root_order - order `k` of the root, must be a power of 2 not less than `size`
 */
template<typename T, typename R>
void fft_cyclic_convolution(T *dataR, T *data1, T *data2, int size, const R& root_base, int root_order) {
	R root = powT(root_base, root_order / size);
	R iroot = powT(root, size - 1); // == root^-1
	// convert to frequency domain
	fft_rec(dataR, data1, size, root); std::swap(data1, dataR);
	fft_rec(dataR, data2, size, root); std::swap(data2, dataR);
	// convolution in time domain is a pointwise
	// multiplication in frequency domain
	for (int i = 0; i < size; i++) dataR[i] = data1[i] * data2[i];
	// inverse transform is same as original transform,
	// but with inverse root and elements divided by size
	T e1 = identityT<R>::of(root);
	std::swap(data1, dataR); fft_rec(dataR, data1, size, iroot);
	T isize = e1 / T(size);
	for (int i = 0; i < size; i++) dataR[i] *= isize;
}

/**
 * FFT Ordinary Convolution of two sequences
 *
 * @param u_begin, u_end - iterators of sequence u; u_size = u_end - u_begin
 * @param v_begin, v_end - iterators of sequence v; v_size = v_end - v_begin
 * @param root_base - a principal k-th root of unity in the ring T
 * @param root_order - order `k` of the root, must be a power of 2 not less than `size`
 * @return - result, array of length size = u_size + v_size - 1
 */
template<typename T, typename R, typename It>
std::vector<T> convolution(It u_begin, It u_end, It v_begin, It v_end, const R& root_base, int root_order) {
	T e1 = identityT<R>::of(root_base), e0 = zeroT<T>::of(e1);
	std::vector<T> r, u(u_begin, u_end), v(v_begin, v_end);
	int n = (int)(u.size() + v.size() - 1);
	int l = 1; while (l < n) l *= 2;
	r.resize(l, e0); u.resize(l, e0); v.resize(l, e0);
	fft_cyclic_convolution(&r[0], &u[0], &v[0], l, root_base, root_order);
	r.resize(n);
	return r;
}

} // math
} // altruct
