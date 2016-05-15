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
template<typename T>
void fft(T *data, int size, T root) {
	T w, e1 = identityT<T>::of(root);
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
		if (j < i) swap(data[i], data[j]);
	}
}

/**
 * Fast Fourier Transform of a sequence
 *
 * Result is stored in `data`. Both `data` and `temp` are modified.
 *
 * @param data - data to transform, array of length `size`
 * @param temp - temporary storage, array of length `size`
 * @param size - number of elements, must be a power of two
 * @param root - a principal n-th root of unity in the ring T
 */
template<typename T>
void fft_rec(T *data, T *temp, int size, const T& root) {
	if (size <= 1) return;
	int h = size / 2;
	// reorder
	for (int i = 0; i < size; i++) {
		temp[i] = data[i];
	}
	for (int i = 0; i < h; i++) {
		data[i + 0] = temp[2 * i + 0];
		data[i + h] = temp[2 * i + 1];
	}
	// recurse
	T root2 = root * root;
	fft_rec(data + 0, temp, h, root2);
	fft_rec(data + h, temp, h, root2);
	// collect
	T u, v, g_i = identityT<T>::of(root);
	for (int i = 0; i < h; i++) {
		u = data[i + 0], v = data[i + h];
		v *= g_i; g_i *= root;
		data[i + 0] = u + v;
		data[i + h] = u - v;
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
 * @param root_order - order k of the root, must be a multiple of `size`
 */
template<typename T>
void fft_cyclic_convolution(T *dataR, T *data1, T *data2, int size, const T& root_base, int root_order) {
	T root = powT(root_base, root_order / size);
	// convert to frequency domain
	fft(data1, size, root);
	fft(data2, size, root);
	// convolution in time domain is a pointwise
	// multiplication in frequency domain
	for (int i = 0; i < size; i++) {
		dataR[i] = data1[i] * data2[i];
	}
	// inverse transform is same as original transform,
	// but with inverse root and elements divided by size
	T e1 = identityT<T>::of(root);
	fft(dataR, size, e1 / root);
	T ni = e1 / T(size);
	for (int i = 0; i < size; i++) {
		dataR[i] *= ni;
	}
}

/**
 * FFT Ordinary Convolution of two sequences
 *
 * @param u_begin, u_end - iterators of sequence u; u_size = u_end - u_begin
 * @param v_begin, v_end - iterators of sequence v; v_size = v_end - v_begin
 * @param root_base - a principal k-th root of unity in the ring T
 * @param root_order - order k of the root, must be a multiple of 2 not less than size
 * @return - result, array of length size = u_size + v_size - 1
 */
template<typename It, typename T = typename std::iterator_traits<It>::value_type>
std::vector<T> convolution(It u_begin, It u_end, It v_begin, It v_end, const T& root_base, int root_order) {
	T e0 = zeroT<T>::of(root_base);
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
