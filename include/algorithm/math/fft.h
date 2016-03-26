#pragma once

#include "base.h"

namespace altruct {
namespace math {

/**
 * Fast Fourier Transform of a sequence
 *
 * Result is stored in `data`.
 * Both `data` and `temp` array are modified.
 *
 * @param data - data to transform, array of length `size`
 * @param temp - temporary storage, array of length `size`
 * @param size - number of elements, must be a power of two
 * @param root - a principal n-th root of unity in the ring T
 */
template<typename T>
void fft(T *data, T *temp, int size, const T& root) {
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
	fft(data + 0, temp, h, root2);
	fft(data + h, temp, h, root2);
	// collect
	T u, v, g_i = 1;
	for (int i = 0; i < h; i++) {
		u = data[i + 0], v = data[i + h];
		v *= g_i; g_i *= root;
		data[i + 0] = u + v;
		data[i + h] = u - v;
	}
}

/**
 * FFT Convolution of two sequences
 *
 * Result is stored in `dataR`.
 * All `dataR`, `data1` and `data2` array are modified.
 *
 * @param dataR - result, array of length `size`
 * @param data1 - data1, array of length `size`
 * @param data2 - data2, array of length `size`
 * @param size - number of elements, must be a power of two
 * @param root_base - a principal k-th root of unity in the ring T
 * @param root_order - order k of the root, must be a multiple of `size`
 */
template<typename T>
void fft_convolve(T *dataR, T *data1, T *data2, int size, const T& root_base, int root_order) {
	T root = powT(root_base, root_order / size);
	// convert to frequency domain
	fft(data1, dataR, size, root);
	fft(data2, dataR, size, root);
	// convolution in time domain is a pointwise
	// multiplication in frequency domain
	for (int i = 0; i < size; i++) {
		dataR[i] = data1[i] * data2[i];
	}
	// inverse transform is same as original transform,
	// but with inverse root and elements divided by size
	T e1 = identityT<T>::of(root);
	fft(dataR, data1, size, e1 / root);
	T ni = e1 / T(size);
	for (int i = 0; i < size; i++) {
		dataR[i] *= ni;
	}
}

} // math
} // altruct
