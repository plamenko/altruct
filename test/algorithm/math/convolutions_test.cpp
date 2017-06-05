#include "altruct/algorithm/math/convolutions.h"
#include "altruct/structure/math/modulo.h"

#include <algorithm>
#include <vector>

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

TEST(convolutions_test, max_convolution) {
	typedef int64_t mod;
	const int n = 10;
	vector<mod> u{ 671, 9230, 3302, 4764, 6135, 7750, 9881, 1189, 411, 8144 };
	vector<mod> v{ 8468, 3944, 4798, 6405, 8016, 8884, 1006, 54, 7066, 3531 };
	vector<mod> r0(n);
	slow_max_convolution(r0.data(), u.data(), v.data(), n);
	vector<mod> r1(n);
	max_convolution(r1.data(), u.data(), v.data(), n);
	EXPECT_EQ(r0, r1);
	// inplace
	vector<mod> w{ 671, 9230, 3302, 4764, 6135, 7750, 9881, 1189, 411, 8144 };
	vector<mod> z0(n);
	slow_max_convolution(z0.data(), w.data(), w.data(), n);
	vector<mod> z1 = w;
	max_convolution(z1.data(), z1.data(), z1.data(), n);
	EXPECT_EQ(z0, z1);
}

TEST(convolutions_test, and_convolution) {
	typedef int64_t mod;
	const int L = 4;
	vector<mod> u{ 671, 9230, 3302, 4764, 6135, 7750, 9881, 1189, 411, 8144, 3258, 4752, 6345, 8756, 6716, 7647 };
	vector<mod> v{ 8468, 3944, 4798, 6405, 8016, 8884, 1006, 54, 7066, 3531, 6778, 9168, 7965, 6873, 6557, 2641 };
	vector<mod> r0(1 << L);
	slow_and_convolution(r0.data(), u.data(), v.data(), L);
	vector<mod> r1(1 << L);
	and_convolution(r1.data(), u.data(), v.data(), L);
	EXPECT_EQ(r0, r1);
	// inplace
	vector<mod> w{ 671, 9230, 3302, 4764, 6135, 7750, 9881, 1189, 411, 8144, 3258, 4752, 6345, 8756, 6716, 7647 };
	vector<mod> z0(1 << L);
	slow_and_convolution(z0.data(), w.data(), w.data(), L);
	vector<mod> z1 = w;
	and_convolution(z1.data(), z1.data(), z1.data(), L);
	EXPECT_EQ(z0, z1);
}

TEST(convolutions_test, or_convolution) {
	typedef int64_t mod;
	const int L = 4;
	vector<mod> u{ 671, 9230, 3302, 4764, 6135, 7750, 9881, 1189, 411, 8144, 3258, 4752, 6345, 8756, 6716, 7647 };
	vector<mod> v{ 8468, 3944, 4798, 6405, 8016, 8884, 1006, 54, 7066, 3531, 6778, 9168, 7965, 6873, 6557, 2641 };
	vector<mod> r0(1 << L);
	slow_or_convolution(r0.data(), u.data(), v.data(), L);
	vector<mod> r1(1 << L);
	or_convolution(r1.data(), u.data(), v.data(), L);
	EXPECT_EQ(r0, r1);
	// inplace
	vector<mod> w{ 671, 9230, 3302, 4764, 6135, 7750, 9881, 1189, 411, 8144, 3258, 4752, 6345, 8756, 6716, 7647 };
	vector<mod> z0(1 << L);
	slow_or_convolution(z0.data(), w.data(), w.data(), L);
	vector<mod> z1 = w;
	or_convolution(z1.data(), z1.data(), z1.data(), L);
	EXPECT_EQ(z0, z1);
}

TEST(convolutions_test, xor_convolution) {
	typedef int64_t mod;
	const int L = 4;
	vector<mod> u{ 671, 9230, 3302, 4764, 6135, 7750, 9881, 1189, 411, 8144, 3258, 4752, 6345, 8756, 6716, 7647 };
	vector<mod> v{ 8468, 3944, 4798, 6405, 8016, 8884, 1006, 54, 7066, 3531, 6778, 9168, 7965, 6873, 6557, 2641 };
	vector<mod> r0(1 << L);
	slow_xor_convolution(r0.data(), u.data(), v.data(), L);
	vector<mod> r1(1 << L);
	xor_convolution(r1.data(), u.data(), v.data(), L);
	EXPECT_EQ(r0, r1);
	// inplace
	vector<mod> w{ 671, 9230, 3302, 4764, 6135, 7750, 9881, 1189, 411, 8144, 3258, 4752, 6345, 8756, 6716, 7647 };
	vector<mod> z0(1 << L);
	slow_xor_convolution(z0.data(), w.data(), w.data(), L);
	vector<mod> z1 = w;
	xor_convolution(z1.data(), z1.data(), z1.data(), L);
	EXPECT_EQ(z0, z1);
}

TEST(convolutions_test, cyclic_convolution) {
	typedef modulo<int, 12289> mod;
	const int n = 16;
	vector<mod> u{ 671, 9230, 3302, 4764, 6135, 7750, 9881, 1189, 411, 8144, 3258, 4752, 6345, 8756, 6716, 7647 };
	vector<mod> v{ 8468, 3944, 4798, 6405, 8016, 8884, 1006, 54, 7066, 3531, 6778, 9168, 7965, 6873, 6557, 2641 };
	vector<mod> e{ 8464, 1567, 1612, 1701, 9738, 11746, 8342, 4708, 10206, 2177, 4098, 5818, 10538, 4795, 3813, 6328 };
	vector<mod> r0(n);
	slow_cyclic_convolution(r0.data(), u.data(), v.data(), n);
	EXPECT_EQ(e, r0);
}
