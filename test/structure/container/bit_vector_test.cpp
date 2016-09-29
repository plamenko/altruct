#include "algorithm/collections/collections.h"
#include "algorithm/math/base.h"
#include "algorithm/random/xorshift.h"
#include "structure/container/bit_vector.h"

#include "structure_test_util.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;
using namespace altruct::collections;
using namespace altruct::container;
using namespace altruct::random;
using namespace altruct::test_util;

namespace {
string vec_to_string(const vector<int>& a) {
	string s(a.begin(), a.end());
	for (char& c : s) c += '0';
	return s;
}

vector<int> random_vec(int size, int seed = 12345) {
	xorshift_64star rng(seed);
	vector<int> a(size);
	for (auto& e : a) {
		e = rng.next() % 2;
	}
	return a;
}

template<typename W>
void test_scan(const bit_vector<W>& v, size_t begin, size_t end) {
	vector<int> expected_values;
	vector<int> expected_positions;
	for (auto it = begin; it != end; ++it) {
		expected_positions.push_back(int(it - begin));
		expected_values.push_back(v.bit_at(it));
	}
	vector<int> actual_values;
	vector<int> actual_positions;
	v.scan(begin, end, [&](W w, size_t pos, int l){
		while (l > 0) {
			actual_positions.push_back(int(pos));
			actual_values.push_back(w & 1);
			pos++; l--; w >>= 1;
		}
		return true;
	});
	ASSERT_EQ(expected_positions, actual_positions);
	ASSERT_EQ(expected_values, actual_values);
}

template<typename W>
void test_scan(const bit_vector<W>& v1, size_t begin1, const bit_vector<W>& v2, size_t begin2, size_t len) {
	vector<int> expected_values1;
	vector<int> expected_values2;
	for (size_t pos = 0; pos < len; pos++) {
		expected_values1.push_back(v1.bit_at(begin1 + pos));
		expected_values2.push_back(v2.bit_at(begin2 + pos));
	}
	vector<int> actual_values1;
	vector<int> actual_values2;
	bit_vector<W>::scan(v1, begin1, v2, begin2, len, [&](W w1, W w2, int l){
		while (l > 0) {
			actual_values1.push_back(w1 & 1);
			actual_values2.push_back(w2 & 1);
			l--; w1 >>= 1; w2 >>= 1;
		}
		return true;
	});
	ASSERT_EQ(expected_values1, actual_values1);
	ASSERT_EQ(expected_values2, actual_values2);
}

template<typename W, typename F>
void test_apply(const vector<int>& a, size_t begin, size_t end, F op) {
	vector<int> a_expected = a;
	for (size_t pos = begin; pos != end; pos++) {
		a_expected[pos] = op(a[pos], pos, 1);
	}
	bit_vector<W> v_expected(a_expected.begin(), a_expected.end());
	bit_vector<W> v_actual(a.begin(), a.end());
	v_actual.apply(begin, end, op);
	ASSERT_EQ(v_expected.words, v_actual.words);
}

template<typename W, typename F>
void test_apply(const vector<int>& a, size_t begin1, size_t begin2, size_t len, F op) {
	vector<int> a_expected = a;
	for (size_t pos = 0; pos < len; pos++) {
		a_expected[begin1 + pos] = op(a[begin1 + pos], a[begin2 + pos], 1);
	}
	bit_vector<W> v_expected(a_expected.begin(), a_expected.end());
	bit_vector<W> v_actual(a.begin(), a.end());
	// use v_actual both as destination and source in order to test iteration direction;
	// since v_actual is both reading and overwriting itself, direction matters;
	bit_vector<W>::apply(v_actual, begin1, v_actual, begin2, len, op);
	ASSERT_EQ(v_expected.words, v_actual.words);
}

}

TEST(bit_vector_test, constructor_64bit) {
	bit_vector<uint64_t> v0;
	EXPECT_EQ(0, v0.size());
	EXPECT_EQ((vector<uint64_t>{0x0000000000000000}), v0.words);

	bit_vector<uint64_t> v1(1);
	EXPECT_EQ(1, v1.size());
	EXPECT_EQ((vector<uint64_t>{0x0000000000000000, 0x0000000000000000}), v1.words);

	bit_vector<uint64_t> v2(63);
	EXPECT_EQ(63, v2.size());
	EXPECT_EQ((vector<uint64_t>{0x0000000000000000, 0x0000000000000000}), v2.words);

	bit_vector<uint64_t> v3(64);
	EXPECT_EQ(64, v3.size());
	EXPECT_EQ((vector<uint64_t>{0x0000000000000000, 0x0000000000000000}), v3.words);

	bit_vector<uint64_t> v4(65);
	EXPECT_EQ(65, v4.size());
	EXPECT_EQ((vector<uint64_t>{0x0000000000000000, 0x0000000000000000, 0x0000000000000000}), v4.words);

	vector<int> a0;
	bit_vector<uint64_t> v5(a0.begin(), a0.end());
	EXPECT_EQ(0, v5.size());
	EXPECT_EQ((vector<uint64_t>{0x0000000000000000}), v5.words);

	vector<int> a1{ 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1 };
	bit_vector<uint64_t> v6(a1.begin(), a1.end());
	EXPECT_EQ(20, v6.size());
	EXPECT_EQ((vector<uint64_t>{0x00000000000AC145, 0x0000000000000000}), v6.words);

	bit_vector<uint64_t> v7{ 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1 };
	EXPECT_EQ(20, v7.size());
	EXPECT_EQ((vector<uint64_t>{0x00000000000AC145, 0x0000000000000000}), v7.words);

	vector<int> p2{ 0, 1, 297 };
	vector<int> a2(298); for (int i : p2) a2[i] = 1;
	bit_vector<uint64_t> v8(a2.begin(), a2.end());
	EXPECT_EQ(298, v8.size());
	EXPECT_EQ((vector<uint64_t>{0x00000000000000003, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000020000000000, 0x0000000000000000}), v8.words);

	auto v9 = v6.clone();
	EXPECT_EQ(20, v9.size());
	EXPECT_EQ((vector<uint64_t>{0x00000000000AC145, 0x0000000000000000}), v9.words);
}

TEST(bit_vector_test, constructor_8bit) {
	bit_vector<uint8_t> v0;
	EXPECT_EQ(0, v0.size());
	EXPECT_EQ((vector<uint8_t>{0x00}), v0.words);

	bit_vector<uint8_t> v1(1);
	EXPECT_EQ(1, v1.size());
	EXPECT_EQ((vector<uint8_t>{0x00, 0x00}), v1.words);

	bit_vector<uint8_t> v2(7);
	EXPECT_EQ(7, v2.size());
	EXPECT_EQ((vector<uint8_t>{0x00, 0x00}), v2.words);

	bit_vector<uint8_t> v3(8);
	EXPECT_EQ(8, v3.size());
	EXPECT_EQ((vector<uint8_t>{0x00, 0x00}), v3.words);

	bit_vector<uint8_t> v4(9);
	EXPECT_EQ(9, v4.size());
	EXPECT_EQ((vector<uint8_t>{0x00, 0x00, 0x00}), v4.words);

	vector<int> a0;
	bit_vector<uint8_t> v5(a0.begin(), a0.end());
	EXPECT_EQ(0, v5.size());
	EXPECT_EQ((vector<uint8_t>{0x00}), v5.words);

	vector<int> a1{ 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1 };
	bit_vector<uint8_t> v6(a1.begin(), a1.end());
	EXPECT_EQ(20, v6.size());
	EXPECT_EQ((vector<uint8_t>{0x45, 0xC1, 0x0A, 0x00}), v6.words);

	bit_vector<uint8_t> v7{ 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1 };
	EXPECT_EQ(20, v7.size());
	EXPECT_EQ((vector<uint8_t>{0x45, 0xC1, 0x0A, 0x00}), v7.words);

	vector<int> p2{ 0, 1, 97 };
	vector<int> a2(98); for (int i : p2) a2[i] = 1;
	bit_vector<uint8_t> v8(a2.begin(), a2.end());
	EXPECT_EQ(98, v8.size());
	EXPECT_EQ((vector<uint8_t>{0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0x00}), v8.words);

	auto v9 = v6.clone();
	EXPECT_EQ(20, v9.size());
	EXPECT_EQ((vector<uint8_t>{0x45, 0xC1, 0x0A, 0x00}), v9.words);
}

TEST(bit_vector_test, resize) {
	bit_vector<uint64_t> v6{ 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1 };
	v6.resize(19);
	EXPECT_EQ(19, v6.size());
	EXPECT_EQ((vector<uint64_t>{0x000000000002C145, 0x0000000000000000}), v6.words);
	v6.resize(17);
	EXPECT_EQ(17, v6.size());
	EXPECT_EQ((vector<uint64_t>{0x000000000000C145, 0x0000000000000000}), v6.words);
	v6.resize(70);
	EXPECT_EQ(70, v6.size());
	EXPECT_EQ((vector<uint64_t>{0x000000000000C145, 0x0000000000000000, 0x0000000000000000}), v6.words);
	v6.words[1] = 0xFFFFFFFFFFFFFFFF;
	EXPECT_EQ(70, v6.size());
	EXPECT_EQ((vector<uint64_t>{0x000000000000C145, 0xFFFFFFFFFFFFFFFF, 0x0000000000000000}), v6.words);
	v6.resize(65);
	EXPECT_EQ(65, v6.size());
	EXPECT_EQ((vector<uint64_t>{0x000000000000C145, 0x0000000000000001, 0x0000000000000000}), v6.words);
	v6.resize(64);
	EXPECT_EQ(64, v6.size());
	EXPECT_EQ((vector<uint64_t>{0x000000000000C145, 0x0000000000000000}), v6.words);
}

TEST(bit_vector_test, reserve) {
	bit_vector<uint64_t> v6{ 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1 };
	v6.reserve(10);
	EXPECT_EQ(20, v6.size());
	EXPECT_EQ((vector<uint64_t>{0x00000000000AC145, 0x0000000000000000}), v6.words);
	v6.reserve(50);
	EXPECT_EQ(50, v6.size());
	EXPECT_EQ((vector<uint64_t>{0x00000000000AC145, 0x0000000000000000}), v6.words);
	v6.reserve(128);
	EXPECT_EQ(128, v6.size());
	EXPECT_EQ((vector<uint64_t>{0x00000000000AC145, 0x0000000000000000, 0x0000000000000000}), v6.words);
}

TEST(bit_vector_test, getters_setters) {
	auto a = random_vec(1000);
	bit_vector<> v(a.begin(), a.end());
	EXPECT_EQ(1000, v.size());
	for (int i = 0; i < v.size(); i++) {
		EXPECT_EQ(a[i], v.bit_at(i)) << i;
	}
	uint64_t w = 0;
	EXPECT_EQ(1000, v.size());
	for (int i = (int)v.size() - 1; i >= 0; i--) {
		w <<= 1, w |= a[i];
		EXPECT_EQ(w, v.word_at(i)) << i;
	}
	bit_vector<> v2(1000);
	EXPECT_EQ(1000, v2.size());
	for (int i = 0; i < v2.size(); i++) {
		v2.set(i, a[i]);
	}
	EXPECT_EQ(1000, v2.size());
	for (int i = 0; i < v2.size(); i++) {
		EXPECT_EQ(a[i], v2.bit_at(i)) << i;
	}
}

TEST(bit_vector_test, bit_proxy) {
	bit_vector<> bv(10);

	bv[4] = 0; // 0 => 0
	EXPECT_EQ(0, (int)bv[4]);
	bv[4] = 1; // 0 => 1
	EXPECT_EQ(1, (int)bv[4]);
	bv[4] = 1; // 1 => 1
	EXPECT_EQ(1, (int)bv[4]);
	bv[4] = 0; // 1 => 0
	EXPECT_EQ(0, (int)bv[4]);
	
	bv[4] ^= 0; // 0 ^ 0 => 0
	EXPECT_EQ(0, (int)bv[4]);
	bv[4] ^= 1; // 0 ^ 1 => 1
	EXPECT_EQ(1, (int)bv[4]);
	bv[4] ^= 0; // 1 ^ 0 => 1
	EXPECT_EQ(1, (int)bv[4]);
	bv[4] ^= 1; // 1 ^ 1 => 0
	EXPECT_EQ(0, (int)bv[4]);
	
	bv[4] |= 0; // 0 | 0 => 0
	EXPECT_EQ(0, (int)bv[4]);
	bv[4] |= 1; // 0 | 1 => 1
	EXPECT_EQ(1, (int)bv[4]);
	bv[4] |= 1; // 1 | 1 => 1
	EXPECT_EQ(1, (int)bv[4]);
	bv[4] |= 0; // 1 | 0 => 1
	EXPECT_EQ(1, (int)bv[4]);

	bv[4] &= 1; // 1 & 1 => 1
	EXPECT_EQ(1, (int)bv[4]);
	bv[4] &= 0; // 1 & 0 => 0
	EXPECT_EQ(0, (int)bv[4]);
	bv[4] &= 1; // 0 & 1 => 0
	EXPECT_EQ(0, (int)bv[4]);
	bv[4] &= 0; // 0 & 0 => 0
	EXPECT_EQ(0, (int)bv[4]);

	bv[5] = 1; bv[6] = 0;
	bv[7] = !bv[5] | bv[6];
	EXPECT_EQ(0, (int)bv[7]);
	bv[7] = !bv[6] | bv[7];
	EXPECT_EQ(1, (int)bv[7]);
	bv[7] = !bv[7] & bv[5];
	EXPECT_EQ(0, (int)bv[7]);
}

TEST(bit_vector_test, to_string) {
	auto a = random_vec(100); auto s = vec_to_string(a);
	bit_vector<uint8_t> v(a.begin(), a.end());
	for (size_t b = 0; b <= a.size(); b++) {
		for (size_t e = b; e <= a.size(); e++) {
			EXPECT_EQ(s.substr(b, e - b), v.to_string(b, e)) << b << " " << e;
		}
	}
}

TEST(bit_vector_test, scan) {
	auto a1 = random_vec(40);
	auto a2 = random_vec(30);
	bit_vector<uint8_t> v1(a1.begin(), a1.end());
	bit_vector<uint8_t> v2(a2.begin(), a2.end());
	for (size_t b1 = 0; b1 <= v1.size(); b1++) {
		for (size_t e1 = b1; e1 <= v1.size(); e1++) {
			test_scan(v1, b1, e1);
		}
	}
	for (size_t b1 = 0; b1 <= v1.size(); b1++) {
		for (size_t b2 = 0; b2 <= v2.size(); b2++) {
			for (size_t len = 0; len <= min(v1.size() - b1, v2.size() - b2); len++) {
				test_scan(v1, b1, v2, b2, len);
			}
		}
	}
}

TEST(bit_vector_test, apply) {
	auto a1 = random_vec(40);
	for (size_t b1 = 0; b1 <= a1.size(); b1++) {
		for (size_t e1 = b1; e1 <= a1.size(); e1++) {
			test_apply<uint8_t>(a1, b1, e1, bit_vector<uint8_t>::op_set0);
			test_apply<uint8_t>(a1, b1, e1, bit_vector<uint8_t>::op_set1);
			test_apply<uint8_t>(a1, b1, e1, bit_vector<uint8_t>::op_flip);
		}
	}
	for (size_t b1 = 0; b1 <= a1.size(); b1++) {
		for (size_t b2 = 0; b2 <= a1.size(); b2++) {
			for (size_t len = 0; len <= min(a1.size() - b1, a1.size() - b2); len++) {
				test_apply<uint8_t>(a1, b1, b2, len, bit_vector<uint8_t>::op_set);
				test_apply<uint8_t>(a1, b1, b2, len, bit_vector<uint8_t>::op_and);
				test_apply<uint8_t>(a1, b1, b2, len, bit_vector<uint8_t>::op_or);
				test_apply<uint8_t>(a1, b1, b2, len, bit_vector<uint8_t>::op_xor);
			}
		}
	}
}

TEST(bit_vector_test, compare) {
	auto a1 = random_vec(40);
	auto a2 = random_vec(30);
	bit_vector<uint8_t> v1(a1.begin(), a1.end());
	bit_vector<uint8_t> v2(a2.begin(), a2.end());
	for (size_t b1 = 0; b1 <= v1.size(); b1++) {
		for (size_t e1 = b1; e1 <= v1.size(); e1++) {
			for (size_t b2 = 0; b2 <= v2.size(); b2++) {
				for (size_t e2 = b2; e2 <= v2.size(); e2++) {
					int r_expected = compare(a1.data() + b1, a1.data() + e1, a2.data() + b2, a2.data() + e2);
					int r_actual = bit_vector<uint8_t>::compare(v1, b1, e1, v2, b2, e2);
					EXPECT_EQ(r_expected, r_actual) << b1 << " " << e1 << " " << b2 << " " << e2;
				}
			}
		}
	}
}

TEST(bit_vector_test, reverse) {
	auto a = random_vec(100);
	bit_vector<uint8_t> v(a.begin(), a.end());
	for (size_t b = 0; b <= a.size(); b++) {
		for (size_t e = b; e <= a.size(); e++) {
			reverse(a.data() + b, a.data() + e);
			bit_vector<uint8_t> v_expected(a.begin(), a.end());
			v.reverse(b, e);
			EXPECT_EQ(v_expected.words, v.words);
		}
	}
}

TEST(bit_vector_test, rotate) {
	auto a = random_vec(40);
	for (size_t b = 0; b <= a.size(); b++) {
		for (size_t m = b; m <= a.size(); m++) {
			for (size_t e = m; e <= a.size(); e++) {
				bit_vector<uint8_t> v(a.begin(), a.end());
				bit_vector<uint8_t> v_expected0(a.begin(), a.end());
				rotate(a.data() + b, a.data() + m, a.data() + e);
				bit_vector<uint8_t> v_expected1(a.begin(), a.end());
				v.rotate_left(b, e, m - b);
				EXPECT_EQ(v_expected1.words, v.words);
				v.rotate_right(b, e, m - b);
				EXPECT_EQ(v_expected0.words, v.words);
			}
		}
	}
}

TEST(bit_vector_test, swap) {
	// TODO
}

TEST(bit_vector_test, hamming) {
	// TODO
}

TEST(bit_vector_test, comparison_operators) {
	bit_vector<uint8_t> v1{ 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1 };
	bit_vector<uint8_t> v2{ 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1 };
	bit_vector<uint8_t> v3{ 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1 };
	ASSERT_COMPARISON_OPERATORS(0, v1, v1);
	ASSERT_COMPARISON_OPERATORS(0, v2, v2);
	ASSERT_COMPARISON_OPERATORS(0, v3, v3);
	ASSERT_COMPARISON_OPERATORS(-1, v1, v2);
	ASSERT_COMPARISON_OPERATORS(+1, v2, v1);
	ASSERT_COMPARISON_OPERATORS(-1, v1, v3);
	ASSERT_COMPARISON_OPERATORS(+1, v3, v1);
	ASSERT_COMPARISON_OPERATORS(+1, v2, v3);
	ASSERT_COMPARISON_OPERATORS(-1, v3, v2);
}

TEST(bit_vector_test, logic_operators) {
	// TODO
}

TEST(bit_vector_test, perf) {
	return; // do not test perf by default
	int n = 50000;
	bit_vector<> v(n);
	uint64_t r = 0;
	auto T0 = clock();
	for (int i = 0; i < 75000; i++) {
		int b = i % 100, c = i % 101, e = n - i % 103;
		//v.apply(b, e, v.op_set0); // 16 ms
		//v.apply(b, e, v.op_set1); // 16 ms
		//v.apply(b, e, v.op_flip); // 66 ms
		//v.reverse(b, e);                  // 315 ms
		//r += v.hamming_distance(v, b, v, c, e - 150); // 460 ms
		//v.rotate_left(b, e, n / 3);       // 630 ms
		//v.rotate_right(b, e, n / 3);      // 630 ms
		//v.swap(b, n/3, n*2/3, e);         // 630 ms
	}
	fprintf(stderr, "%d ms    %llx\n", clock() - T0, v.words[0] + r);
}
