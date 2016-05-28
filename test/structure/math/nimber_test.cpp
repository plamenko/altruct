#include "structure/math/nimber.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

typedef nimber<int> nim;

TEST(nimber_test, constructor) {
	nim n1;
	EXPECT_EQ(0, n1.v);
	nim n2(5);
	EXPECT_EQ(5, n2.v);
	nim n3(n2);
	EXPECT_EQ(5, n3.v);
}

TEST(nimber_test, operators_comparison) {
	EXPECT_TRUE(nim(2) <  nim(5));
	EXPECT_TRUE(nim(2) <= nim(5));
	EXPECT_FALSE(nim(2) >  nim(5));
	EXPECT_FALSE(nim(2) >= nim(5));
	EXPECT_FALSE(nim(2) == nim(5));
	EXPECT_TRUE(nim(2) != nim(5));

	EXPECT_FALSE(nim(5) <  nim(5));
	EXPECT_TRUE(nim(5) <= nim(5));
	EXPECT_FALSE(nim(5) >  nim(5));
	EXPECT_TRUE(nim(5) >= nim(5));
	EXPECT_TRUE(nim(5) == nim(5));
	EXPECT_FALSE(nim(5) != nim(5));

	EXPECT_FALSE(nim(5) <  nim(2));
	EXPECT_FALSE(nim(5) <= nim(2));
	EXPECT_TRUE(nim(5) >  nim(2));
	EXPECT_TRUE(nim(5) >= nim(2));
	EXPECT_FALSE(nim(5) == nim(2));
	EXPECT_TRUE(nim(5) != nim(2));
}

TEST(nimber_test, operators_arithmetic) {
	const nim n1(3);
	const nim n2(10);
	EXPECT_EQ(nim(9), n1 + n2);
	EXPECT_EQ(nim(9), n1 - n2);
	EXPECT_EQ(nim(3), -n1);
	EXPECT_EQ(nim(5), n1 * n2);
	EXPECT_EQ(nim(4), n1 / n2);
	EXPECT_EQ(nim(0), n1 % n2);
	EXPECT_EQ(nim(9), n2 + n1);
	EXPECT_EQ(nim(9), n2 - n1);
	EXPECT_EQ(nim(10), -n2);
	EXPECT_EQ(nim(5), n2 * n1);
	EXPECT_EQ(nim(15), n2 / n1);
}

TEST(nimber_test, operators_inplace) {
	const nim n1(3);
	const nim n2(10);
	nim nr;
	nr = n1; nr += n2;
	EXPECT_EQ(nim(9), nr);
	nr = n1; nr -= n2;
	EXPECT_EQ(nim(9), nr);
	nr = n1; nr *= n2;
	EXPECT_EQ(nim(5), nr);
	nr = n1; nr /= n2;
	EXPECT_EQ(nim(4), nr);
	nr = n1; nr %= n2;
	EXPECT_EQ(nim(0), nr);
	nr = n2; nr += n1;
	EXPECT_EQ(nim(9), nr);
	nr = n2; nr -= n1;
	EXPECT_EQ(nim(9), nr);
	nr = n2; nr *= n1;
	EXPECT_EQ(nim(5), nr);
	nr = n2; nr /= n1;
	EXPECT_EQ(nim(15), nr);
	nr = n2; nr %= n1;
	EXPECT_EQ(nim(0), nr);
}

TEST(nimber_test, operators_inplace_self) {
	const nim n1(13);
	nim nr;
	nr = n1; nr += nr;
	EXPECT_EQ(nim(0), nr);
	nr = n1; nr -= nr;
	EXPECT_EQ(nim(0), nr);
	nr = n1; nr *= nr;
	EXPECT_EQ(nim(10), nr);
	nr = n1; nr /= nr;
	EXPECT_EQ(nim(1), nr);
	nr = n1; nr %= nr;
	EXPECT_EQ(nim(0), nr);
}

TEST(nimber_test, identity) {
	const nim n1(5);
	const nim e0 = zeroT<nim>::of(n1);
	const nim e1 = identityT<nim>::of(n1);
	EXPECT_EQ(0, e0.v);
	EXPECT_EQ(1, e1.v);
}

TEST(nimber_test, inverse) {
	vector<nim> vi; for (int v = 0; v <= 30; v++) vi.push_back(nim(v).inverse());
	EXPECT_EQ((vector<nim>{0, 1, 3, 2, 15, 12, 9, 11, 10, 6, 8, 7, 5, 14, 13, 4, 170, 160, 109, 107, 131, 139, 116, 115, 228, 234, 92, 89, 73, 77, 220}), vi);
	//for (int v = 1; v <= 1000000; v++) {
	//  nim i = nim(v).inverse();
	//	EXPECT_EQ(nim(1), nim(v) * i) << v;
	//}
}

TEST(nimber_test, sqrt) {
	vector<nim> vq; for (int v = 0; v <= 30; v++) vq.push_back(nim(v).sqrt());
	EXPECT_EQ((vector<nim>{0, 1, 3, 2, 7, 6, 4, 5, 14, 15, 13, 12, 9, 8, 10, 11, 30, 31, 29, 28, 25, 24, 26, 27, 16, 17, 19, 18, 23, 22, 20}), vq);
	//for (int v = 1; v <= 1000000; v++) {
	//	nim q = nim(v).sqrt();
	//	EXPECT_EQ(nim(v), q * q) << v;
	//}
}

TEST(nimber_test, mul_perf) {
	return; // skip perf test
	int T0 = clock();
	int sz = 10000;
	for (int a = 0; a <= sz; a++) {
		for (int b = 0; b <= sz; b++) {
			nim m1 = nim(a) * nim(b);
			nim m2 = nim::mul2(nim(a), nim(b));
			EXPECT_EQ(m1.v, m2.v) << "a: " << a << ", b: " << b;
		}
	}
	printf("%d ms\n", clock() - T0);

	int T1 = clock();
	nim m1;
	for (int a = 0; a <= sz; a++) {
		for (int b = 0; b <= sz; b++) {
			m1 = nim(a) * nim(b);
		}
	}
	printf("m1: %d ms  %d\n", clock() - T1, m1);

	int T2 = clock();
	nim m2;
	for (int a = 0; a <= sz; a++) {
		for (int b = 0; b <= sz; b++) {
			m2 = nim::mul2(nim(a), nim(b));
		}
	}
	printf("m2: %d ms  %d\n", clock() - T2, m2);
}

TEST(nimber_test, nim8) {
	return; // skip slow test
	typedef nimber<int8_t> nim8;
	for (int a = 0; a < 256; a++) {
		for (int b = 0; b < 256; b++) {
			nim m0 = nim(a) * nim(b);
			nim8 m1 = nim8(a) * nim8(b);
			EXPECT_EQ(m0.v, (int)m1.v) << "a: " << a << ", b: " << b;
			nim8 m2 = nim8::mul2(nim8(a), nim8(b));
			EXPECT_EQ(m0.v, (int)m2.v) << "a: " << a << ", b: " << b;
			nim8 d1 = nim8(a) / nim8(b);
			if (b != 0) EXPECT_EQ(nim8(a), d1 * nim8(b)) << "a: " << a << ", b: " << b;
			nim8 d2 = nim8(b) / nim8(a);
			if (a != 0) EXPECT_EQ(nim8(b), d2 * nim8(a)) << "a: " << a << ", b: " << b;
		}
		nim8 ai = nim8(a).inverse();
		if (a != 0) EXPECT_EQ(nim8(1), ai * nim8(a)) << "a: " << a;
		nim8 aq = nim8(a).sqrt();
		EXPECT_EQ(nim8(a), aq * aq) << "a: " << a;
	}
}
