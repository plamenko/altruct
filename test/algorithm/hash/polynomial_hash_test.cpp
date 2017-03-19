#include "algorithm/hash/polynomial_hash.h"

#include "gtest/gtest.h"

#include <string>

using namespace std;
using namespace altruct::hash;

typedef polynomial_hash1<758986603, 36759071, 366621061> phash1;

typedef polynomial_hash<1> phashK1;
template<> int32_t phashK1::M[1] = { 758986603 };
template<> int32_t phashK1::B[1] = { 36759071 };
template<> int32_t phashK1::BI[1] = { 366621061 };

typedef polynomial_hash<2> phash2;
template<> int32_t phash2::M[2] = { 758986603, 1000000007 };
template<> int32_t phash2::B[2] = { 36759071, 32547971 };
template<> int32_t phash2::BI[2] = { 366621061, 624567078 };

TEST(polynomial_hash_test, polynomial_hash_constructor) {
	phash2 h0;
	EXPECT_EQ(0, h0.h[0]);
	EXPECT_EQ(0, h0.h[1]);

	phash2 h1{ 12 };
	EXPECT_EQ(12, h1.h[0]);
	EXPECT_EQ(0, h1.h[1]);

	phash2 h2{ 123, 456 };
	EXPECT_EQ(123, h2.h[0]);
	EXPECT_EQ(456, h2.h[1]);
	
	phash2 h3{ 1230, 4560, 7890 };
	EXPECT_EQ(1230, h3.h[0]);
	EXPECT_EQ(4560, h3.h[1]);
}

TEST(polynomial_hash_test, polynomial_hash_ensure) {
	phash2 h;
	EXPECT_EQ(0, h.W[0].size());
	EXPECT_EQ(0, h.W[1].size());
	h.add(1111, 10);
	EXPECT_EQ(10 + 1, h.W[0].size());
	EXPECT_EQ(10 + 1, h.W[1].size());
	h.add(2222, 5);
	EXPECT_EQ(10 + 1, h.W[0].size());
	EXPECT_EQ(10 + 1, h.W[1].size());
	h.add(3333, 20);
	EXPECT_EQ(20 + 1, h.W[0].size());
	EXPECT_EQ(20 + 1, h.W[1].size());
	h.add(4444, 19);
	EXPECT_EQ(20 + 1, h.W[0].size());
	EXPECT_EQ(20 + 1, h.W[1].size());
}

TEST(polynomial_hash_test, polynomial_hash_add) {
	phash2 h, h2;
	h.add(1111, 10);
	EXPECT_EQ(363095428, h.h[0]);
	EXPECT_EQ(500424796, h.h[1]);
	h.add(2222, 5);
	EXPECT_EQ(94197494, h.h[0]);
	EXPECT_EQ(596589649, h.h[1]);

	h.add(h2, 7);
	EXPECT_EQ(94197494, h.h[0]);
	EXPECT_EQ(596589649, h.h[1]);
	
	h2.add(3333, 20);
	EXPECT_EQ(461023273, h2.h[0]);
	EXPECT_EQ(648151220, h2.h[1]);
	h2.add(4444, 19);
	EXPECT_EQ(400641131, h2.h[0]);
	EXPECT_EQ(639934933, h2.h[1]);

	h.add(h2, 7);
	EXPECT_EQ(437201974, h.h[0]);
	EXPECT_EQ(641593066, h.h[1]);
}

TEST(polynomial_hash_test, polynomial_hash_sub_shr) {
	phash2 h1{ 2414915, 934336517 };
	phash2 h2 = h1;
	h2.add(44, 4);
	h2.add(55, 5);
	h2.sub_shr(h1, 4);
	EXPECT_EQ(503775743, h2.h[0]);
	EXPECT_EQ(790138442, h2.h[1]);

	phash2 h;
	h.add(44, 0);
	h.add(55, 1);
	EXPECT_EQ(503775743, h.h[0]);
	EXPECT_EQ(790138442, h.h[1]);
}

TEST(polynomial_hash_test, polynomial_hash_comparison) {
	EXPECT_FALSE((phash2{ 123, 456 }) == (phash2{ 10, 234 }));
	EXPECT_FALSE((phash2{ 123, 456 }) == (phash2{ 10, 456 }));
	EXPECT_FALSE((phash2{ 123, 456 }) == (phash2{ 10, 789 }));
	EXPECT_FALSE((phash2{ 123, 456 }) == (phash2{ 123, 234 }));
	EXPECT_TRUE((phash2{ 123, 456 }) == (phash2{ 123, 456 }));
	EXPECT_FALSE((phash2{ 123, 456 }) == (phash2{ 123, 789 }));
	EXPECT_FALSE((phash2{ 123, 456 }) == (phash2{ 999, 234 }));
	EXPECT_FALSE((phash2{ 123, 456 }) == (phash2{ 999, 456 }));
	EXPECT_FALSE((phash2{ 123, 456 }) == (phash2{ 999, 789 }));

	EXPECT_FALSE((phash2{ 123, 456 }) < (phash2{ 10, 234 }));
	EXPECT_FALSE((phash2{ 123, 456 }) < (phash2{ 10, 456 }));
	EXPECT_FALSE((phash2{ 123, 456 }) < (phash2{ 10, 789 }));
	EXPECT_FALSE((phash2{ 123, 456 }) < (phash2{ 123, 234 }));
	EXPECT_FALSE((phash2{ 123, 456 }) < (phash2{ 123, 456 }));
	EXPECT_TRUE((phash2{ 123, 456 }) < (phash2{ 123, 789 }));
	EXPECT_TRUE((phash2{ 123, 456 }) < (phash2{ 999, 234 }));
	EXPECT_TRUE((phash2{ 123, 456 }) < (phash2{ 999, 456 }));
	EXPECT_TRUE((phash2{ 123, 456 }) < (phash2{ 999, 789 }));
}

TEST(polynomial_hash_test, polynomial_hash_operators) {
	phash2 h1{ 1000000, 2000000 };
	phash2 h2{ 100000, 200000 };
	EXPECT_EQ((phash2{ 1100000, 2200000 }).h, (h1 + h2).h);
	EXPECT_EQ((phash2{ 1300000, 2300000 }).h, (h1 + 300000).h);
	EXPECT_EQ((phash2{ 900000, 1800000 }).h, (h1 - h2).h);
	EXPECT_EQ((phash2{ 700000, 1700000 }).h, (h1 - 300000).h);
	EXPECT_EQ((phash2{ 572755007, 999997207 }).h, (h1 * h2).h);
	EXPECT_EQ((phash2{ 200291815, 999995807 }).h, (h1 * 300000).h);
	EXPECT_EQ((phash2{ 416902469, 242891385 }).h, (h1 << 5).h);
	EXPECT_EQ((phash2{ 667139015, 214851611 }).h, (h1 >> 5).h);

	phash2 r;
	r = h1; r += h2;
	EXPECT_EQ((phash2{ 1100000, 2200000 }).h, r.h);
	r = h1; r += 300000;
	EXPECT_EQ((phash2{ 1300000, 2300000 }).h, r.h);
	r = h1; r -= h2;
	EXPECT_EQ((phash2{ 900000, 1800000 }).h, r.h);
	r = h1; r -= 300000;
	EXPECT_EQ((phash2{ 700000, 1700000 }).h, r.h);
	r = h1; r *= h2;
	EXPECT_EQ((phash2{ 572755007, 999997207 }).h, r.h);
	r = h1; r *= 300000;
	EXPECT_EQ((phash2{ 200291815, 999995807 }).h, r.h);
	r = h1; r <<= 5;
	EXPECT_EQ((phash2{ 416902469, 242891385 }).h, r.h);
	r = h1; r >>= 5;
	EXPECT_EQ((phash2{ 667139015, 214851611 }).h, r.h);

	r = h1; r += r;
	EXPECT_EQ((phash2{ 2000000, 4000000 }).h, r.h);
	r = h1; r -= r;
	EXPECT_EQ((phash2{ 0, 0 }).h, r.h);
	r = h1; r *= r;
	EXPECT_EQ((phash2{ 414643849, 999972007 }).h, r.h);

	r = h1; r.add(h2, 100);
	EXPECT_EQ((h1 + (h2 << 100)).h, r.h);
	r = h1; r.sub_shr(h2, 100);
	EXPECT_EQ(((h1 - h2) >> 100).h, r.h);
	r = h1; r.add(300000, 100);
	EXPECT_EQ((h1 + (phash2{ 300000, 300000 } << 100)).h, r.h);
	r = h1; r.sub_shr(300000, 100);
	EXPECT_EQ(((h1 - phash2{ 300000, 300000 }) >> 100).h, r.h);
}

TEST(polynomial_hash_test, cumulative_hash_constructor) {
	string s = "banana";
	cumulative_hash<phash2> h1(s.begin(), s.end());
	cumulative_hash<phash2> h2;
	for (const auto& c : s) h2.push_back(c);
	for (size_t b = 0; b < s.size(); b++) {
		for (size_t e = b; e <= s.size(); e++) {
			EXPECT_EQ(h1.get(b, e), h2.get(b, e)) << "[" << b << ", " << e << ")";
		}
	}
}

TEST(polynomial_hash_test, cumulative_hash_online) {
	cumulative_hash<phash2> h;
	EXPECT_EQ(0, h.size());
	EXPECT_EQ(phash2{}, h.get(0, 0));

	h.push_back(123);
	EXPECT_EQ(1, h.size());
	EXPECT_EQ(phash2{}, h.get(0, 0));
	EXPECT_EQ(phash2{}, h.get(1, 1));
	EXPECT_EQ((phash2{ 123, 123 }), h.get(0, 1));

	h.push_back(456);
	EXPECT_EQ(2, h.size());
	EXPECT_EQ(phash2{}, h.get(0, 0));
	EXPECT_EQ(phash2{}, h.get(1, 1));
	EXPECT_EQ(phash2{}, h.get(2, 2));
	EXPECT_EQ((phash2{ 123, 123 }), h.get(0, 1));
	EXPECT_EQ((phash2{ 456, 456 }), h.get(1, 2));
	EXPECT_EQ((phash2{ 64431233, 841874801 }), h.get(0, 2));

	h.push_back(789);
	EXPECT_EQ(3, h.size());
	EXPECT_EQ(phash2{}, h.get(0, 0));
	EXPECT_EQ(phash2{}, h.get(1, 1));
	EXPECT_EQ(phash2{}, h.get(2, 2));
	EXPECT_EQ(phash2{}, h.get(3, 3));
	EXPECT_EQ((phash2{ 123, 123 }), h.get(0, 1));
	EXPECT_EQ((phash2{ 456, 456 }), h.get(1, 2));
	EXPECT_EQ((phash2{ 789, 789 }), h.get(2, 3));
	EXPECT_EQ((phash2{ 64431233, 841874801 }), h.get(0, 2));
	EXPECT_EQ((phash2{ 161416561, 680349400 }), h.get(1, 3));
	EXPECT_EQ((phash2{ 90981281, 386059579 }), h.get(0, 3));

	h.pop_back();
	EXPECT_EQ(2, h.size());
	EXPECT_EQ(phash2{}, h.get(0, 0));
	EXPECT_EQ(phash2{}, h.get(1, 1));
	EXPECT_EQ(phash2{}, h.get(2, 2));
	EXPECT_EQ((phash2{ 123, 123 }), h.get(0, 1));
	EXPECT_EQ((phash2{ 456, 456 }), h.get(1, 2));
	EXPECT_EQ((phash2{ 64431233, 841874801 }), h.get(0, 2));

	h.pop_back();
	EXPECT_EQ(1, h.size());
	EXPECT_EQ(phash2{}, h.get(0, 0));
	EXPECT_EQ(phash2{}, h.get(1, 1));
	EXPECT_EQ((phash2{ 123, 123 }), h.get(0, 1));

	h.pop_back();
	EXPECT_EQ(0, h.size());
	EXPECT_EQ(phash2{}, h.get(0, 0));
}

TEST(polynomial_hash_test, polynomial_hash1) {
    vector<int> v;
    int x = 31;
    for (int i = 0; i < 1000; i++) {
        v.push_back(x = x * int64_t(x) % 997230937);
    }
    cumulative_hash<phash1> ch1(v.begin(), v.end());
    cumulative_hash<phashK1> chk1(v.begin(), v.end());
    for (int e = 0; e <= v.size(); e++) {
        for (int b = 0; b <= e; b++) {
            EXPECT_EQ(chk1.get(b, e).hash(), ch1.get(b, e).hash()) << " b = " << b << ", e = " << e;
        }
    }
}
