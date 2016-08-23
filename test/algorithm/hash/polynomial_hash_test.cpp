#include "algorithm/hash/polynomial_hash.h"

#include "gtest/gtest.h"

#include <string>

using namespace std;
using namespace altruct::hash;

typedef polynomial_hash<2, int64_t> phash2;
template<> int64_t phash2::M[2] = { 758986603, 1000000007 };
template<> int64_t phash2::B[2] = { 36759071, 32547971 };
template<> int64_t phash2::BI[2] = { 366621061, 624567078 };

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
	h.append(1111, 10);
	EXPECT_EQ(10 + 1, h.W[0].size());
	EXPECT_EQ(10 + 1, h.W[1].size());
	h.append(2222, 5);
	EXPECT_EQ(10 + 1, h.W[0].size());
	EXPECT_EQ(10 + 1, h.W[1].size());
	h.append(3333, 20);
	EXPECT_EQ(20 + 1, h.W[0].size());
	EXPECT_EQ(20 + 1, h.W[1].size());
	h.append(4444, 19);
	EXPECT_EQ(20 + 1, h.W[0].size());
	EXPECT_EQ(20 + 1, h.W[1].size());
}

TEST(polynomial_hash_test, polynomial_hash_append) {
	phash2 h, h2;
	h.append(1111, 10);
	EXPECT_EQ(363095428, h.h[0]);
	EXPECT_EQ(500424796, h.h[1]);
	h.append(2222, 5);
	EXPECT_EQ(94197494, h.h[0]);
	EXPECT_EQ(596589649, h.h[1]);

	h.append(h2, 7);
	EXPECT_EQ(94197494, h.h[0]);
	EXPECT_EQ(596589649, h.h[1]);
	
	h2.append(3333, 20);
	EXPECT_EQ(461023273, h2.h[0]);
	EXPECT_EQ(648151220, h2.h[1]);
	h2.append(4444, 19);
	EXPECT_EQ(400641131, h2.h[0]);
	EXPECT_EQ(639934933, h2.h[1]);

	h.append(h2, 7);
	EXPECT_EQ(437201974, h.h[0]);
	EXPECT_EQ(641593066, h.h[1]);
}

TEST(polynomial_hash_test, polynomial_hash_subtract) {
	phash2 h1{ 2414915, 934336517 };
	phash2 h2 = h1;
	h2.append(44, 4);
	h2.append(55, 5);
	h2.subtract(h1, 4);
	EXPECT_EQ(503775743, h2.h[0]);
	EXPECT_EQ(790138442, h2.h[1]);

	phash2 h;
	h.append(44, 0);
	h.append(55, 1);
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

TEST(polynomial_hash_test, cumulative_hash) {
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
