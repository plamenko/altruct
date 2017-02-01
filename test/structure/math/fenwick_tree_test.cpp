#include "algorithm/math/base.h"
#include "structure/math/fenwick_tree.h"

#include "gtest/gtest.h"

#include <vector>

using namespace std;
using namespace altruct::math;

TEST(fenwick_tree_test, fenwick_tree_sum) {
	int n = 10;
	fenwick_tree<int> f(n + 1, [](int r, int v){return r + v; });
	for (int i = 0; i <= n; i++) {
		f.add(i, sqT(i + 1));
	}
	f.add(5, 100);
	f.add(2, -10);
	vector<int> va(n + 1);
	for (int i = 0; i <= n; i++) {
		va[i] = f.get_sum(i);
	}
	EXPECT_EQ((vector<int>{1, 5, 4, 20, 45, 181, 230, 294, 375, 475, 596}), va);

    f.reset();
    for (int i = 0; i <= n; i++) {
        va[i] = f.get_sum(i);
    }
    EXPECT_EQ((vector<int>(n + 1, 0)), va);
}

TEST(fenwick_tree_test, fenwick_tree_max) {
    int n = 10;
    fenwick_tree<int> f(n + 1, [](int r, int v){return std::max(r, v); }, -1000000000);
    for (int i = 0; i <= n; i++) {
        f.add(i, sqT(i + 1) * ((i % 2) ? +1 : -1));
    }
    f.add(5, 50);
    f.add(2, -10);
    vector<int> va(n + 1);
    for (int i = 0; i <= n; i++) {
        va[i] = f.get_sum(i, -1000000000);
    }
    EXPECT_EQ((vector<int>{-1, 4, 4, 16, 16, 50, 50, 64, 64, 100, 100}), va);

    f.reset(-1000000000);
    for (int i = 0; i <= n; i++) {
        va[i] = f.get_sum(i, -1000000000);
    }
    EXPECT_EQ((vector<int>(n + 1, -1000000000)), va);
}
