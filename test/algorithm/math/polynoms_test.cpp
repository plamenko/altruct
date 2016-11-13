#include "algorithm/math/polynoms.h"
#include "structure/math/fraction.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

double round(double val, double mul) {
	return round(val / mul) * mul;
}

vector<double> round(vector<double> vz, double mul) {
	for (auto &z : vz) {
		z = round(z, mul);
	}
	return vz;
}

TEST(polynoms_test, search) {
	const polynom<int> p1{ 7, -5, -13, 4 };
	EXPECT_DOUBLE_EQ(-0.8262501959871101, monotonic_search(p1, -1e20, -0.177, 0.0));
	EXPECT_DOUBLE_EQ(+0.6112574125565371, monotonic_search(p1, -0.177, +2.344, 0.0));
	EXPECT_DOUBLE_EQ(+3.4649927834305730, monotonic_search(p1, +2.344, +1e20, 0.0));
}

TEST(polynoms_test, zeros) {
	const polynom<int> p1{ 7, -5, -13, 4 };
	const polynom<int> p2{ 70, -5, -13, 4 };
	const polynom<int> p3{ -12, 16, -7, 1 };
	vector<double> vz1 = round(find_zeros(p1, 1e20, 1e-12), 1e-12);
	EXPECT_EQ((vector<double>{-0.826250195987, 0.611257412557, 3.464992783431}), vz1);
	vector<double> vz2 = round(find_zeros(p2, 1e20, 1e-12), 1e-12);
	EXPECT_EQ((vector<double>{-1.957184056592}), vz2);
	vector<double> vz3 = round(find_zeros(p3, 1e10, 1e-12), 1e-12);
	EXPECT_EQ((vector<double>{2.0, 2.0, 3.0}), vz3);
}

TEST(polynoms_test, polynom_sum) {
	typedef fraction<int> frac;
	EXPECT_EQ((polynom<frac>{0, 1} / frac(1)), polynom_sum(polynom<frac>{ 1 }));
	EXPECT_EQ((polynom<frac>{0, 1, 1} / frac(2)), polynom_sum(polynom<frac>{ 0, 1 }));
	EXPECT_EQ((polynom<frac>{0, 1, 3, 2} / frac(6)), polynom_sum(polynom<frac>{ 0, 0, 1 }));
	EXPECT_EQ((polynom<frac>{0, 19, 15, 14} / frac(6)), polynom_sum(polynom<frac>{ 3, -2, 7 }));
}
