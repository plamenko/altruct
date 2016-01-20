#include "structure/math/matrix.h"
#include "structure/math/modulo.h"

#include "gtest/gtest.h"

#include <vector>

using namespace std;
using namespace altruct::math;

typedef modulo<int, 1000000007> mod;

TEST(matrix_test, constructor) {
	matrix<int> m1;
	EXPECT_EQ(0, m1.rows());
	EXPECT_EQ(0, m1.cols());
	
	matrix<int> m2(3);
	EXPECT_EQ(3, m2.rows());
	EXPECT_EQ(3, m2.cols());
	EXPECT_EQ((vector<vector<int>>{ { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } }), m2.a);
	
	matrix<int> m3(3, 2);
	EXPECT_EQ(3, m3.rows());
	EXPECT_EQ(2, m3.cols());
	EXPECT_EQ((vector<vector<int>>{ { 0, 0 }, { 0, 0 }, { 0, 0 } }), m3.a);

	matrix<int> m4({
		{ 00, 01, 02, 03, 04, 05 },
		{ 10, 11, 12, 13, 14, 15 },
		{ 20, 21, 22, 23, 24, 25 },
		{ 30, 31, 32, 33, 34, 35 },
		{ 40, 41, 42, 43, 44, 45 }});
	vector<vector<int>> a4 {
		{ 00, 01, 02, 03, 04, 05 },
		{ 10, 11, 12, 13, 14, 15 },
		{ 20, 21, 22, 23, 24, 25 },
		{ 30, 31, 32, 33, 34, 35 },
		{ 40, 41, 42, 43, 44, 45 } };
	EXPECT_EQ(5, m4.rows());
	EXPECT_EQ(6, m4.cols());
	EXPECT_EQ(a4, m4.a);

	matrix<int> m5(m4);
	EXPECT_EQ(5, m5.rows());
	EXPECT_EQ(6, m5.cols());
	EXPECT_EQ(a4, m5.a);

	matrix<int> m6(m4, 2, 1);
	vector<vector<int>> a6{
		{ 21, 22, 23, 24, 25 },
		{ 31, 32, 33, 34, 35 },
		{ 41, 42, 43, 44, 45 } };
	EXPECT_EQ(3, m6.rows());
	EXPECT_EQ(5, m6.cols());
	EXPECT_EQ(a6, m6.a);

	matrix<int> m7(m4, 2, 1, 2, 3);
	EXPECT_EQ(2, m7.rows());
	EXPECT_EQ(3, m7.cols());
	EXPECT_EQ((vector<vector<int>>{{ 21, 22, 23 }, { 31, 32, 33 } }), m7.a);
}

TEST(modulo_test, swap) {
	matrix<int> m1({ { 1, 2, 3 }, { 4, 5, 6 } });
	matrix<int> m2({ { 7, 8 }, { 9, 0 }, { 1, 2 } });
	m1.swap(m2);
	EXPECT_EQ((vector<vector<int>>{{ 7, 8 }, { 9, 0 }, { 1, 2 } }), m1.a);
	EXPECT_EQ((vector<vector<int>>{ { 1, 2, 3 }, { 4, 5, 6 } }), m2.a);
}

TEST(modulo_test, brackets) {
	const matrix<int> m1({ { 1, 2, 3 }, { 4, 5, 6 } });
	EXPECT_EQ((vector<int>{ 1, 2, 3 }), m1[0]);
	EXPECT_EQ(4, m1[1][0]);
	
	matrix<int> m2({ { 1, 2, 3 }, { 4, 5, 6 } });
	vector<int> &row2 = m2[1];
	EXPECT_EQ((vector<int>{ 4, 5, 6 }), row2);
	row2[1] = 7;
	EXPECT_EQ((vector<int>{ 4, 7, 6 }), row2);
	m2[1][0] = 8;
	EXPECT_EQ((vector<vector<int>>{{ 1, 2, 3 }, { 8, 7, 6 } }), m2.a);
}

TEST(matrix_test, operators_comparison) {
	const matrix<int> m1({ { 1, 2, 3 }, { 4, 5, 6 } });
	const matrix<int> m2({ { 7, 8 }, { 9, 0 }, { 1, 2 } });
	EXPECT_EQ(false, m1 == m2);
	EXPECT_EQ(true, m1 != m2);
	EXPECT_EQ(true, m1 < m2);
	EXPECT_EQ(false, m1 > m2);
	EXPECT_EQ(true, m1 <= m2);
	EXPECT_EQ(false, m1 >= m2);
	EXPECT_EQ(false, m2 == m1);
	EXPECT_EQ(true, m2 != m1);
	EXPECT_EQ(false, m2 < m1);
	EXPECT_EQ(true, m2 > m1);
	EXPECT_EQ(false, m2 <= m1);
	EXPECT_EQ(true, m2 >= m1);
	EXPECT_EQ(true, m2 == m2);
	EXPECT_EQ(false, m2 != m2);
	EXPECT_EQ(false, m2 < m2);
	EXPECT_EQ(false, m2 > m2);
	EXPECT_EQ(true, m2 <= m2);
	EXPECT_EQ(true, m2 >= m2);
}

TEST(matrix_test, operators_arithmetic) {
	const matrix<int> m1({ { 1, 2, 3 }, { 4, 5, 6 } });
	const matrix<int> m2({ { 7, 8 }, { 9, 0 }, { 1, 2 } });
	const matrix<int> m3({ { 2, 7 }, { 4, 1 }, { 3, 5 } });
	const matrix<int> m4({ { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9} });

	EXPECT_EQ((matrix<int>{ { 9, 15 }, { 13, 1 }, { 4, 7 } }), m2 + m3);
	EXPECT_EQ((matrix<int>{ { 5, 1 }, { 5, -1 }, { -2, -3 } }), m2 - m3);
	EXPECT_EQ((matrix<int>{ { 28, 14 }, { 79, 44 } }), m1 * m2);
	EXPECT_EQ((matrix<int>{ { 39, 54, 69 }, { 9, 18, 27 }, { 9, 12, 15 } }), m2 * m1);
	EXPECT_EQ((matrix<int>{{ 10, 20, 30 }, { 40, 50, 60 } }), m1 * 10);
	
	matrix <int> mr;
	mr = m2; mr += m3;
	EXPECT_EQ((matrix<int>{ { 9, 15 }, { 13, 1 }, { 4, 7 } }), mr);
	mr = m2; mr -= m3;
	EXPECT_EQ((matrix<int>{ { 5, 1 }, { 5, -1 }, { -2, -3 } }), mr);
	mr = m1; mr *= m2;
	EXPECT_EQ((matrix<int>{ { 28, 14 }, { 79, 44 } }), mr);
	mr = m1; mr *= 10;
	EXPECT_EQ((matrix<int>{{ 10, 20, 30 }, { 40, 50, 60 } }), mr);

	mr = m4; mr += mr;
	EXPECT_EQ((matrix<int>{ {2, 4, 6 }, { 8, 10, 12 }, { 14, 16, 18 } }), mr);
	mr = m4; mr -= mr;
	EXPECT_EQ((matrix<int>{ { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } }), mr);
	mr = m4; mr *= mr;
	EXPECT_EQ((matrix<int>{ { 30, 36, 42 }, { 66, 81, 96 }, { 102, 126, 150 }}), mr);
}

TEST(matrix_test, operators_inverse) {
	const matrix<mod> m1({ { 20, 30, 40 }, { 50, 60, 70 } });
	const matrix<mod> m2({ { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } });
	const matrix<mod> m3({ { 2, 3, 5 }, { 7, 11, 13 }, { 17, 19, 23 } });
	
	EXPECT_EQ(mod(0), m2.det());
	EXPECT_EQ((matrix<mod>{ { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0} }), m2.inverse());
	EXPECT_EQ(mod(-78), m3.det());
	EXPECT_EQ((matrix<mod>{ { 6, 26, -16 }, { 60, -39, 9 }, { -54, 13, 1 }}) / mod(-78), m3.inverse());
	
	EXPECT_EQ((matrix<mod>{ { 2, 3, 4 }, { 5, 6, 7 } }), m1 / mod(10));
	EXPECT_EQ((matrix<mod>{ { -36, -13, 5 }, { 0, -13, -13 }, { 36, -13, -31 }}) / mod(-78), m2 / m3);
	
	matrix <mod> mr;
	mr = m1; mr /= mod(10);
	EXPECT_EQ((matrix<mod>{{ 2, 3, 4 }, { 5, 6, 7 } }), mr);
	mr = m2; mr /= m3;
	EXPECT_EQ((matrix<mod>{{ 36, 13, -5 }, { 0, 13, 13 }, { -36, 13, 31 }}) / mod(78), mr);
}

TEST(matrix_test, power) {
	const matrix<mod> m1({ { 2, 3, 5 }, { 7, 11, 13 }, { 17, 19, 23 } });
	EXPECT_EQ((matrix<mod>{ { 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1} }), m1.pow(0));
	EXPECT_EQ((matrix<mod>{ { 2, 3, 5 }, { 7, 11, 13 }, { 17, 19, 23 } }), m1.pow(1));
	EXPECT_EQ((matrix<mod>{ { 3946, 4920, 6064 }, { 11456, 14278, 17588 }, { 20632, 25700, 31654 } }), m1.pow(3));
	EXPECT_EQ((matrix<mod>{ { -55788, 107120, -48832 }, { 247392, -205764, 66936 }, { -164496, 97240, -22532 } }) / mod(-78 * 78 * 78), m1.pow(-3));
}

TEST(matrix_test, transpose) {
	const matrix<int> m1({ { 1, 2, 3 }, { 4, 5, 6 } });
	EXPECT_EQ((matrix<int>{ { 1, 4 }, { 2, 5 }, { 3, 6 }}), m1.transpose());
}

TEST(matrix_test, identity) {
	EXPECT_EQ((matrix<int>{ { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }}), matrix<int>::identity(3));
}
