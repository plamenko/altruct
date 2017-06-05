#include "altruct/algorithm/math/fractions.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::math;

typedef fraction<int> frac;

TEST(fractions_test, farey_sequence_prev_inc) {
	vector<frac> vf;
	int n = 8;
	frac qp = frac(-1, 0), q = frac(0, 1), qn;
	for (; q < frac(1, 1); qp = q, q = qn) {
		vf.push_back(q);
		qn = farey_neighbour(n, qp, q);
	}
	vf.push_back(q);
	EXPECT_EQ((vector<frac>{
		{ 0, 1 }, { 1, 8 }, { 1, 7 }, { 1, 6 }, { 1, 5 }, { 1, 4 }, { 2, 7 }, { 1, 3 }, { 3, 8 }, { 2, 5 }, { 3, 7 }, { 1, 2 },
		{ 4, 7 }, { 3, 5 }, { 5, 8 }, { 2, 3 }, { 5, 7 }, { 3, 4 }, { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 8 }, { 1, 1 }}), vf);
}

TEST(fractions_test, farey_sequence_prev_dec) {
	vector<frac> vf;
	int n = 5;
	frac qp = frac(+1, 0), q = frac(1, 1), qn;
	for (; q > frac(0, 1); qp = q, q = qn) {
		vf.push_back(q);
		qn = farey_neighbour(n, qp, q);
	}
	vf.push_back(q);
	EXPECT_EQ((vector<frac>{{1, 1}, { 4, 5 }, { 3, 4 }, { 2, 3 }, { 3, 5 }, { 1, 2 }, { 2, 5 }, { 1, 3 }, { 1, 4 }, { 1, 5 }, { 0, 1 }}), vf);
}

TEST(fractions_test, farey_sequence_inc) {
	vector<frac> vf;
	for (frac q = frac(0, 1); q <= frac(1, 1); q = farey_neighbour<int>(5, { -1, 0 }, q)) {
		vf.push_back(q);
	}
	EXPECT_EQ((vector<frac>{{ 0, 1 }, { 1, 5 }, { 1, 4 }, { 1, 3 }, { 2, 5 }, { 1, 2 }, { 3, 5 }, { 2, 3 }, { 3, 4 }, { 4, 5 }, { 1, 1 }}), vf);
}

TEST(fractions_test, farey_sequence_dec) {
	vector<frac> vf;
	for (frac q = frac(1, 1); q >= frac(0, 1); q = farey_neighbour<int>(5, { +1, 0 }, q)) {
		vf.push_back(q);
	}
	EXPECT_EQ((vector<frac>{{ 1, 1 }, { 4, 5 }, { 3, 4 }, { 2, 3 }, { 3, 5 }, { 1, 2 }, { 2, 5 }, { 1, 3 }, { 1, 4 }, { 1, 5 }, { 0, 1 }}), vf);
}
