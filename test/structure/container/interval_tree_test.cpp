#include "structure/container/interval_tree.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::container;

int inf = numeric_limits<int>::max();

struct atom_min {
	int val;
	bool pending;

	atom_min(int val = inf) : val(val), pending(false) {}

	static void resolve_up(atom_min& parent, const atom_min& left, const atom_min& right) {
		parent.val = min(left.val, right.val);
	}

	static void resolve_down(atom_min& parent, atom_min& left, atom_min& right) {
		if (!parent.pending) return;
		left = right = parent;
		parent.pending = false;
	}

	static std::function<void(atom_min&)> set_functor(int val) {
		return [val](atom_min& t){ t.val = val; t.pending = true; };
	}
};

template<typename T, typename F>
T slow_get(const vector<T>& v, size_t begin, size_t end, F f, const T& id = T()) {
	T t = id;
	for (size_t i = begin; i < end; i++) {
		f(t, t, v[i]);
	}
	return t;
}

template<typename T, typename F>
void verify_all(interval_tree<T>& st, const vector<T>& v, F f, const T& id = T()) {
	for (size_t begin = 0; begin < v.size(); begin++) {
		for (size_t end = begin; end < v.size(); end++) {
			EXPECT_EQ(slow_get(v, begin, end, f, id).val, st.get(begin, end).val) << " unexpected result of get(" << begin << ", " << end << ")";
		}
	}
}

TEST(interval_tree_test, build_int_min) {
	vector<atom_min> v{ 2, -3, 4, 6, 11, 1, 0, -5, 7, -3 };

	interval_tree<atom_min> st1(v.size(), atom_min::resolve_up, atom_min::resolve_down);
	EXPECT_EQ(16, st1.size());
	for (size_t i = 0; i < v.size(); i++) st1.update(i, i + 1, atom_min::set_functor(v[i].val));
	verify_all(st1, v, atom_min::resolve_up);

	interval_tree<atom_min> st2(v.begin(), v.end(), atom_min::resolve_up, atom_min::resolve_down);
	EXPECT_EQ(16, st2.size());
	verify_all(st2, v, atom_min::resolve_up);
}

TEST(interval_tree_test, modify_int_min) {
	vector<atom_min> v{ 2, -3, 4, 6, 11, 1, 0, -5, 7, -3 };
	// make modifications both on verification vector v1
	// and the actual component under test st1;
	// set elements at random indices
	vector<int> beg{ 5, 1, 3, 8, 7, 9, 6, 2, 0, 4 };
	vector<int> end{ 9, 7, 4, 10, 8, 10, 9, 4, 10, 7 };
	vector<atom_min> v1(v.size(), inf);
	interval_tree<atom_min> st1(v.size(), atom_min::resolve_up, atom_min::resolve_down);
	for (size_t i = 0; i < v.size(); i++) {
		int b = beg[i], e = end[i];
		fill(v1.begin() + b, v1.begin() + e, v[b]);
		st1.update(b, e, atom_min::set_functor(v[b].val));
		verify_all(st1, v1, atom_min::resolve_up);
	}
}

TEST(interval_tree_test, modify_rebuild) {
	vector<atom_min> v{ 2, -3, 4, 6, 11, 1, 0, -5, 7, -3 };
	interval_tree<atom_min> st(v.begin(), v.end(), atom_min::resolve_up, atom_min::resolve_down);
	st[3].val = v[3].val = 9;
	st[6].val = v[6].val = 2;
	st[8].val = v[8].val = -7;
	st.rebuild();
	verify_all(st, v, atom_min::resolve_up);
}
