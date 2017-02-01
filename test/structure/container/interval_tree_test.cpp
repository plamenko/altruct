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

	static std::function<bool(atom_min&)> set_functor(int val) {
        return [val](atom_min& t){ t.val = val; t.pending = true; return true; };
	}
};

struct atom_sum {
    int size;
    int val;
    int pending;

    atom_sum() : val(0), pending(0), size(0) {}
    atom_sum(int val) : val(val), pending(0), size(1) {}
    atom_sum(int size, int val, int pending) : val(val), pending(pending), size(size) {}

    bool operator == (const atom_sum& rhs) const {
        return size == rhs.size && val == rhs.val && pending == rhs.pending;
    }

    void add(int v) {
        val += v * size;
        pending += v;
    }

    static void resolve_up(atom_sum& parent, const atom_sum& left, const atom_sum& right) {
        parent.size = left.size + right.size;
        parent.val = left.val + right.val;
    }

    static void resolve_down(atom_sum& parent, atom_sum& left, atom_sum& right) {
        left.add(parent.pending);
        right.add(parent.pending);
        parent.pending = 0;
    }

    static std::function<bool(atom_sum&)> add_functor(int val, int max_size = 1000000000) {
        return [val, max_size](atom_sum& t){
            if (t.size > max_size) return false;
            t.add(val); return true;
        };
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

TEST(interval_tree_test, modify_range_rebuild) {
    vector<atom_min> v{ 2, -3, 4, 6, 11, 1, 0, -5, 7, -3 };
    interval_tree<atom_min> st(v.begin(), v.end(), atom_min::resolve_up, atom_min::resolve_down);
    st[6].val = v[6].val = 2;
    st[8].val = v[8].val = -7;
    st.rebuild(6, 8 + 1);
    verify_all(st, v, atom_min::resolve_up);
}

TEST(interval_tree_test, deep_modify) {
    vector<atom_sum> v{ 2, -3, 4, 6, 11, 1, 0, -5, 7, -3, 3, 1, -4, 1, 5, 9, -2, 6, -5, 3, 1 };
    interval_tree<atom_sum> st(v.begin(), v.end(), atom_sum::resolve_up, atom_sum::resolve_down);

    vector<atom_sum> e{
        { 0, 0, 0 },
        { 21, 38, 0 },
        { 16, 35, 0 }, { 5, 3, 0 },
        { 8, 16, 0 }, { 8, 19, 0 }, { 5, 3, 0 }, { 0, 0, 0 },
        { 4, 9, 0 }, { 4, 7, 0 }, { 4, 8, 0 }, { 4, 11, 0 }, { 4, 2, 0 }, { 1, 1, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
        { 2, -1, 0 }, { 2, 10, 0 }, { 2, 12, 0 }, { 2, -5, 0 }, { 2, 4, 0 }, { 2, 4, 0 }, { 2, -3, 0 }, { 2, 14, 0 }, { 2, 4, 0 }, { 2, -2, 0 }, { 1, 1, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
        { 1, 2, 0 }, { 1, -3, 0 }, { 1, 4, 0 }, { 1, 6, 0 }, { 1, 11, 0 }, { 1, 1, 0 }, { 1, 0, 0 }, { 1, -5, 0 }, { 1, 7, 0 }, { 1, -3, 0 }, { 1, 3, 0 }, { 1, 1, 0 }, { 1, -4, 0 }, { 1, 1, 0 }, { 1, 5, 0 }, { 1, 9, 0 }, { 1, -2, 0 }, { 1, 6, 0 }, { 1, -5, 0 }, { 1, 3, 0 }, { 1, 1, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, };
    EXPECT_EQ(e, st.v);
    
    auto st0 = st; st0.update(5, 21, atom_sum::add_functor(8));
    vector<atom_sum> e0{
        { 0, 0, 0 },
        { 21, 166, 0 },
        { 16, 123, 0 }, { 5, 43, 0 },
        { 8, 40, 0 }, { 8, 83, 8 }, { 5, 43, 0 }, { 0, 0, 0 },
        { 4, 9, 0 }, { 4, 31, 0 }, { 4, 8, 0 }, { 4, 11, 0 }, { 4, 34, 8 }, { 1, 9, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
        { 2, -1, 0 }, { 2, 10, 0 }, { 2, 20, 0 }, { 2, 11, 8 }, { 2, 4, 0 }, { 2, 4, 0 }, { 2, -3, 0 }, { 2, 14, 0 }, { 2, 4, 0 }, { 2, -2, 0 }, { 1, 9, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
        { 1, 2, 0 }, { 1, -3, 0 }, { 1, 4, 0 }, { 1, 6, 0 }, { 1, 11, 0 }, { 1, 9, 8 }, { 1, 0, 0 }, { 1, -5, 0 }, { 1, 7, 0 }, { 1, -3, 0 }, { 1, 3, 0 }, { 1, 1, 0 }, { 1, -4, 0 }, { 1, 1, 0 }, { 1, 5, 0 }, { 1, 9, 0 }, { 1, -2, 0 }, { 1, 6, 0 }, { 1, -5, 0 }, { 1, 3, 0 }, { 1, 9, 8 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, };
    EXPECT_EQ(e0, st0.v);
    
    auto st1 = st0; st1.update(3, 17, atom_sum::add_functor(10, 1));
    vector<atom_sum> e1{
        { 0, 0, 0 },
        { 21, 306, 0 },
        { 16, 253, 0 }, { 5, 53, 0 },
        { 8, 90, 0 }, { 8, 163, 0 }, { 5, 53, 0 }, { 0, 0, 0 },
        { 4, 19, 0 }, { 4, 71, 0 }, { 4, 80, 0 }, { 4, 83, 0 }, { 4, 44, 0 }, { 1, 9, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
        { 2, -1, 0 }, { 2, 20, 0 }, { 2, 40, 0 }, { 2, 31, 0 }, { 2, 40, 0 }, { 2, 40, 0 }, { 2, 33, 0 }, { 2, 50, 0 }, { 2, 30, 0 }, { 2, 14, 8 }, { 1, 9, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
        { 1, 2, 0 }, { 1, -3, 0 }, { 1, 4, 0 }, { 1, 16, 10 }, { 1, 21, 10 }, { 1, 19, 18 }, { 1, 18, 18 }, { 1, 13, 18 }, { 1, 25, 18 }, { 1, 15, 18 }, { 1, 21, 18 }, { 1, 19, 18 }, { 1, 14, 18 }, { 1, 19, 18 }, { 1, 23, 18 }, { 1, 27, 18 }, { 1, 16, 18 }, { 1, 14, 8 }, { 1, -5, 0 }, { 1, 3, 0 }, { 1, 9, 8 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, };
    EXPECT_EQ(e1, st1.v);

    auto st2 = st0; st2.update(3, 17, atom_sum::add_functor(10, 2));
    vector<atom_sum> e2{
        { 0, 0, 0 },
        { 21, 306, 0 },
        { 16, 253, 0 }, { 5, 53, 0 },
        { 8, 90, 0 }, { 8, 163, 0 }, { 5, 53, 0 }, { 0, 0, 0 },
        { 4, 19, 0 }, { 4, 71, 0 }, { 4, 80, 0 }, { 4, 83, 0 }, { 4, 44, 0 }, { 1, 9, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
        { 2, -1, 0 }, { 2, 20, 0 }, { 2, 40, 10 }, { 2, 31, 18 }, { 2, 40, 18 }, { 2, 40, 18 }, { 2, 33, 18 }, { 2, 50, 18 }, { 2, 30, 0 }, { 2, 14, 8 }, { 1, 9, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
        { 1, 2, 0 }, { 1, -3, 0 }, { 1, 4, 0 }, { 1, 16, 10 }, { 1, 11, 0 }, { 1, 9, 8 }, { 1, 0, 0 }, { 1, -5, 0 }, { 1, 7, 0 }, { 1, -3, 0 }, { 1, 3, 0 }, { 1, 1, 0 }, { 1, -4, 0 }, { 1, 1, 0 }, { 1, 5, 0 }, { 1, 9, 0 }, { 1, 16, 18 }, { 1, 14, 8 }, { 1, -5, 0 }, { 1, 3, 0 }, { 1, 9, 8 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, };
    EXPECT_EQ(e2, st2.v);

    auto st3 = st0; st3.restore(7, 14);
    vector<atom_sum> e3{
        { 0, 0, 0 },
        { 21, 166, 0 },
        { 16, 123, 0 }, { 5, 43, 0 },
        { 8, 40, 0 }, { 8, 83, 0 }, { 5, 43, 0 }, { 0, 0, 0 },
        { 4, 9, 0 }, { 4, 31, 0 }, { 4, 40, 0 }, { 4, 43, 0 }, { 4, 34, 8 }, { 1, 9, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
        { 2, -1, 0 }, { 2, 10, 0 }, { 2, 20, 0 }, { 2, 11, 0 }, { 2, 20, 0 }, { 2, 20, 0 }, { 2, 13, 0 }, { 2, 30, 8 }, { 2, 4, 0 }, { 2, -2, 0 }, { 1, 9, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
        { 1, 2, 0 }, { 1, -3, 0 }, { 1, 4, 0 }, { 1, 6, 0 }, { 1, 11, 0 }, { 1, 9, 8 }, { 1, 8, 8 }, { 1, 3, 8 }, { 1, 15, 8 }, { 1, 5, 8 }, { 1, 11, 8 }, { 1, 9, 8 }, { 1, 4, 8 }, { 1, 9, 8 }, { 1, 5, 0 }, { 1, 9, 0 }, { 1, -2, 0 }, { 1, 6, 0 }, { 1, -5, 0 }, { 1, 3, 0 }, { 1, 9, 8 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, };
    EXPECT_EQ(e3, st3.v);
}
