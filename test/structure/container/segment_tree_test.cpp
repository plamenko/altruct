#include "altruct/structure/container/segment_tree.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::container;

template<typename T, typename F>
T slow_get(const vector<T>& v, size_t begin, size_t end, F f, const T& id = T()) {
    T t = id;
    for (size_t i = begin; i < end; i++) {
        t = f(t, v[i]);
    }
    return t;
}

template<typename T, typename F>
void verify_all(const segment_tree<T>& st, const vector<T>& v, F f, const T& id = T()) {
    for (size_t begin = 0; begin < v.size(); begin++) {
        EXPECT_EQ(v[begin], st.get(begin)) << " unexpected result of get(" << begin << ")";
        for (size_t end = begin; end < v.size(); end++) {
            EXPECT_EQ(slow_get(v, begin, end, f, id), st.get(begin, end)) << " unexpected result of get(" << begin << ", " << end << ")";
        }
    }
}

TEST(segment_tree_test, build_str_cat) {
    auto concat = [](const string& s1, const string& s2){ return s1 + s2; };
    vector<string> v{ "aaa", "b", "cc", "dddd", "ee", "ff", "g", "hhhhh" };

    segment_tree<string> st1(v.size(), concat);
    EXPECT_EQ(8, st1.size());
    for (size_t i = 0; i < v.size(); i++) st1.set(i, v[i]);
    verify_all(st1, v, concat);

    segment_tree<string> st2(v.begin(), v.end(), concat);
    EXPECT_EQ(8, st2.size());
    verify_all(st2, v, concat);
}

TEST(segment_tree_test, build_int_min) {
    int inf = numeric_limits<int>::max();
    auto min_f = [](int v1, int v2){ return std::min(v1, v2); };
    vector<int> v{ 2, -3, 4, 6, 11, 1, 0, -5, 7, -3 };

    segment_tree<int> st1(v.size(), min_f, inf);
    EXPECT_EQ(16, st1.size());
    for (size_t i = 0; i < v.size(); i++) st1.set(i, v[i]);
    verify_all(st1, v, min_f, inf);

    segment_tree<int> st2(v.begin(), v.end(), min_f, inf);
    EXPECT_EQ(16, st2.size());
    verify_all(st2, v, min_f, inf);
}

TEST(segment_tree_test, modify_int_min) {
    int inf = numeric_limits<int>::max();
    auto min_f = [](int v1, int v2){ return std::min(v1, v2); };
    vector<int> v{ 2, -3, 4, 6, 11, 1, 0, -5, 7, -3 };

    // make modifications both on verification vector v1
    // and the actual component under test st1;
    // set elements at random indices
    vector<int> p{ 5, 1, 3, 8, 7, 9, 6, 2, 0, 4 };
    vector<int> v1(v.size(), inf);
    segment_tree<int> st1(v.size(), min_f, inf);
    for (size_t i = 0; i < v.size(); i++) {
        int j = p[i];
        v1[j] = v[j];
        st1.set(j, v[j]);
        verify_all(st1, v1, min_f, inf);
    }
}

TEST(segment_tree_test, modify_rebuild) {
    int inf = numeric_limits<int>::max();
    auto min_f = [](int v1, int v2){ return std::min(v1, v2); };
    vector<int> v{ 2, -3, 4, 6, 11, 1, 0, -5, 7, -3 };

    segment_tree<int> st(v.begin(), v.end(), min_f, inf);
    st[3] = v[3] = 9;
    st[6] = v[6] = 2;
    st[8] = v[8] = -7;
    st.rebuild();
    verify_all(st, v, min_f, inf);
}

TEST(segment_tree_test, modify_range_rebuild) {
    int inf = numeric_limits<int>::max();
    auto min_f = [](int v1, int v2){ return std::min(v1, v2); };
    vector<int> v{ 2, -3, 4, 6, 11, 1, 0, -5, 7, -3 };

    segment_tree<int> st(v.begin(), v.end(), min_f, inf);
    st[6] = v[6] = 2;
    st[8] = v[8] = -7;
    st.rebuild(6, 8 + 1);
    verify_all(st, v, min_f, inf);
}
