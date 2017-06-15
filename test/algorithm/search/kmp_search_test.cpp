#include "altruct/algorithm/search/kmp_search.h"

#include "gtest/gtest.h"

#include <string>

using namespace std;
using namespace altruct::search;

namespace {
template<typename RanIt>
void test(const vector<size_t>& expected_matches, size_t expected_ret, RanIt t, size_t n, RanIt p, size_t m, bool all, const string& message) {
    vector<size_t> actual_matches;
    size_t actual_ret = kmp_search(t, n, p, m, [&](size_t pos) {
        actual_matches.push_back(pos);
        return all;
    });
    EXPECT_EQ(expected_matches, actual_matches) << message;
    EXPECT_EQ(expected_ret, actual_ret) << message;
}

void test(const vector<size_t>& expected_matches, size_t expected_ret, const string& t, const string& p) {
    string message = "'" + t + "', '" + p + "'";
    vector<size_t> expected_first_match;
    size_t expected_first_ret = expected_ret;
    if (expected_matches.size() > 0) {
        expected_first_match.push_back(expected_matches[0]);
        expected_first_ret = expected_matches[0];
    }
    test(expected_matches, expected_ret, t.begin(), t.size(), p.begin(), p.size(), true, message + ", iterator, all");
    test(expected_first_match, expected_first_ret, t.begin(), t.size(), p.begin(), p.size(), false, message + ", iterator, first");
    test(expected_matches, expected_ret, t.c_str(), t.size(), p.c_str(), p.size(), true, message + ", pointer, all");
    test(expected_first_match, expected_first_ret, t.c_str(), t.size(), p.c_str(), p.size(), false, message + ", pointer, first");
}
}

TEST(kmp_search_test, pointer) {
    // both empty
    test((vector<size_t>{}), 0, "", "");

    // empty pattern
    test((vector<size_t>{}), 0, "a", "");
    test((vector<size_t>{}), 0, "aa", "");
    test((vector<size_t>{}), 0, "abc", "");

    // empty text
    test((vector<size_t>{}), 0, "", "a");
    test((vector<size_t>{}), 0, "", "aa");
    test((vector<size_t>{}), 0, "", "abc");

    // pattern longer than text
    test((vector<size_t>{}), 1, "a", "aa");
    test((vector<size_t>{}), 1, "a", "aaa");
    test((vector<size_t>{}), 2, "aa", "aaa");
    test((vector<size_t>{}), 2, "aa", "aaaa");

    // single letter pattern
    test((vector<size_t>{0}), 0, "a", "a");
    test((vector<size_t>{0,1}), 1, "aa", "a");
    test((vector<size_t>{0,1,2}), 2, "aaa", "a");

    // repeated pattern
    test((vector<size_t>{0}), 0, "ab", "ab");
    test((vector<size_t>{0, 2}), 2, "abab", "ab");
    test((vector<size_t>{0, 2, 4}), 4, "ababab", "ab");

    // repeating pattern
    test((vector<size_t>{0}), 0, "abcabc", "abcabc");
    test((vector<size_t>{0, 3}), 3, "abcabcabc", "abcabc");
    test((vector<size_t>{0, 3, 6}), 6, "abcabcabcabc", "abcabc");

    // complex pattern
    test((vector<size_t>{15, 24}), 24, "ABC ABCDAB ABCDABCDABDE ABCDABD", "ABCDABD");
    test((vector<size_t>{}), 31, "ABC ABCDAB ABCDABCDABDE ABCDABD", "ABCDABDX");
}
