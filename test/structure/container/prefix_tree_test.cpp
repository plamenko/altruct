#include "altruct/structure/container/prefix_tree.h"

#include "altruct/algorithm/collections/collections.h"
#include "altruct/io/iostream_overloads.h"
#include "structure_test_util.h"

#include <functional>
#include <vector>
#include <set>
#include <map>
#include <chrono>

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::container;
using namespace altruct::collections;
using namespace altruct::test_util;

TEST(prefix_tree_test, constructor_and_size) {
    prefix_tree<> pt1;
    EXPECT_EQ(0, pt1.num_words());
    EXPECT_EQ(0, pt1.num_letters());
}

TEST(prefix_tree_test, update_and_query) {
    prefix_tree<> pt1;
    // add & append
    EXPECT_EQ(1, pt1.add("ananas", pt1.ordinal_lower_alpha));     // new word 1
    EXPECT_EQ(1, pt1.num_words());
    EXPECT_EQ(6, pt1.num_letters());
    EXPECT_EQ(2, pt1.add("ana", pt1.ordinal_lower_alpha));        // new word 2, existing prefix of the full length 3
    EXPECT_EQ(2, pt1.num_words());
    EXPECT_EQ(6, pt1.num_letters());
    EXPECT_EQ(1, pt1.add("ananas", pt1.ordinal_lower_alpha));     // existing word 1
    EXPECT_EQ(2, pt1.num_words());
    EXPECT_EQ(6, pt1.num_letters());
    EXPECT_EQ(3, pt1.append(2, "kin", pt1.ordinal_lower_alpha));  // new word 3, ana + kin
    EXPECT_EQ(3, pt1.num_words());
    EXPECT_EQ(9, pt1.num_letters());
    EXPECT_EQ(3, pt1.find("anakin", pt1.ordinal_lower_alpha));    // existing word 3
    EXPECT_EQ(3, pt1.num_words());
    EXPECT_EQ(9, pt1.num_letters());
    EXPECT_EQ(0, pt1.find("anak", pt1.ordinal_lower_alpha));      // prefix, but not a word
    EXPECT_EQ(3, pt1.num_words());
    EXPECT_EQ(9, pt1.num_letters());
    EXPECT_EQ(3, pt1.add("anakin", pt1.ordinal_lower_alpha));     // wxisting word 3
    EXPECT_EQ(3, pt1.num_words());
    EXPECT_EQ(9, pt1.num_letters());
    EXPECT_EQ(4, pt1.add("anakonda", pt1.ordinal_lower_alpha));   // new word 4, existing prefix of length 4
    EXPECT_EQ(4, pt1.num_words());
    EXPECT_EQ(13, pt1.num_letters());
    EXPECT_EQ(5, pt1.append(0, "blah", pt1.ordinal_lower_alpha)); // new word 5, blah
    EXPECT_EQ(5, pt1.num_words());
    EXPECT_EQ(17, pt1.num_letters());

    // get
    EXPECT_EQ("", pt1.get(0, pt1.letter_lower_alpha));
    EXPECT_EQ("ananas", pt1.get(1, pt1.letter_lower_alpha));
    EXPECT_EQ("ana", pt1.get(2, pt1.letter_lower_alpha));
    EXPECT_EQ("anakin", pt1.get(3, pt1.letter_lower_alpha));
    EXPECT_EQ("anakonda", pt1.get(4, pt1.letter_lower_alpha));
    EXPECT_EQ("blah", pt1.get(5, pt1.letter_lower_alpha));
    EXPECT_EQ("BLAH", pt1.get(5, pt1.letter_upper_alpha));

    // for each
    vector<pair<string, int>> vs1;
    pt1.for_each([&](const string& s, int i){ vs1.push_back({ s, i }); }, pt1.letter_lower_alpha);
    EXPECT_EQ((vector<pair<string, int>>{{ "ana", 2 }, { "anakin", 3 }, { "anakonda", 4 }, { "ananas", 1 }, { "blah", 5 }}), vs1);
    vector<pair<string, int>> vs2;
    pt1.for_each(2, [&](const string& s, int i){ vs2.push_back({ s, i }); }, pt1.letter_lower_alpha);
    EXPECT_EQ((vector<pair<string, int>>{{ "ana", 2 }, { "anakin", 3 }, { "anakonda", 4 }, { "ananas", 1 }}), vs2);

    // erase
    pt1.clear();
    EXPECT_EQ(0, pt1.num_words());
    EXPECT_EQ(0, pt1.num_letters());
    EXPECT_EQ(1, pt1.add("ban", pt1.ordinal_lower_alpha));
    EXPECT_EQ(1, pt1.num_words());
    EXPECT_EQ(3, pt1.num_letters());
    EXPECT_EQ(2, pt1.add("banana", pt1.ordinal_lower_alpha));
    EXPECT_EQ(2, pt1.num_words());
    EXPECT_EQ(6, pt1.num_letters());
    EXPECT_EQ(3, pt1.add("bager", pt1.ordinal_lower_alpha));
    EXPECT_EQ(3, pt1.num_words());
    EXPECT_EQ(9, pt1.num_letters());
    EXPECT_EQ(1, pt1.erase("ban", pt1.ordinal_lower_alpha)); // "bager" now becomes a word 1
    EXPECT_EQ(1, pt1.find("bager", pt1.ordinal_lower_alpha));
    EXPECT_EQ(2, pt1.num_words());
    EXPECT_EQ(9, pt1.num_letters());
    EXPECT_EQ(2, pt1.erase("banana", pt1.ordinal_lower_alpha));
    EXPECT_EQ(1, pt1.num_words());
    EXPECT_EQ(5, pt1.num_letters());
    EXPECT_EQ(0, pt1.erase("bag", pt1.ordinal_lower_alpha)); // "bag" doesn't exist
    EXPECT_EQ(1, pt1.num_words());
    EXPECT_EQ(5, pt1.num_letters());
    EXPECT_EQ(1, pt1.erase("bager", pt1.ordinal_lower_alpha));
    EXPECT_EQ(0, pt1.num_words());
    EXPECT_EQ(0, pt1.num_letters());

    // linear tree
    pt1.add("a", pt1.ordinal_lower_alpha);
    pt1.add("aaa", pt1.ordinal_lower_alpha);
    pt1.add("aa", pt1.ordinal_lower_alpha);
    pt1.add("aaaa", pt1.ordinal_lower_alpha);
    pt1.add("aaaaa", pt1.ordinal_lower_alpha);
    EXPECT_EQ(5, pt1.num_words());
    EXPECT_EQ(5, pt1.num_letters());
    vector<pair<string, int>> vs3;
    pt1.for_each([&](const string& s, int i){ vs3.push_back({ s, i }); }, pt1.letter_lower_alpha);
    EXPECT_EQ((vector<pair<string, int>>{{ "a", 1 }, { "aa", 3 }, { "aaa", 2 }, { "aaaa", 4 }, { "aaaaa", 5 }}), vs3);

    // empty word operations are no-op
    EXPECT_EQ(0, pt1.add("", pt1.ordinal_lower_alpha));
    EXPECT_EQ(5, pt1.num_words());
    EXPECT_EQ(5, pt1.num_letters());
    EXPECT_EQ(0, pt1.append(0, "", pt1.ordinal_lower_alpha));
    EXPECT_EQ(5, pt1.num_words());
    EXPECT_EQ(5, pt1.num_letters());
    EXPECT_EQ(3, pt1.append(3, "", pt1.ordinal_lower_alpha));
    EXPECT_EQ(5, pt1.num_words());
    EXPECT_EQ(5, pt1.num_letters());
    EXPECT_EQ(0, pt1.find("", pt1.ordinal_lower_alpha));
    EXPECT_EQ(5, pt1.num_words());
    EXPECT_EQ(5, pt1.num_letters());
    EXPECT_EQ("", pt1.get(0, pt1.letter_lower_alpha));
    EXPECT_EQ(5, pt1.num_words());
    EXPECT_EQ(5, pt1.num_letters());
    EXPECT_EQ(0, pt1.erase("", pt1.ordinal_lower_alpha));
    EXPECT_EQ(5, pt1.num_words());
    EXPECT_EQ(5, pt1.num_letters());
}
