#include "altruct/io/stream_tokenizer.h"

#include "gtest/gtest.h"

#include <sstream>

using namespace altruct::io::stream_tokenizer;

TEST(stream_tokenizer_test, token_allowed_p) {
    std::stringstream ss("21112101211");
    char _;
    token<allowed_p<'1', '2'>> t1, t2;
    ss >> t1 >> _ >> t2;
    EXPECT_EQ("211121", t1.s);
    EXPECT_EQ("1211", t2.s);
}

TEST(stream_tokenizer_test, token_delimited_p) {
    std::stringstream ss("211121,11;1211");
    char _;
    token<delimited_p<',', ';'>> t1, t2, t3;
    ss >> t1 >> _ >> t2 >> _ >> t3;;
    EXPECT_EQ("211121", t1.s);
    EXPECT_EQ("11", t2.s);
    EXPECT_EQ("1211", t3.s);
}

TEST(stream_tokenizer_test, tokens_int_space) {
    std::stringstream ss("123 45  678   9");
    tokens<int, char> t;
    ss >> t;
    EXPECT_EQ((std::vector<int>{123, 45, 678, 9}), t.v);
}

TEST(stream_tokenizer_test, tokens_int_comma) {
    std::stringstream ss("123,45,678,9");
    tokens<int, char> t;
    ss >> t;
    EXPECT_EQ((std::vector<int>{123, 45, 678, 9}), t.v);
}

TEST(stream_tokenizer_test, tokens_allowed_alphanum) {
    std::stringstream ss("aaa:b-cc ddd;e.ffff");
    tokens_allowed<alphanum_p> t;
    ss >> t;
    EXPECT_EQ((std::vector<std::string>{"aaa", "b", "cc", "ddd", "e", "ffff"}), unbox_tokens(t));
}

TEST(stream_tokenizer_test, tokens_delimited) {
    std::stringstream ss("aaaXbXccYdddXeYffff");
    tokens_delimited<'X', 'Y'> t;
    ss >> t;
    EXPECT_EQ((std::vector<std::string>{"aaa", "b", "cc", "ddd", "e", "ffff"}), unbox_tokens(t));
}

TEST(stream_tokenizer_test, token_alphanum) {
    std::stringstream ss("abc:d-ef");
    char _;
    token_alphanum t1, t2, t3;
    ss >> t1 >> _ >> t2 >> _ >> t3;
    EXPECT_EQ("abc", t1.s);
    EXPECT_EQ("d", t2.s);
    EXPECT_EQ("ef", t3.s);
}

TEST(stream_tokenizer_test, token_binary) {
    std::stringstream ss("10112102011");
    char _;
    token_binary t1, t2, t3;
    ss >> t1 >> _ >> t2 >> _ >> t3;
    EXPECT_EQ("1011", t1.s);
    EXPECT_EQ("10", t2.s);
    EXPECT_EQ("011", t3.s);
}

TEST(stream_tokenizer_test, token_delimited_comma) {
    std::stringstream ss("1234,ab,+-!");
    char _;
    token_delimited_comma t1, t2, t3;
    ss >> t1 >> _ >> t2 >> _ >> t3;
    EXPECT_EQ("1234", t1.s);
    EXPECT_EQ("ab", t2.s);
    EXPECT_EQ("+-!", t3.s);
}

TEST(stream_tokenizer_test, token_delimited_semicolon) {
    std::stringstream ss("1234;ab;+-!");
    char _;
    token_delimited_semicolon t1, t2, t3;
    ss >> t1 >> _ >> t2 >> _ >> t3;
    EXPECT_EQ("1234", t1.s);
    EXPECT_EQ("ab", t2.s);
    EXPECT_EQ("+-!", t3.s);
}

TEST(stream_tokenizer_test, tokens_delimited_space) {
    std::stringstream ss("1234 ab +-!");
    tokens_delimited_space t;
    ss >> t;
    EXPECT_EQ((std::vector<std::string>{"1234", "ab", "+-!"}), unbox_tokens(t));
}

TEST(stream_tokenizer_test, tokens_delimited_comma) {
    std::stringstream ss("1234,ab,+-!");
    tokens_delimited_comma t;
    ss >> t;
    EXPECT_EQ((std::vector<std::string>{"1234", "ab", "+-!"}), unbox_tokens(t));
}

TEST(stream_tokenizer_test, tokens_delimited_semicolon) {
    std::stringstream ss("1234;ab;+-!");
    tokens_delimited_semicolon t;
    ss >> t;
    EXPECT_EQ((std::vector<std::string>{"1234", "ab", "+-!"}), unbox_tokens(t));
}

TEST(stream_tokenizer_test, ints_delimited_space) {
    std::stringstream ss("1234 65 789");
    ints_delimited_space t;
    ss >> t;
    EXPECT_EQ((std::vector<int>{1234, 65, 789}), t.v);
}

TEST(stream_tokenizer_test, int64s_delimited_space) {
    std::stringstream ss("1234000000001 65000000002 789000000003");
    int64s_delimited_space t;
    ss >> t;
    EXPECT_EQ((std::vector<int64_t>{1234000000001, 65000000002, 789000000003}), t.v);
}

TEST(stream_tokenizer_test, unbox_tokens_hex) {
    std::stringstream ss("ff abcd 100");
    tokens_delimited_space t;
    ss >> t;
    auto hex_to_int = [](const std::string& s) { return std::stoi(s, nullptr, 16); };
    auto va = unbox_tokens<int>(t, hex_to_int);
    EXPECT_EQ((std::vector<int>{255, 43981, 256}), va);
}
