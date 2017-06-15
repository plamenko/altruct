#include "altruct/io/reader.h"

#include "gtest/gtest.h"

#include <fstream>
#include <sstream>
#include <functional>

#ifdef WIN32
    #define snprintf sprintf_s
#endif

using namespace std;
using namespace altruct::io;

namespace {
const char* data = "The quick brown fox jumps over the lazy dog. The End.";
}

void test_reader(reader& rin) {
    EXPECT_EQ('T', rin.read_char());
    char buff[100] = {};
    EXPECT_EQ(12, rin.read(buff, 12));
    EXPECT_EQ("he quick bro", string(buff, buff + 12));
    EXPECT_EQ('w', rin.read_char());
    EXPECT_EQ('n', rin.read_char());
    EXPECT_EQ(' ', rin.read_char());
    EXPECT_EQ(37, rin.read(buff, 100));
    EXPECT_EQ("fox jumps over the lazy dog. The End.", string(buff, buff + 37));
    EXPECT_EQ(0, rin.read(buff, 100));
    EXPECT_EQ(-1, rin.read_char());
}

void test_eof1(reader& rin) {
    char buff[100] = {};
    EXPECT_EQ(53, rin.read(buff, 53));
    EXPECT_EQ("The quick brown fox jumps over the lazy dog. The End.", string(buff, buff + 53));
    EXPECT_TRUE(bool(rin)) << "last read was successful, expecting true";
    EXPECT_EQ(0, rin.read(buff, 1));
    EXPECT_FALSE(bool(rin)) << "last read reached EOF, expecting false";
}

void test_eof2(reader& rin) {
    char buff[100] = {};
    EXPECT_EQ(53, rin.read(buff, 100));
    EXPECT_EQ("The quick brown fox jumps over the lazy dog. The End.", string(buff, buff + 53));
    EXPECT_FALSE(bool(rin)) << "last read reached EOF, expecting false";
}

void test_eof3(reader& rin) {
    char buff[100] = {};
    EXPECT_EQ(52, rin.read(buff, 52));
    EXPECT_EQ("The quick brown fox jumps over the lazy dog. The End", string(buff, buff + 52));
    EXPECT_TRUE(bool(rin)) << "last read was successful, expecting true";
    EXPECT_EQ('.', rin.read_char());
    EXPECT_TRUE(bool(rin)) << "last read was successful, expecting true";
    EXPECT_EQ(-1, rin.read_char());
    EXPECT_FALSE(bool(rin)) << "last read reached EOF, expecting false";
}

template<typename F>
void test_all(F runner) {
    runner(test_reader);
    runner(test_eof1);
    runner(test_eof2);
    runner(test_eof3);
}

TEST(reader_test, file_reader) {
    const char* name = "reader_test_temp_file";

    FILE* file = fopen(name, "w");
    fwrite(data, 1, strlen(data), file);
    fclose(file);

    test_all([&](std::function<void(reader& rin)> test_func) {
        file = fopen(name, "r");
        file_reader rin(file);
        test_func(rin);
        fclose(file);
    });

    remove(name);
}

TEST(reader_test, fstream_reader) {
    const char* name = "reader_test_temp_file";

    ofstream os(name);
    os << data;
    os.close();

    test_all([&](std::function<void(reader& rin)> test_func) {
        ifstream is(name);
        stream_reader rin(is);
        test_func(rin);
        is.close();
    });

    remove(name);
}

TEST(reader_test, sstream_reader) {
    test_all([&](std::function<void(reader& rin)> test_func) {
        istringstream is(data);
        stream_reader rin(is);
        test_func(rin);
    });
}

TEST(reader_test, string_reader) {
    test_all([&](std::function<void(reader& rin)> test_func) {
        string_reader rin(data);
        test_func(rin);
    });
}

TEST(reader_test, buffered_reader) {
    test_all([&](std::function<void(reader& rin)> test_func) {
        string_reader sin(data);
        buffered_reader rin(sin, 1000);
        test_func(rin);
    });

    string_reader sin(data);
    buffered_reader rin(sin, 10);

    // read_char + unread_char
    EXPECT_EQ('T', rin.read_char());
    EXPECT_EQ('h', rin.read_char());
    rin.unread_char();
    EXPECT_EQ('h', rin.read_char());
    rin.unread_char();
    EXPECT_EQ('h', rin.read_char());
    EXPECT_EQ('e', rin.read_char());
    EXPECT_EQ(' ', rin.read_char());

    // reserve + refill
    EXPECT_EQ(6, rin.reserve(0));
    EXPECT_EQ(6, rin.reserve(6));
    // needs refill
    EXPECT_EQ(10, rin.reserve(7));

    // read + refill
    char buff[100] = {};
    EXPECT_EQ(3, rin.read(buff, 3));
    EXPECT_EQ("qui", string(buff, buff + 3));
    // needs refill
    EXPECT_EQ(9, rin.read(buff, 9));
    EXPECT_EQ("ck brown ", string(buff, buff + 9));
    // needs refill, buffer_size characters read
    EXPECT_EQ(10, rin.read(buff, 50));
    EXPECT_EQ("fox jumps ", string(buff, buff + 10));

    // read_char + refill
    EXPECT_EQ(0, rin.reserve(0));
    // needs refill
    EXPECT_EQ('o', rin.read_char());
    EXPECT_EQ('v', rin.read_char());
    EXPECT_EQ(8, rin.reserve(0));

    // data + refill
    EXPECT_STREQ("er the l", rin.data());
    EXPECT_STREQ("er the l", rin.data(8));
    // needs refill
    EXPECT_STREQ("er the laz", rin.data(9));

    // skip
    EXPECT_EQ(4, rin.skip(4));
    EXPECT_EQ(6, rin.reserve(0));
    // only 6 available
    EXPECT_EQ(6, rin.skip(8));
    EXPECT_EQ(0, rin.reserve(0));

    // refill + idempotent
    EXPECT_EQ(10, rin.refill());
    EXPECT_EQ(10, rin.refill());
    EXPECT_STREQ("y dog. The", rin.data());

    // counter + advance + idempotent
    rin.counter() = 7;
    EXPECT_EQ(7, rin.counter());
    rin.advance();
    EXPECT_EQ(0, rin.counter());
    EXPECT_STREQ("The", rin.data());
    rin.advance();
    EXPECT_EQ(0, rin.counter());
    EXPECT_STREQ("The", rin.data());

    // refill + eof
    EXPECT_EQ(8, rin.refill());
    EXPECT_EQ(8, rin.refill());
    EXPECT_STREQ("The End.", rin.data());
}

TEST(reader_test, simple_reader) {
    string_reader in(
        " x yz\n  \n "
        "  \t  w 0  "
        "12   -789 "
        "34   -5678"
        "90   +1234"
        "5678901234"
        "56789_abc "
        " da bud di"
        "n don camp"
        "alon\n come"
        "t alpha \td"
        "ef  0 12.3"
        "45 -1.45e7"
        "+1.45e-81 "
        "-1.45e+81-"
        "skip a1b2c"
        "de!!!     "
        "a 42 3.14 "
        "1234567890"
        "1234567   "
        "1.23456789"
        "012345 lon"
        "g_string_s"
        " 12:45.00 "
        " 11;23;"
    );
    simple_reader rin(in, 10);
    simple_reader_stream sin(rin);

    EXPECT_EQ(' ', rin.read_char());
    EXPECT_EQ('x', rin.read_char());
    EXPECT_EQ(' ', rin.read_char());
    EXPECT_EQ('y', rin.read_char());
    EXPECT_EQ('z', rin.read_char());
    EXPECT_EQ('\n', rin.read_char());
    rin.skip_whitespaces();
    EXPECT_EQ('w', rin.read_char());

    EXPECT_EQ(0, rin.read_int());
    EXPECT_EQ(12, rin.read_int());
    EXPECT_EQ(-789, rin.read_int());
    EXPECT_EQ(34, rin.read_int());
    EXPECT_EQ(-567890, rin.read_int());

    EXPECT_EQ(1234567890123456789LL, rin.read_ll());

    EXPECT_EQ("_abc", rin.read_string());
    EXPECT_EQ("  da bud din don campalon", rin.read_line());
    EXPECT_EQ(" comet alpha ", rin.read_line('\t'));
    EXPECT_EQ("def", rin.read_string());

    EXPECT_NEAR(0, rin.read_float(), 1e-6);
    EXPECT_NEAR(12.345, rin.read_float(), 1e-6);
    EXPECT_NEAR(-1.45e7, rin.read_double(), 1e1 * 1e-14);
    EXPECT_NEAR(+1.45e-81, rin.read_double(), 1e-81 * 1e-14);
    EXPECT_NEAR(-1.45e+81, rin.read_double(), 1e+81 * 1e-14);

    EXPECT_EQ('-', rin.read_char());
    EXPECT_EQ("skip", rin.read_string());

    int x;
    sscanf(rin.data(8), "%x%n", &x, &rin.counter());
    rin.advance();
    EXPECT_EQ(0xa1b2cde, x);

    EXPECT_EQ("!!!", rin.read_string());

    rin.skip_whitespaces();
    char c; sin >> c;
    EXPECT_EQ('a', c);
    int i; sin >> i;
    EXPECT_EQ(42, i);
    float f; sin >> f;
    EXPECT_NEAR(3.14f, f, 1e-6f);
    long long l; sin >> l;
    EXPECT_EQ(12345678901234567LL, l);
    double d; sin >> d;
    EXPECT_NEAR(1.23456789012345, d, 1e-15);
    string s; sin >> s;
    EXPECT_EQ("long_string_s", s);

    int j; sin >> j;
    EXPECT_EQ(12, j);
    EXPECT_EQ(':', rin.read_char());
    int k; sin >> k;
    EXPECT_EQ(45, k);
    EXPECT_EQ('.', rin.read_char());
    int t; sin >> t;
    EXPECT_EQ(0, t);

    int y, z;
    sscanf(rin.data(8), "%d;%d%n", &y, &z, &rin.counter());
    rin.advance();
    EXPECT_EQ(11, y);
    EXPECT_EQ(23, z);
    EXPECT_EQ(';', rin.read_char());
    EXPECT_EQ(-1, rin.read_char());
    EXPECT_EQ("", rin.read_line());
}
