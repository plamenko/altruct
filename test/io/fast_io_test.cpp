#include "io/fast_io.h"

#include "gtest/gtest.h"

#include <cstdlib>

using namespace std;
using namespace altruct::io;

class temp_file {
public:
	FILE* file;
	const char* name;
	temp_file(const char* name, const char* data) : name(name) {
		fopen_s(&file, name, "w");
		fwrite(data, 1, strlen(data), file);
		fclose(file);
		fopen_s(&file, name, "r");
	}
	~temp_file() {
		fclose(file);
		remove(name);
	}
};

TEST(fast_io_test, fast_read) {
	temp_file tmp("fast_io_test_temp_file",
		" x yz\n  \n "
		"  \t  w    "
		"12   -789 "
		"34   -5678"
		"90   +1234"
		"5678901234"
		"56789_abc "
		"defgh 12.3"
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
	);
	auto fin = fast_read(tmp.file, 10);
	
	EXPECT_EQ(' ', fin.read_char());
	EXPECT_EQ('x', fin.read_char());
	EXPECT_EQ(' ', fin.read_char());
	EXPECT_EQ('y', fin.read_char());
	EXPECT_EQ('z', fin.read_char());
	EXPECT_EQ('\n', fin.read_char());
	fin.skip_whitespaces();
	EXPECT_EQ('w', fin.read_char());

	EXPECT_EQ(12, fin.read_int());
	EXPECT_EQ(-789, fin.read_int());
	EXPECT_EQ(34, fin.read_int());
	EXPECT_EQ(-567890, fin.read_int());

	EXPECT_EQ(1234567890123456789LL, fin.read_ll());

	EXPECT_EQ("_abc", fin.read_string());
	EXPECT_EQ("defgh", fin.read_string());

	EXPECT_NEAR(12.345, fin.read_float(), 1e-6);
	EXPECT_NEAR(-1.45e7, fin.read_double(), 1e1 * 1e-14);
	EXPECT_NEAR(+1.45e-81, fin.read_double(), 1e-81 * 1e-14);
	EXPECT_NEAR(-1.45e+81, fin.read_double(), 1e+81 * 1e-14);
	
	EXPECT_EQ('-', fin.read_char());
	EXPECT_EQ("skip", fin.read_string());

	int x;
	fin.reserve(8);
	sscanf_s(fin.ptr, "%x%n", &x, &fin.cnt);
	fin.ptr += fin.cnt;
	EXPECT_EQ(0xa1b2cde, x);

	EXPECT_EQ("!!!", fin.read_string());

	fin.skip_whitespaces();
	char c; fin >> c;
	EXPECT_EQ('a', c);
	int i; fin >> i;
	EXPECT_EQ(42, i);
	float f; fin >> f;
	EXPECT_NEAR(3.14f, f, 1e-6f);
	long long l; fin >> l;
	EXPECT_EQ(12345678901234567LL, l);
	double d; fin >> d;
	EXPECT_NEAR(1.23456789012345, d, 1e-15);
	string s; fin >> s;
	EXPECT_EQ("long_string_s", s);
}
