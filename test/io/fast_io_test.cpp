#include "io/fast_io.h"

#include "gtest/gtest.h"

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <functional>

using namespace std;
using namespace altruct::io;

string read_file(const char* name) {
	ifstream is(name);
	stringstream ss;
	ss << is.rdbuf();
	return ss.str();
}

class temp_read_file {
public:
	FILE* file;
	const char* name;
	temp_read_file(const char* name, const char* data) : name(name) {
		fopen_s(&file, name, "w");
		fwrite(data, 1, strlen(data), file);
		fclose(file);
		fopen_s(&file, name, "r");
	}
	~temp_read_file() {
		fclose(file);
		remove(name);
	}
};

class temp_write_file {
public:
	FILE* file;
	const char* name;
	temp_write_file(const char* name) : name(name) {
		fopen_s(&file, name, "w");
	}
	~temp_write_file() {
		fclose(file);
		remove(name);
	}
	void close() {
		fclose(file);
	}
};

TEST(fast_io_test, fast_read) {
	temp_read_file tmp("fast_io_test_temp_file",
		" x yz\n  \n "
		"  \t  w 0  "
		"12   -789 "
		"34   -5678"
		"90   +1234"
		"5678901234"
		"56789_abc "
		"def 0 12.3"
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

	EXPECT_EQ(0, fin.read_int());
	EXPECT_EQ(12, fin.read_int());
	EXPECT_EQ(-789, fin.read_int());
	EXPECT_EQ(34, fin.read_int());
	EXPECT_EQ(-567890, fin.read_int());

	EXPECT_EQ(1234567890123456789LL, fin.read_ll());

	EXPECT_EQ("_abc", fin.read_string());
	EXPECT_EQ("def", fin.read_string());

	EXPECT_NEAR(0, fin.read_float(), 1e-6);
	EXPECT_NEAR(12.345, fin.read_float(), 1e-6);
	EXPECT_NEAR(-1.45e7, fin.read_double(), 1e1 * 1e-14);
	EXPECT_NEAR(+1.45e-81, fin.read_double(), 1e-81 * 1e-14);
	EXPECT_NEAR(-1.45e+81, fin.read_double(), 1e+81 * 1e-14);
	
	EXPECT_EQ('-', fin.read_char());
	EXPECT_EQ("skip", fin.read_string());

	int x;
	sscanf_s(fin.ptr, "%x%n", &x, fin.reserve_cnt(8));
	fin.advance();
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

string do_write(function<void(fast_write&)> writer) {
	temp_write_file tmp("fast_io_test_temp_file");
	auto fout = fast_write(tmp.file, 10);
	writer(fout);
	fout.flush();
	tmp.close();
	return read_file(tmp.name);
}

string keep(string s, int num) {
	for (char& c : s) {
		if ('0' <= c && c <= '9' && --num < 0) {
			c = '#';
		}
	}
	return s;
}

TEST(fast_io_test, fast_write) {
	EXPECT_EQ("x", do_write([](fast_write& fout){ fout.write_char('x'); }));
	EXPECT_EQ("test", do_write([](fast_write& fout){ fout.write_string("test"); }));
	EXPECT_EQ("string", do_write([](fast_write& fout){ fout.write_string(string("string")); }));
	
	EXPECT_EQ("0", do_write([](fast_write& fout){ fout.write_int(0); }));
	EXPECT_EQ("42", do_write([](fast_write& fout){ fout.write_int(42); }));
	EXPECT_EQ("-42", do_write([](fast_write& fout){ fout.write_int(-42); }));
	EXPECT_EQ("12345678901234567", do_write([](fast_write& fout){ fout.write_ll(12345678901234567LL); }));
	EXPECT_EQ("-12345678901234567", do_write([](fast_write& fout){ fout.write_ll(-12345678901234567LL); }));

	EXPECT_EQ("0.000000", do_write([](fast_write& fout){ fout.write_float(0.f); }));
	EXPECT_EQ("0.00001234567", do_write([](fast_write& fout){ fout.write_float(0.00001234567f, 11); }));
	EXPECT_EQ("-0.00001234567", do_write([](fast_write& fout){ fout.write_float(-0.00001234567f, 11); }));
	EXPECT_EQ("0.230000", do_write([](fast_write& fout){ fout.write_float(0.230000f); }));
	EXPECT_EQ("-0.230000", do_write([](fast_write& fout){ fout.write_float(-0.230000f); }));
	EXPECT_EQ("1.230000", do_write([](fast_write& fout){ fout.write_float(1.230000f); }));
	EXPECT_EQ("-1.230000", do_write([](fast_write& fout){ fout.write_float(-1.230000f); }));
	EXPECT_EQ("123.4500", do_write([](fast_write& fout){ fout.write_float(123.45f, 4); }));
	EXPECT_EQ("-123.4500", do_write([](fast_write& fout){ fout.write_float(-123.45f, 4); }));
	EXPECT_EQ("123.6", do_write([](fast_write& fout){ fout.write_float(123.6499f, 1); }));
	EXPECT_EQ("-123.6", do_write([](fast_write& fout){ fout.write_float(-123.6499f, 1); }));
	EXPECT_EQ("123.7", do_write([](fast_write& fout){ fout.write_float(123.6500f, 1); }));
	EXPECT_EQ("-123.7", do_write([](fast_write& fout){ fout.write_float(-123.6500f, 1); }));
	
	EXPECT_EQ("1234500################.######", keep(do_write([](fast_write& fout){ fout.write_float(123.45e20f); }), 7));
	EXPECT_EQ("-1234500################.######", keep(do_write([](fast_write& fout){ fout.write_float(-123.45e20f); }), 7));
	EXPECT_EQ("12345678901234###################.######", keep(do_write([](fast_write& fout){ fout.write_double(123.456789012345e30); }), 14));
	EXPECT_EQ("-12345678901234###################.######", keep(do_write([](fast_write& fout){ fout.write_double(-123.456789012345e30); }), 14));
	
	EXPECT_EQ("0.0000000000000e+000", do_write([](fast_write& fout){ fout.write_double(0.0, 13, true); }));
	EXPECT_EQ("1.2345678901234e+000", do_write([](fast_write& fout){ fout.write_double(1.2345678901234, 13, true); }));
	EXPECT_EQ("-1.2345678901234e+000", do_write([](fast_write& fout){ fout.write_double(-1.2345678901234, 13, true); }));
	EXPECT_EQ("1.2345678901234e-020", do_write([](fast_write& fout){ fout.write_double(1.2345678901234e-20, 13, true); }));
	EXPECT_EQ("-1.2345678901234e-020", do_write([](fast_write& fout){ fout.write_double(-1.2345678901234e-20, 13, true); }));
	EXPECT_EQ("1.2345678901234e+020", do_write([](fast_write& fout){ fout.write_double(1.2345678901234e+20, 13, true); }));
	EXPECT_EQ("-1.2345678901234e+020", do_write([](fast_write& fout){ fout.write_double(-1.2345678901234e+20, 13, true); }));
	EXPECT_EQ("1.2345678901234e-300", do_write([](fast_write& fout){ fout.write_double(1.2345678901234e-300, 13, true); }));
	EXPECT_EQ("-1.2345678901234e-300", do_write([](fast_write& fout){ fout.write_double(-1.2345678901234e-300, 13, true); }));

	EXPECT_EQ("1.236e+002", do_write([](fast_write& fout){ fout.write_float(123.6499f, 3, true); }));
	EXPECT_EQ("-1.236e+002", do_write([](fast_write& fout){ fout.write_float(-123.6499f, 3, true); }));
	EXPECT_EQ("1.237e+002", do_write([](fast_write& fout){ fout.write_float(123.6500f, 3, true); }));
	EXPECT_EQ("-1.237e+002", do_write([](fast_write& fout){ fout.write_float(-123.6500f, 3, true); }));

	EXPECT_EQ("yz", do_write([](fast_write& fout){ fout << 'y' << 'z'; }));
	EXPECT_EQ("concat", do_write([](fast_write& fout){ fout << "con" << "cat"; }));
	EXPECT_EQ("421", do_write([](fast_write& fout){ fout << 42 << 1; }));
	EXPECT_EQ("12345678901234567 aaa", do_write([](fast_write& fout){ fout << 12345678901234567 << " aaa"; }));
	EXPECT_EQ("2.718282 3.141593", do_write([](fast_write& fout){ fout << 2.7182818f << " " << 3.1415926535897932384626433832795; }));

	EXPECT_EQ("random_prefix_xyz 12345678 suffix", do_write([](fast_write& fout){ 
		fout << "random_prefix_xyz ";
		sprintf_s(fout.ptr, fout.reserve(9), "%x", 0x12345678);
		fout.advance();
		fout << " suffix";
	}));
}
