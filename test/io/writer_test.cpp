#include "io/writer.h"

#include "gtest/gtest.h"

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <functional>

using namespace std;
using namespace altruct::io;

namespace {
string read_file(const char* name) {
	ifstream is(name);
	stringstream ss;
	ss << is.rdbuf();
	return ss.str();
}

string keep(string s, int num) {
	for (char& c : s) {
		if ('0' <= c && c <= '9' && --num < 0) {
			c = '#';
		}
	}
	return s;
}

string do_write(function<void(simple_writer&, simple_writer_stream&)> writer) {
	ostringstream os;
	stream_writer out(os);
	simple_writer wout(out, 10);
	simple_writer_stream sout(wout);
	writer(wout, sout);
	wout.flush();
	return os.str();
}

const char* data = "The quick brown fox jumps over the lazy dog. The End.";
}

void do_write(writer& wout) {
	wout.write_char('T');
	string data1 = "he quick brown fox jumps over the lazy dog.";
	wout.write(data1.c_str(), data1.size());
	wout.write_char(' ');
	string data2 = "The End.";
	wout.write(data2.c_str(), data2.size());
}

TEST(writer_test, file_writer) {
	const char* name = "writer_test_temp_file";

	FILE* file = nullptr;
	fopen_s(&file, name, "w");
	{
		file_writer wout(file);
		do_write(wout);
	}
	fclose(file);

	EXPECT_EQ(data, read_file(name));
	
	remove(name);
}

TEST(writer_test, fstream_writer) {
	const char* name = "writer_test_temp_file";

	ofstream os(name);
	{
		stream_writer wout(os);
		do_write(wout);
	}
	os.close();

	EXPECT_EQ(data, read_file(name));

	remove(name);
}

TEST(writer_test, sstream_writer) {
	ostringstream os;
	stream_writer wout(os);
	do_write(wout);

	EXPECT_EQ(data, os.str());
}

TEST(writer_test, string_writer) {
	vector<char> buff(100);
	string_writer wout(buff.data(), buff.size());
	do_write(wout);

	EXPECT_STREQ(data, buff.data());
}

TEST(writer_test, buffered_writer) {
	vector<char> buff(100);
	string_writer out(buff.data(), buff.size());
	{
		buffered_writer wout(out, 10);

		// buffered write
		wout.write_char('T');
		string data1 = "he quick ";
		wout.write(data1.c_str(), data1.size());
		EXPECT_STREQ("", buff.data());

		// auto flush
		wout.write_char('b');
		EXPECT_STREQ("The quick ", buff.data());
		string data2 = "rown ";
		wout.write(data2.c_str(), data2.size());
		EXPECT_STREQ("The quick ", buff.data());

		// manual flush
		wout.flush();
		EXPECT_STREQ("The quick brown ", buff.data());
		// flush idempotent
		wout.flush();
		EXPECT_STREQ("The quick brown ", buff.data());

		// available
		EXPECT_EQ(10, wout.available());
		string data3 = "fox ";
		wout.write(data3.c_str(), data3.size());
		EXPECT_STREQ("The quick brown ", buff.data());
		EXPECT_EQ(6, wout.available());

		// reserve
		wout.reserve(5);
		EXPECT_EQ(6, wout.available());
		wout.reserve(6);
		EXPECT_EQ(6, wout.available());
		// reserve auto flush
		wout.reserve(7);
		EXPECT_EQ(10, wout.available());
		EXPECT_STREQ("The quick brown fox ", buff.data());

		// write multiple buffers at once
		EXPECT_EQ(10, wout.available());
		string data4 = "jumps over the lazy dog. ";
		wout.write(data4.c_str(), data4.size());
		EXPECT_STREQ("The quick brown fox jumps over the lazy ", buff.data());
		EXPECT_EQ(5, wout.available());

		// data + advance
		wout.reserve(4);
		sprintf_s(wout.data(), wout.available(), "The ");
		wout.advance();
		EXPECT_STREQ("The quick brown fox jumps over the lazy ", buff.data());
		EXPECT_EQ(1, wout.available());
		// reserve flush
		wout.reserve(4);
		sprintf_s(wout.data(), wout.available(), "End.");
		wout.advance();
		EXPECT_STREQ("The quick brown fox jumps over the lazy dog. The ", buff.data());
		EXPECT_EQ(6, wout.available());
	}
	// destructor
	EXPECT_STREQ("The quick brown fox jumps over the lazy dog. The End.", buff.data());
}

TEST(writer_test, simple_writer) {
	EXPECT_EQ("x", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_char('x'); }));
	EXPECT_EQ("test", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_string("test"); }));
	EXPECT_EQ("string", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_string(string("string")); }));
	
	EXPECT_EQ("0", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_int(0); }));
	EXPECT_EQ("42", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_int(42); }));
	EXPECT_EQ("-42", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_int(-42); }));
	EXPECT_EQ("12345678901234567", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_ll(12345678901234567LL); }));
	EXPECT_EQ("-12345678901234567", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_ll(-12345678901234567LL); }));

	EXPECT_EQ("0.000000", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_float(0.f); }));
	EXPECT_EQ("0.00001234567", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_float(0.00001234567f, 11); }));
	EXPECT_EQ("-0.00001234567", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_float(-0.00001234567f, 11); }));
	EXPECT_EQ("0.230000", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_float(0.230000f); }));
	EXPECT_EQ("-0.230000", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_float(-0.230000f); }));
	EXPECT_EQ("1.230000", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_float(1.230000f); }));
	EXPECT_EQ("-1.230000", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_float(-1.230000f); }));
	EXPECT_EQ("123.4500", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_float(123.45f, 4); }));
	EXPECT_EQ("-123.4500", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_float(-123.45f, 4); }));
	EXPECT_EQ("123.6", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_float(123.6499f, 1); }));
	EXPECT_EQ("-123.6", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_float(-123.6499f, 1); }));
	EXPECT_EQ("123.7", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_float(123.6500f, 1); }));
	EXPECT_EQ("-123.7", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_float(-123.6500f, 1); }));
	
	EXPECT_EQ("1234500################.######", keep(do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_float(123.45e20f); }), 7));
	EXPECT_EQ("-1234500################.######", keep(do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_float(-123.45e20f); }), 7));
	EXPECT_EQ("12345678901234###################.######", keep(do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_double(123.456789012345e30); }), 14));
	EXPECT_EQ("-12345678901234###################.######", keep(do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_double(-123.456789012345e30); }), 14));
	
	EXPECT_EQ("0.0000000000000e+000", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_double(0.0, 13, true); }));
	EXPECT_EQ("1.2345678901234e+000", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_double(1.2345678901234, 13, true); }));
	EXPECT_EQ("-1.2345678901234e+000", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_double(-1.2345678901234, 13, true); }));
	EXPECT_EQ("1.2345678901234e-020", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_double(1.2345678901234e-20, 13, true); }));
	EXPECT_EQ("-1.2345678901234e-020", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_double(-1.2345678901234e-20, 13, true); }));
	EXPECT_EQ("1.2345678901234e+020", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_double(1.2345678901234e+20, 13, true); }));
	EXPECT_EQ("-1.2345678901234e+020", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_double(-1.2345678901234e+20, 13, true); }));
	EXPECT_EQ("1.2345678901234e-300", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_double(1.2345678901234e-300, 13, true); }));
	EXPECT_EQ("-1.2345678901234e-300", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_double(-1.2345678901234e-300, 13, true); }));

	EXPECT_EQ("1.236e+002", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_float(123.6499f, 3, true); }));
	EXPECT_EQ("-1.236e+002", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_float(-123.6499f, 3, true); }));
	EXPECT_EQ("1.237e+002", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_float(123.6500f, 3, true); }));
	EXPECT_EQ("-1.237e+002", do_write([](simple_writer& wout, simple_writer_stream& sout){ wout.write_float(-123.6500f, 3, true); }));

	EXPECT_EQ("yz", do_write([](simple_writer& wout, simple_writer_stream& sout){ sout << 'y' << 'z'; }));
	EXPECT_EQ("concat", do_write([](simple_writer& wout, simple_writer_stream& sout){ sout << "con" << "cat"; }));
	EXPECT_EQ("421", do_write([](simple_writer& wout, simple_writer_stream& sout){ sout << 42 << 1; }));
	EXPECT_EQ("12345678901234567 aaa", do_write([](simple_writer& wout, simple_writer_stream& sout){ sout << 12345678901234567 << " aaa"; }));
	EXPECT_EQ("2.718282 3.141593", do_write([](simple_writer& wout, simple_writer_stream& sout){ sout << 2.7182818f << " " << 3.1415926535897932384626433832795; }));

	EXPECT_EQ("random_prefix_xyz 12345678 suffix", do_write([](simple_writer& wout, simple_writer_stream& sout){
		sout << "random_prefix_xyz ";
		wout.reserve(9);
		sprintf_s(wout.data(), wout.available(), "%x", 0x12345678);
		wout.advance();
		sout << " suffix";
	}));
}
