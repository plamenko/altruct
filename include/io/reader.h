#pragma once

#include <cmath>
#include <cstdio>
#include <string>
#include <iostream>
#include <algorithm>

namespace altruct {
namespace io {

/**
 * A character reader interface.
 */
class reader {
public:
	// Reads and returns the next character.
	// Returns -1 if no characters are available.
	virtual int read_char() = 0;
	
	// Tries to read `count` characters into `buffer`.
	// Returns the number of characters read.
	virtual size_t read(char* buffer, size_t count) = 0;
	
	// Returns false if an error flags is set, and true otherwise.
	// Error flag is usually set if EOF occured during a read operation.
	// Note that the value returned by this function depends on the last
	// operation performed on the stream (and not on the next).
	virtual operator bool() const = 0;

	virtual ~reader() {}
	reader() {}
	reader(const reader&) = delete;
	reader& operator = (const reader&) = delete;
};

/**
 * A character reader implementation backed by a file.
 */
class file_reader : public reader {
	FILE* in;
	bool read_failed;
	
public:
	// Constructs a new `file_reader` on top of `in` file.
	// `in` pointer is owned by the caller!
	file_reader(FILE* in) : in(in), read_failed(false) {}

	int read_char() override {
		int res = fgetc(in);
		if (res == -1) read_failed = true;
		return res;
	}

	size_t read(char* buffer, size_t count) override {
		size_t res = fread(buffer, 1, count, in);
		if (res < count) read_failed = true;
		return res;
	}

	operator bool() const override {
		return !read_failed;
	}
};

/**
 * A character reader implementation backed by an input stream.
 */
class stream_reader : public reader {
	std::istream& in;
	
public:
	// Constructs a new `stream_reader` on top of `in` stream.
	// The caller needs to keep the underlying stream alive!
	stream_reader(std::istream& in) : in(in) {}

	int read_char() override {
		return in.get();
	}

	size_t read(char* buffer, size_t count) override {
		in.read(buffer, count);
		return (size_t)in.gcount();
	}
	
	operator bool() const override {
		return bool(in);
	}
};

/**
 * A character reader implementation backed by a string.
 */
class string_reader : public reader {
	const char* in;
	size_t len;
	size_t pos;
	bool read_failed;

public:
	// Constructs a new `string_reader` on top of `in` string.
	// `in` pointer is owned by the caller!
	string_reader(const char* in) : string_reader(in, strlen(in)) {}
	string_reader(const char* in, size_t len) : in(in), len(len), pos(0), read_failed(false){}

	int read_char() override {
		if (pos >= len) {
			read_failed = true;
			return -1;
		}
		return in[pos++];
	}

	size_t read(char* buffer, size_t count) override {
		size_t available = len - pos;
		if (available < count) read_failed = true;
		count = std::min(count, available);
		std::copy(in + pos, in + pos + count, buffer);
		pos += count;
		return count;
	}
	
	operator bool() const override {
		return !read_failed;
	}
};

/**
 * A class that buffers the underlying reader.
 *
 * Example of formatted input:
 *   sscanf_s(reader.data(100), "%x %d %f %n", &v1, &v2, &v3, &reader.counter());
 *   reader.advance();
 */
class buffered_reader : public reader {
	reader& in;
	char* buff;
	char* ptr;
	char* end;
	size_t capacity;
	int cnt; // used for formatted input
	bool read_failed;

public:
	// Constructs a new `buffered_reader` on top of `in` reader.
	// The caller needs to keep the underlying reader alive!
	buffered_reader(reader& in, size_t buffer_size = 1 << 20) : in(in), cnt(0), read_failed(false) {
		capacity = buffer_size;
		buff = new char[capacity + 1];
		buff[capacity] = 0;
		ptr = end = buff;
	}

	~buffered_reader() {
		delete[] buff;
	}

	operator bool() const {
		return !read_failed;
	}

	// Refills the buffer.
	// Returns the number of available characters.
	size_t refill() {
		// move the available characters to the beginning of the buffer
		size_t available = end - ptr;
		std::copy(ptr, end, buff);
		ptr = buff;
		end = buff + available;
		// refill
		while (available < capacity) {
			size_t len = in.read(end, capacity - available);
			if (len == 0) break;
			available += len;
			end += len;
		}
		*end = 0;
		return available;
	}

	// Tries to reserve at least `count` characters.
	// Returns the number of available characters.
	size_t reserve(size_t count = 1024) {
		size_t available = end - ptr;
		return (available >= count) ? available : refill();
	}

	// Returns a pointer to the raw character buffer,
	// provided that there are at least `count` characters
	// available, or `nullptr` otherwise.
	const char* data(size_t count = 0) {
		size_t available = reserve(count);
		if (available < count) {
			read_failed = true;
			return nullptr;
		}
		return ptr;
	}

	// Returns the reference to `cnt` variable.
	// This is convenient for `sscanf` usages.
	int& counter() {
		return cnt;
	}

	// Skips `cnt` characters and resets `cnt` to zero.
	// This is convenient for `sscanf` usages.
	void advance() {
		skip(cnt);
		cnt = 0;
	}

	// Tries to skip `count` characters.
	// Returns the number of characters consumed.
	size_t skip(size_t count) {
		size_t available = end - ptr;
		count = std::min(count, available);
		ptr += count;
		return count;
	}

	// Unreads the last read character.
	// Allowed to be called only immediately
	// after a sucessfull `read_char()`.
	void unread_char() {
		--ptr;
	}

	int read_char() override {
		size_t available = reserve(1);
		if (available < 1) {
			read_failed = true;
			return -1;
		}
		return *ptr++;
	}

	size_t read(char* buffer, size_t count) override {
		size_t available = reserve(count);
		if (available < count) read_failed = true;
		count = std::min(count, available);
		std::copy(ptr, ptr + count, buffer);
		ptr += count;
		return count;
	}
};

/**
 * A class that buffers the underlying reader and provides basic input facilities.
 *
 * Example of simple input:
 *   int x = reader.read_int();
 *   string s = reader.read_string();
 *
 * Example of formatted input:
 *   sscanf_s(reader.data(100), "%x %d %f %n", &v1, &v2, &v3, &reader.counter());
 *   reader.advance();
 */
class simple_reader : public buffered_reader {
public:
	// Constructs a new `simple_reader` on top of `in` reader.
	// The caller needs to keep the underlying reader alive!
	simple_reader(reader& in, size_t buffer_size = 1 << 20) : buffered_reader(in, buffer_size) {
	}

	// Reads a string of characters until a delimiter.
	// Delimiter is consumed, but not returned.
	std::string read_line(char delimiter = '\n') {
		std::string s;
		int c = read_char();
		while (c != delimiter && c != -1) {
			s.push_back((char)c);
			c = read_char();
		}
		return s;
	}
	
	// Reads a string of characters until a whitespace.
	std::string read_string() {
		skip_whitespaces();
		std::string s;
		int c = read_char();
		while (c > 0x20) {
			s.push_back((char)c);
			c = read_char();
		}
		if (c != -1) unread_char();
		return s;
	}

	// Skips whitespaces.
	void skip_whitespaces() {
		int c = read_char();
		while (0x00 <= c && c <= 0x20) { // \t \n ' ' ...
			c = read_char();
		}
		if (c != -1) unread_char();
	}

	// Reads a sign.
	// Character is consumed only if '+' or '-'.
	int read_sign() {
		int c = read_char();
		if (c == '-') return -1;
		if (c == '+') return +1;
		if (c != -1) unread_char();
		return +1;
	}

	// Reads decimal digits.
	template<typename T>
	T read_digits(T d = 0, int* cnt = nullptr) {
		int l = 0;
		int c = read_char();
		while ('0' <= c && c <= '9') {
			d = d * 10 + (c - '0');
			l++;
			c = read_char();
		}
		if (c != -1) unread_char();
		if (cnt) *cnt = l;
		return d;
	}

	// Reads a single signed decimal integer of type I.
	template<typename I>
	I read_integral_number() {
		skip_whitespaces();
		int s = read_sign();
		I d = read_digits<I>();
		return (s < 0) ? -d : d;
	}

	// Reads a single signed decimal floating-point number of type F.
	template<typename F>
	F read_floating_point_number() {
		skip_whitespaces();
		int s = read_sign();
		F f = read_digits<F>();
		int e = 0;
		int c = read_char();
		if (c == '.') {
			int cnt = 0;
			f = read_digits<F>(f, &cnt);
			e -= cnt;
			c = read_char();
		}
		if (c == 'e' || c == 'E') {
			int z = read_sign();
			e += z * read_digits<int>();
			c = read_char();
		}
		if (c != -1) unread_char();
		f *= pow(F(10), e);
		return (s < 0) ? -f : f;
	}

	// Reads a single int.
	int read_int() {
		return read_integral_number<int>();
	}

	// Reads a single long.
	long read_l() {
		return read_integral_number<long>();
	}

	// Reads a single long long.
	long long read_ll() {
		return read_integral_number<long long>();
	}

	// Reads a single float.
	float read_float() {
		return read_floating_point_number<float>();
	}

	// Reads a single double.
	double read_double() {
		return read_floating_point_number<double>();
	}
};

/**
 * A class that provides streaming interface on top of `simple_reader`.
 *
 * Example of streamed input:
 *   in >> x >> s;
 */
class simple_reader_stream {
	simple_reader& in;

public:
	// Constructs a new `simple_reader_stream` on top of `in` reader.
	// The caller needs to keep the underlying reader alive!
	simple_reader_stream(simple_reader& in) : in(in) {}

	operator bool() const { return bool(in); }

	std::string read_line(char delimiter = '\n') { return in.read_line(delimiter); }
	
	simple_reader_stream& operator >> (char& v) { v = in.read_char(); return *this; }
	simple_reader_stream& operator >> (std::string& s) { s = in.read_string(); return *this; }
	simple_reader_stream& operator >> (int& v) { v = in.read_int(); return *this; }
	simple_reader_stream& operator >> (long& v) { v = in.read_l(); return *this; }
	simple_reader_stream& operator >> (long long& v) { v = in.read_ll(); return *this; }
	simple_reader_stream& operator >> (float& v) { v = in.read_float(); return *this; }
	simple_reader_stream& operator >> (double& v) { v = in.read_double(); return *this; }
};

template<typename V = void> simple_reader& simple_in_reader() {
	static file_reader fin(stdin);
	static simple_reader sin(fin);
	return sin;
}
template<typename V = void> simple_reader_stream& simple_in_stream() {
	static simple_reader_stream sin(simple_in_reader());
	return sin;
}

}
}
