#pragma once

#include <cstdio>
#include <vector>
#include <string>
#include <iostream>

namespace altruct {
namespace io {

class fast_read {
public:
	FILE* fin;
	char* buff;
	char* ptr;
	char* end;
	size_t capacity;
	int cnt; // used for formatted input

	fast_read(FILE* fin = stdin, size_t buffer_size = 1 << 20) : fin(fin) {
		capacity = buffer_size;
		buff = new char[capacity];
		ptr = end = buff;
	}

	~fast_read() {
		delete[] buff;
	}

	// Refills the buffer.
	// Returns the number of available characters.
	size_t refill() {
		if (ptr >= end) {
			ptr = end = buff;
		}
		end += fread(end, 1, capacity - (end - buff), fin);
		return end - ptr;
	}

	// Tries to reserve at least `size` characters.
	// Returns the number of available characters.
	size_t reserve(size_t size = 1024) {
		size_t available = end - ptr;
		if (available >= size) return available;
		// move the available characters to the
		// beginning of the buffer and refill
		memcpy(buff, ptr, available);
		ptr = buff;
		end = buff + available;
		return refill();
	}
	
	// Reads and returns the next character.
	// Returns -1 if there are no more characters.
	int read_char() {
		if (ptr >= end) {
			refill();
			if (ptr >= end) return -1;
		}
		return *ptr++;
	}

	// Skips whitespaces.
	void skip_whitespaces() {
		int c = read_char();
		while (0x00 <= c && c <= 0x20) { // \t \n ' ' ...
			c = read_char();
		}
		if (c != -1) --ptr; // unread_char();
	}

	// Reads a sign.
	// Character is consumed only if '+' or '-'.
	int read_sign() {
		int c = read_char();
		if (c == '-') return -1;
		if (c == '+') return +1;
		if (c != -1) --ptr; // unread_char();
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
		if (c != -1) --ptr; // unread_char();
		if (cnt) *cnt = l;
		return d;
	}

	// Reads a single signed decimal integer of type I.
	template<typename I>
	I read_integral() {
		skip_whitespaces();
		int s = read_sign();
		I d = read_digits<I>();
		return (s < 0) ? -d : d;
	}

	// Reads a single signed decimal floating-point number of type F.
	template<typename F>
	F read_floating_point() {
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
		if (c != -1) --ptr; // unread_char();
		f *= pow(F(10), e);
		return (s < 0) ? -f : f;
	}

	// Reads a single int.
	int read_int() {
		return read_integral<int>();
	}

	// Reads a single long long.
	long long read_ll() {
		return read_integral<long long>();
	}

	// Reads a single float.
	float read_float() {
		return read_floating_point<float>();
	}

	// Reads a single double.
	double read_double() {
		return read_floating_point<double>();
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
		if (c != -1) --ptr; // unread_char();
		return s;
	}

	// stream overloads
	fast_read& operator >> (char& v) { v = read_char(); return *this; }
	fast_read& operator >> (int& v) { v = read_int(); return *this; }
	fast_read& operator >> (long long& v) { v = read_ll(); return *this; }
	fast_read& operator >> (float& v) { v = read_float(); return *this; }
	fast_read& operator >> (double& v) { v = read_double(); return *this; }
	fast_read& operator >> (std::string& s) { s = read_string(); return *this; }
};

fast_read& fast_in() { static fast_read fin(stdin); return fin; }

}
}
