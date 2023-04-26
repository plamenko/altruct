#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <iostream>

namespace altruct {
namespace io {

/**
 * A class that buffers the input file and provides basic input facilities.
 *
 * Example of simple input:
 *   int x = fin.read_int();
 *   string s = fin.read_string();
 *
 * Example of streamed input:
 *   fin >> x >> s;
 *
 * Example of formatted input:
 *   sscanf_s(fin.data(100), "%x %d %f %n", &v1, &v2, &v3, &fin.counter());
 *   fin.advance();
 */
class fast_read {
    FILE* fin;
    char* buff;
    char* ptr;
    char* end;
    size_t capacity;
    int cnt; // used for formatted input

public:
    fast_read(FILE* fin = stdin, size_t buffer_size = 1 << 20) : fin(fin) {
        capacity = buffer_size;
        buff = new char[capacity+1];
        buff[capacity] = 0;
        ptr = end = buff;
        cnt = 0;
    }

    ~fast_read() {
        delete[] buff;
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
            size_t len = fread(end, 1, capacity - available, fin);
            if (len == 0) break;
            available += len;
            end += len;
        }
        return available;
    }

    // Tries to reserve at least `count` characters.
    // Returns the number of available characters.
    size_t reserve(size_t count = 1024) {
        size_t available = end - ptr;
        return (available >= count) ? available : refill();
    }

    // Returns the number of available characters.
    size_t available() {
        return end - ptr;
    }

    // Returns a pointer to the raw character buffer,
    // provided that there are at least `count` characters
    // available, or `nullptr` otherwise.
    const char* data(size_t count = 0) {
        size_t available = reserve(count);
        return (available >= count) ? ptr : nullptr;
    }

    // Returns the reference to `cnt` variable.
    // This is convenient for `sscanf` usages.
    int& counter() {
        return cnt;
    }

    // Skips `cnt` characters.
    // This is convenient for `sscanf` usages.
    void advance() {
        skip(cnt);
        cnt = 0;
    }

    // Skips `count` characters.
    void skip(size_t count) {
        ptr += count;
    }

    // Unreads the last read character.
    // Allowed to be called only immediately
    // after a sucessfull `read_char()`.
    void unread_char() {
        --ptr;
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

    // stream overloads
    fast_read& operator >> (char& v) { v = read_char(); return *this; }
    fast_read& operator >> (std::string& s) { s = read_string(); return *this; }
    fast_read& operator >> (int& v) { v = read_int(); return *this; }
    fast_read& operator >> (long& v) { v = read_l(); return *this; }
    fast_read& operator >> (long long& v) { v = read_ll(); return *this; }
    fast_read& operator >> (float& v) { v = read_float(); return *this; }
    fast_read& operator >> (double& v) { v = read_double(); return *this; }
};

/**
 * A class that buffers the output file and provides basic output facilities.
 *
 * Example of simple output:
 *   fout.write_int(42);
 *   fout.write_string("abc");
 *
 * Example of streamed output:
 *   fout << 42 << "abc";
 *
 * Example of formatted output:
 *   fout.reserve(100);
 *   sprintf_s(fout.data(), fout.available(), "some numbers: %x %d %f", 0x123, 42, 123.45f);
 *   fout.advance();
 */
class fast_write {
    FILE* fout;
    char* buff;
    char* ptr;
    char* end;
    size_t capacity;
    int cnt; // used for formatted input

public:
    fast_write(FILE* fout = stdout, size_t buffer_size = 1 << 20) : fout(fout) {
        capacity = buffer_size;
        buff = new char[capacity + 1];
        buff[capacity] = 0;
        ptr = buff;
        end = buff + capacity;
    }

    ~fast_write() {
        flush();
        delete[] buff;
    }

    // Flushes the buffer and the underlying file.
    size_t flush() {
        fwrite(buff, 1, ptr - buff, fout);
        fflush(fout);
        ptr = buff;
        return end - ptr;
    }

    // Tries to reserve space for at least `count` characters.
    // Returns the number of characters that can be written.
    size_t reserve(size_t count = 1024) {
        size_t available = end - ptr;
        return (available >= count) ? available : flush();
    }

    // Returns the number of characters that can be written.
    size_t available() {
        return end - ptr;
    }

    // Returns a pointer to the raw character buffer,
    // provided that there is space for at least `count` characters
    // available to be written, or `nullptr` otherwise.
    char* data(size_t count = 0) {
        size_t available = reserve(count);
        return (available >= count) ? ptr : nullptr;
    }

    // Advances the write pointer till the '\0'.
    // This is convenient for `sprintf` usages.
    fast_write& advance() {
        ptr += strlen(ptr);
        return *this;
    }

    fast_write& write_char(int c) {
        if (ptr >= end) {
            flush();
        }
        *ptr++ = (char)c;
        return *this;
    }

    fast_write& write_string(const char *s) {
        for (int i = 0; s[i]; i++) {
            write_char(s[i]);
        }
        return *this;
    }

    fast_write& write_string(const std::string& s) {
        return write_string(s.c_str());
    }

    template<typename I>
    fast_write& write_integral_number(I d, int len = 1) {
        int i;
        static char tc[50];
        if (d < 0) {
            d = -d;
            write_char('-');
        }
        for (i = 0; (d > 0 || i < len) && i < sizeof(tc); i++) {
            tc[i] = (d % 10) + '0';
            d /= 10;
        }
        while (i > 0) {
            write_char(tc[--i]);
        }
        return *this;
    }

    template<typename F>
    fast_write& write_floating_point_number(F f, int precision = 6, bool scientific = false) {
        if (f < 0) {
            f = -f;
            write_char('-');
        }
        int e = 1, g = exponent(f);
        if (scientific) {
            f /= pow(F(10), g);
        } else {
            e = std::max(g, 0) + 1;
        }
        f += pow(F(10), -precision) / 2; // round
        while (e > -precision) {
            if (e == 0) {
                write_char('.');
            }
            e--;
            F w = pow(F(10), e);
            int d = (int)(f / w);
            f -= w * d;
            write_char(d + '0');
        };
        if (scientific) {
            write_char('e');
            write_char((g < 0) ? '-' : '+');
            write_integral_number<int>((g < 0) ? -g : +g, 3);
        }
        return *this;
    }

    template<typename F>
    int exponent(F f) {
        return (f == 0) ? 0 : (int)floor(log10(f));
    }

    fast_write& write_int(int d) {
        return write_integral_number<int>(d);
    }

    fast_write& write_ll(long long d) {
        return write_integral_number<long long>(d);
    }

    fast_write& write_float(float f, int precision = 6, bool scientific = false) {
        // use double internally
        return write_floating_point_number<double>(f, precision, scientific);
    }

    fast_write& write_double(double f, int precision = 6, bool scientific = false) {
        return write_floating_point_number<double>(f, precision, scientific);
    }

    // stream overloads
    fast_write& operator << (const char& v) { write_char(v); return *this; }
    fast_write& operator << (const std::string& s) { write_string(s); return *this; }
    fast_write& operator << (const int& v) { write_int(v); return *this; }
    fast_write& operator << (const long& v) { write_ll(v); return *this; }
    fast_write& operator << (const long long& v) { write_ll(v); return *this; }
    fast_write& operator << (const float& v) { write_float(v); return *this; }
    fast_write& operator << (const double& v) { write_double(v); return *this; }
};

template<typename V = void> fast_read& fast_in() { static fast_read fin(stdin); return fin; }
template<typename V = void> fast_write& fast_out() { static fast_write fout(stdout); return fout; }
template<typename V = void> fast_write& fast_err() { static fast_write ferr(stderr); return ferr; }

}
}
