#pragma once

#include <cmath>
#include <cstdio>
#include <string>
#include <iostream>
#include <algorithm>

namespace altruct {
namespace io {

/**
 * A character writer interface.
 */
class writer {
public:
    // Writes a character.
    virtual writer& write_char(int c) = 0;

    // Writes `count` characters from `buffer`.
    virtual writer& write(const char *buffer, size_t count) = 0;

    // Flushes the written characters.
    virtual writer& flush() = 0;

    virtual ~writer() {}
    writer() {}
    writer(const writer&) = delete;
    writer& operator = (const writer&) = delete;
};

/**
 * A character writer implementation backed by a file.
 */
class file_writer : public writer {
    FILE* out;

public:
    // Constructs a new `file_writer` on top of `out` file.
    // `out` pointer is owned by the caller!
    file_writer(FILE* out) : out(out) {}

    ~file_writer() {
        flush();
    }

    writer& write_char(int c) override {
        fputc(c, out);
        return *this;
    }

    writer& write(const char *buffer, size_t count) override {
        fwrite(buffer, 1, count, out);
        return *this;
    }

    writer& flush() override {
        fflush(out);
        return *this;
    }
};

/**
 * A character writer implementation backed by an output stream.
 */
class stream_writer : public writer {
    std::ostream& out;

public:
    // Constructs a new `stream_writer` on top of `out` stream.
    // The caller needs to keep the underlying stream alive!
    stream_writer(std::ostream& out) : out(out) {}

    ~stream_writer() {
        flush();
    }

    stream_writer& write_char(int c) override {
        out.put(c);
        return *this;
    }

    stream_writer& write(const char *buffer, size_t count) override {
        out.write(buffer, count);
        return *this;
    }

    stream_writer& flush() override {
        out.flush();
        return *this;
    }
};

/**
 * A character writer implementation backed by a char array.
 */
class string_writer : public writer {
    char* out;
    size_t capacity;
    size_t pos;

public:
    // Constructs a new `string_writer` on top of `out` string.
    // `out` pointer is owned by the caller!
    string_writer(char* out, size_t capacity) : out(out), capacity(capacity), pos(0) {}

    string_writer& write_char(int c) override {
        if (pos >= capacity) return *this;
        out[pos++] = (char)c;
        return *this;
    }

    string_writer& write(const char *buffer, size_t count) override {
        if (pos >= capacity) return *this;
        count = std::min(count, capacity - pos);
        std::copy(buffer, buffer + count, out + pos);
        pos += count;
        return *this;
    }

    string_writer& flush() override {
        return *this;
    }
};

/**
 * A class that buffers the underlying writer.
 *
 * Example of formatted output:
 *   fout.reserve(100);
 *   sprintf_s(fout.data(), fout.available(), "some numbers: %x %d %f", 0x123, 42, 123.45f);
 *   fout.advance();
 */
class buffered_writer : public writer {
    writer& out;
    char* buff;
    char* ptr;
    char* end;
    size_t capacity;
    int cnt; // used for formatted input

public:
    // Constructs a new `buffered_writer` on top of `out` writer.
    // The caller needs to keep the underlying writer alive!
    buffered_writer(writer& out, size_t buffer_size = 1 << 20) : out(out), cnt(0) {
        capacity = buffer_size;
        buff = new char[capacity + 1];
        buff[capacity] = 0;
        ptr = buff;
        end = buff + capacity;
    }

    ~buffered_writer() {
        flush();
        delete[] buff;
    }

    // Flushes the buffer and the underlying file.
    buffered_writer& flush() override {
        out.write(buff, ptr - buff);
        out.flush();
        ptr = buff;
        *ptr = 0;
        return *this;
    }

    // Tries to reserve space for at least `count` characters.
    // Returns the number of characters that can be written.
    size_t reserve(size_t count = 1024) {
        size_t available = end - ptr;
        return (available >= count) ? available : flush().available();
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
    buffered_writer& advance() {
        ptr += strlen(ptr);
        return *this;
    }

    buffered_writer& write_char(int c) override {
        if (ptr >= end) {
            flush();
        }
        *ptr++ = (char)c;
        return *this;
    }

    buffered_writer& write(const char *buffer, size_t count) override {
        while (count > 0) {
            size_t available = end - ptr;
            size_t len = std::min(available, count);
            std::copy(buffer, buffer + len, ptr);
            ptr += len;
            buffer += len;
            count -= len;
            if (count > 0) flush();
        }
        return *this;
    }
};

/**
 * A class that buffers the underlying writer and provides basic output facilities.
 *
 * Example of simple output:
 *   fout.write_int(42);
 *   fout.write_string("abc");
 *
 * Example of formatted output:
 *   fout.reserve(100);
 *   sprintf_s(fout.data(), fout.available(), "some numbers: %x %d %f", 0x123, 42, 123.45f);
 *   fout.advance();
 */
class simple_writer : public buffered_writer {
public:
    // Constructs a new `simple_writer` on top of `out` writer.
    // The caller needs to keep the underlying writer alive!
    simple_writer(writer& out, size_t buffer_size = 1 << 20) : buffered_writer(out, buffer_size) {
    }

    simple_writer& write_string(const char *s) {
        write(s, strlen(s));
        return *this;
    }

    simple_writer& write_string(const std::string& s) {
        write(s.c_str(), s.size());
        return *this;
    }

    template<typename I>
    simple_writer& write_integral_number(I d, int len = 1) {
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
    simple_writer& write_floating_point_number(F f, int precision = 6, bool scientific = false) {
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

    simple_writer& write_int(int d) {
        return write_integral_number<int>(d);
    }

    simple_writer& write_ll(long long d) {
        return write_integral_number<long long>(d);
    }

    simple_writer& write_float(float f, int precision = 6, bool scientific = false) {
        // use double internally
        return write_floating_point_number<double>(f, precision, scientific);
    }

    simple_writer& write_double(double f, int precision = 6, bool scientific = false) {
        return write_floating_point_number<double>(f, precision, scientific);
    }
};

/**
 * A class that provides streaming interface on top of `simple_writer`.
 *
 * Example of streamed output:
 *   out << x << s;
 */
class simple_writer_stream {
    simple_writer& out;

public:
    // Constructs a new `simple_writer_stream` on top of `out` writer.
    // The caller needs to keep the underlying writer alive!
    simple_writer_stream(simple_writer& out) : out(out) {}

    simple_writer_stream& operator << (const char& v) { out.write_char(v); return *this; }
    simple_writer_stream& operator << (const std::string& s) { out.write_string(s); return *this; }
    simple_writer_stream& operator << (const int& v) { out.write_int(v); return *this; }
    simple_writer_stream& operator << (const long& v) { out.write_ll(v); return *this; }
    simple_writer_stream& operator << (const long long& v) { out.write_ll(v); return *this; }
    simple_writer_stream& operator << (const float& v) { out.write_float(v); return *this; }
    simple_writer_stream& operator << (const double& v) { out.write_double(v); return *this; }
};

template<typename V = void> simple_writer& simple_out_writer() {
    static file_writer fout(stdout);
    static simple_writer sout(fout);
    return sout;
}
template<typename V = void> simple_writer_stream& simple_out_stream() {
    static simple_writer_stream sout(simple_out_writer());
    return sout;
}

template<typename V = void> simple_writer& simple_err_writer() {
    static file_writer ferr(stderr);
    static simple_writer serr(ferr);
    return serr;
}
template<typename V = void> simple_writer_stream& simple_err_stream() {
    static simple_writer_stream serr(simple_err_writer());
    return serr;
}

}
}
