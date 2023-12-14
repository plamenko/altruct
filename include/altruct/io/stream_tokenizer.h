#pragma once

#include <string>
#include <iostream>
#include <sstream>
#include <vector>

namespace altruct {
namespace io {
namespace stream_tokenizer {

template<typename P>
struct token {
    std::string s;
};

template<typename T, typename D>
struct tokens {
    std::vector<T> v;
};

template<typename T, typename P, typename D, typename F>
std::vector<T> unbox_tokens(const tokens<token<P>, D>& ts, F unbox_f) {
    std::vector<T> vt;
    for (const auto& t : ts.v) {
        vt.push_back(unbox_f(t.s));
    }
    return vt;
}
    
template<typename P, typename D>
std::vector<std::string> unbox_tokens(const tokens<token<P>, D>& ts) {
    auto identity_f = [](const std::string& s) { return s; };
    return unbox_tokens<std::string>(ts, identity_f);
}

template<char... Ds>
struct delimited_p {
    static bool eval(char c) { return (... && (c != Ds)); }
};

template<char... As>
struct allowed_p {
    static bool eval(char c) { return (... || (c == As)); }
};

template<typename P>
struct not_p {
    static bool eval(char c) { return !P::eval(c); }
};

struct alphanum_p {
    static bool eval(char c) { return std::isalnum(c); }
};

template<typename P>
using tokens_allowed = tokens<token<P>, token<not_p<P>>>;

template<char... Ds>
using tokens_delimited = tokens<token<delimited_p<Ds...>>, token<allowed_p<Ds...>>>;
//using tokens_delimited = tokens_allowed<delimited_p<Ds...>>;

using token_alphanum = token<alphanum_p>;

using token_binary = token<allowed_p<'0', '1'>>;

using token_delimited_comma = token<delimited_p<','>>;
using token_delimited_semicolon = token<delimited_p<';'>>;

using tokens_delimited_space = tokens_delimited<' '>;
using tokens_delimited_comma = tokens_delimited<','>;
using tokens_delimited_semicolon = tokens_delimited<';'>;

using ints_delimited_space = tokens<int, char>;
using int64s_delimited_space = tokens<int64_t, char>;

} // stream_tokenizer
} // io
} // altruct

template<typename P>
std::istream& operator>>(std::istream& is, altruct::io::stream_tokenizer::token<P>& t) {
    char ch;
    t.s.clear();
    is >> std::ws;
    while (is.get(ch)) {
        if (P::eval(ch)) {
            t.s.push_back(ch);
        } else {
            is.unget();
            break;
        }
    }
    return is;
}

template<typename T, typename D>
std::istream& operator>>(std::istream& is, altruct::io::stream_tokenizer::tokens<T, D>& ts) {
    T t; D d;
    ts.v.clear();
    while (is) {
        is >> t;
        ts.v.push_back(t);
        bool wasSkipWs = is.flags() & std::ios_base::skipws;
        is >> std::noskipws >> d;
        if (wasSkipWs) is >> std::skipws;
    }
    return is;
}
