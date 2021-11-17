#pragma once

#include "altruct/algorithm/math/base.h"

#include <string>
#include <sstream>

namespace altruct {
namespace math {

class symbolic {
public:
    std::string v;

    symbolic(std::string name = "?") : v(std::move(name)) {}

    symbolic& operator += (const symbolic& s) { return *this = *this + s; }
    symbolic& operator -= (const symbolic& s) { return *this = *this - s; }
    symbolic& operator *= (const symbolic& s) { return *this = *this * s; }
    symbolic& operator /= (const symbolic& s) { return *this = *this / s; }
    symbolic& operator %= (const symbolic& s) { return *this = *this % s; }
    symbolic& operator &= (const symbolic& s) { return *this = *this & s; }
    symbolic& operator |= (const symbolic& s) { return *this = *this | s; }
    symbolic& operator ^= (const symbolic& s) { return *this = *this ^ s; }
    symbolic& operator <<= (const symbolic& s) { return *this = *this << s; }
    symbolic& operator >>= (const symbolic& s) { return *this = *this >> s; }

    symbolic  operator +  ()                  const { return symbolic("(+" + v + ")"); }
    symbolic  operator -  ()                  const { return symbolic("(-" + v + ")"); }
    symbolic  operator +  (const symbolic& s) const { return symbolic("(" + v + "+" + s.v + ")"); }
    symbolic  operator -  (const symbolic& s) const { return symbolic("(" + v + "-" + s.v + ")"); }
    symbolic  operator *  (const symbolic& s) const { return symbolic("(" + v + "*" + s.v + ")"); }
    symbolic  operator /  (const symbolic& s) const { return symbolic("(" + v + "/" + s.v + ")"); }
    symbolic  operator %  (const symbolic& s) const { return symbolic("(" + v + "%" + s.v + ")"); }

    symbolic  operator ~  ()                  const { return symbolic("(~" + v + ")"); }
    symbolic  operator &  (const symbolic& s) const { return symbolic("(" + v + "&" + s.v + ")"); }
    symbolic  operator |  (const symbolic& s) const { return symbolic("(" + v + "|" + s.v + ")"); }
    symbolic  operator ^  (const symbolic& s) const { return symbolic("(" + v + "^" + s.v + ")"); }
    symbolic  operator << (const symbolic& s) const { return symbolic("(" + v + "<<" + s.v + ")"); }
    symbolic  operator >> (const symbolic& s) const { return symbolic("(" + v + ">>" + s.v + ")"); }

    symbolic  operator !  ()                  const { return symbolic("(!" + v + ")"); }
    symbolic  operator && (const symbolic& s) const { return symbolic("(" + v + "&&" + s.v + ")"); }
    symbolic  operator || (const symbolic& s) const { return symbolic("(" + v + "||" + s.v + ")"); }
};

template<typename T>
struct castT<symbolic, T> {
    static symbolic of(const T& v) {
        std::stringstream ss; ss << v;
        return symbolic(ss.str());
    }
    static symbolic of(const symbolic& ref, const T& v) {
        return of(v);
    }
};
template<>
struct castT<symbolic, symbolic> : nopCastT<symbolic> {};

template<>
struct identityT<symbolic> {
    static symbolic of(const symbolic& s) {
        return symbolic("1");
    }
};

template<>
struct zeroT<symbolic> {
    static symbolic of(const symbolic& s) {
        return symbolic("0");
    }
};

} // math
} // altruct
