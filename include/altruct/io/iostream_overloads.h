#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>

#include "altruct/algorithm/math/triples.h"
#include "altruct/structure/math/fraction.h"
#include "altruct/structure/math/modulo.h"
#include "altruct/structure/math/polynom.h"
#include "altruct/structure/math/matrix.h"
#include "altruct/structure/math/clifford3.h"
#include "altruct/structure/math/pga.h"
#include "altruct/structure/math/symbolic.h"

/** Forward declarations */
template<typename T1, typename T2>
std::ostream& operator << (std::ostream& os, const std::pair<T1, T2>& rhs);
template<typename T, typename A>
std::ostream& operator << (std::ostream& os, const std::vector<T, A>& container);
template<typename T, typename P, typename A>
std::ostream& operator << (std::ostream& os, const std::set<T, P, A>& container);
template<typename K, typename V, typename P, typename A>
std::ostream& operator << (std::ostream& os, const std::map<K, V, P, A>& container);

template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::triple<T>& rhs);
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::fraction<T>& rhs);
template<typename T, int ID, int STORAGE_TYPE>
std::ostream& operator << (std::ostream& os, const altruct::math::modulo<T, ID, STORAGE_TYPE>& rhs);
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::polynom<T>& rhs);
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::matrix<T>& rhs);
template<typename T=std::string>
std::ostream& operator << (std::ostream& os, const altruct::math::symbolic& s);

template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::cl3::rotor<T>& r);
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::cl3::vector<T>& v);
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::cl3::multivector<T>& m);

template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::zero<T>& z);
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::blade0<T>& b0);
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::blade1<T>& b1);
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::blade2E<T>& b2E);
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::blade2e<T>& b2e);
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::blade3<T>& b3);
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::blade4<T>& b4);
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::blade02<T>& b02);
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::blade24<T>& b24);
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::blade024<T>& b);
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::blade13<T>& b);
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::multivector<T>& m);

/** std::ostream manipulator base */
struct altruct_io_manipulator_base {};
template<typename T = void>
std::ostream& operator << (std::ostream& os, const altruct_io_manipulator_base& f) { return os; }
/** std::ostream manipulator macro */
#define ALTRUCT_IO_MANIPULATOR(manipulator_name, field, default_val) \
struct manipulator_name : public altruct_io_manipulator_base { \
    decltype(field) old_val; \
    manipulator_name(decltype(field) val) : old_val(field) { field = val; } \
    ~manipulator_name() { field = old_val; } \
    static bool set_default() { field = default_val; return true; } \
};
//bool manipulator_name##_initialized = manipulator_name::set_default();


/** std::ostream specialization for std::pair */
template<typename T1, typename T2>
std::ostream& operator << (std::ostream& os, const std::pair<T1, T2>& rhs) {
    return os << "{" << rhs.first << ", " << rhs.second << "}";
}

/** std::ostream specialization for containers */
template<typename C>
std::ostream& output_container(std::ostream& os, const C& container, const std::string& delimiter = ", ") {
    os << "{";
    bool first = true;
    for (const auto& elem : container) {
        if (!first) os << delimiter;
        os << elem;
        first = false;
    }
    os << "}";
    return os;
}

/** std::ostream specialization for std::vector */
template<typename T, typename A>
std::ostream& operator << (std::ostream& os, const std::vector<T, A>& container) { return output_container(os, container); }
/** std::ostream specialization for std::set */
template<typename T, typename P, typename A>
std::ostream& operator << (std::ostream& os, const std::set<T, P, A>& container) { return output_container(os, container); }
/** std::ostream specialization for std::map */
template<typename K, typename V, typename P, typename A>
std::ostream& operator << (std::ostream& os, const std::map<K, V, P, A>& container) { return output_container(os, container); }

/** std::ostream specialization for altruct::math::triple */
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::triple<T>& t) {
    return os << "{" << t.a << ", " << t.b << ", " << t.c << "}";
}

/** std::ostream specialization for altruct::math::fraction */
struct iostream_fraction_state_t {
    bool always_output_denominator;
    bool output_as_pair;
};
extern iostream_fraction_state_t iostream_fraction_state;
ALTRUCT_IO_MANIPULATOR(io_fraction_denominator, iostream_fraction_state.always_output_denominator, false);
ALTRUCT_IO_MANIPULATOR(io_fraction_as_pair, iostream_fraction_state.output_as_pair, false);
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::fraction<T>& rhs) {
    if (iostream_fraction_state.output_as_pair) {
        return os << std::make_pair(rhs.p, rhs.q);
    } else if (iostream_fraction_state.always_output_denominator || !(rhs.q == altruct::math::identityOf(rhs.q))) {
        return os << rhs.p << "/" << rhs.q;
    } else {
        return os << rhs.p;
    }
}

/** std::ostream specialization for altruct::math::modulo */
struct iostream_modulo_state_t {
    bool output_modulus;
    bool output_as_pair;
};
extern iostream_modulo_state_t iostream_modulo_state;
ALTRUCT_IO_MANIPULATOR(io_modulo_modulus, iostream_modulo_state.output_modulus, false);
ALTRUCT_IO_MANIPULATOR(io_modulo_as_pair, iostream_modulo_state.output_as_pair, false);
template<typename T, int ID, int STORAGE_TYPE>
std::ostream& operator << (std::ostream& os, const altruct::math::modulo<T, ID, STORAGE_TYPE>& rhs) {
    if (iostream_modulo_state.output_as_pair) {
        return os << std::make_pair(rhs.v, rhs.M());
    } else if (iostream_modulo_state.output_modulus) {
       return os << rhs.v << " (mod " << rhs.M() << ")";
    } else {
        return os << rhs.v;
    }
}

/** std::ostream specialization for altruct::math::polynom */
struct iostream_polynom_state_t {
    bool output_as_vector;
};
extern iostream_polynom_state_t iostream_polynom_state;
ALTRUCT_IO_MANIPULATOR(io_polynom_as_vector, iostream_polynom_state.output_as_vector, true);
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::polynom<T>& rhs) {
    if (iostream_polynom_state.output_as_vector) {
        return os << rhs.c;
    } else {
        int h = rhs.deg(), l = rhs.lowest();
        for (int i = h; i >= l; i--) {
            if (i < h && rhs[i] == rhs.ZERO_COEFF) continue;
            if (i < h) os << " + ";
            os << rhs[i];
            if (i > 0) os << " x";
            if (i > 1) os << "^" << i;
        }
        return os;
    }
}

/** std::ostream specialization for altruct::math::matrix */
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::matrix<T>& rhs) {
    return os << rhs.a;
}

/** std::ostream specialization for altruct::math::symbolic */
template<typename T = std::string>
std::ostream& operator << (std::ostream& os, const altruct::math::symbolic& s) {
    return os << s.v;
}

/** std::ostream specialization for altruct::math::cl3 */
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::cl3::rotor<T>& r) {
    return os << "{" << r.s << ", " << r.yz << ", " << r.zx << ", " << r.xy << "}";
}
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::cl3::vector<T>& v) {
    return os << "{" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << "}";
}
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::cl3::multivector<T>& m) {
    return os << "{" << m.r << ", " << m.v << "}";
}

/** std::ostream specializations for altruct::math::pga */
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::blade0<T>& b0) {
    return os << b0.s << " id";
}
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::zero<T>& z) {
    return os << "0";
}
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::blade1<T>& b1) {
    return os << b1.e0 << " e0" << " + " << b1.v.x << " e1" << " + " << b1.v.y << " e2" << " + " << b1.v.z << " e3";
}
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::blade2E<T>& b2E) {
    return os << b2E.biE.x << " e23" << " + " << b2E.biE.y << " e31" << " + " << b2E.biE.z << " e12";
}
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::blade2e<T>& b2e) {
    return os << b2e.bie.x << " e01" << " + " << b2e.bie.y << " e02" << " + " << b2e.bie.z << " e03";
}
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::blade3<T>& b3) {
    return os << b3.e123 << " e123" << " + " << b3.triP.x << " e032" << " + " << b3.triP.y << " e013" << " + " << b3.triP.z << " e021";
}
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::blade4<T>& b4) {
    return os << b4.e0123 << " e0123";
}
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::blade02<T>& b02) {
    return os << b02.s << " id" << " + " << b02.biE.x << " e23" << " + " << b02.biE.y << " e31" << " + " << b02.biE.z << " e12";
}
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::blade24<T>& b24) {
    return os << b24.bie.x << " e01" << " + " << b24.bie.y << " e02" << " + " << b24.bie.z << " e03" << " + " << b24.e0123 << " e0123";
}
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::blade024<T>& b) {
    return os << b.b02 << " + " << b.b24;
}
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::blade13<T>& b) {
    return os << b.b1 << " + " << b.b3;
}
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::pga::multivector<T>& m) {
    return os << m.b13 << " + " << m.b024;
}
