#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>

#include "structure/math/fraction.h"
#include "structure/math/modulo.h"
#include "structure/math/polynom.h"
#include "structure/math/matrix.h"

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


/** std::ostream specialization for altruct::math::fraction */
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::fraction<T>& rhs) {
    return os << rhs.p << "/" << rhs.q;
}

/** std::ostream specialization for altruct::math::modulo */
template<typename T, int ID, int STORAGE_TYPE>
std::ostream& operator << (std::ostream& os, const altruct::math::modulo<T, ID, STORAGE_TYPE>& rhs) {
    return os << rhs.v;
}

/** std::ostream specialization for altruct::math::polynom */
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::polynom<T>& rhs) {
    return os << rhs.c;
}

/** std::ostream specialization for altruct::math::matrix */
template<typename T>
std::ostream& operator << (std::ostream& os, const altruct::math::matrix<T>& rhs) {
    return os << rhs.a;
}
