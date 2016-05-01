#pragma once

#include <cstddef>

namespace std {

#ifndef hash_combine
// Code from Boost
// Reciprocal of the golden ratio helps spread entropy and handles duplicates.
// See Mike Seymour in magic-numbers-in-boosthash-combine:
// http://stackoverflow.com/questions/4948780
template <class T>
inline void hash_combine(std::size_t& seed, T const& v) {
	seed ^= std::hash<T>()(v)+0x9e3779b9 + (seed << 6) + (seed >> 2);
}
#endif // hash_combine

} // std
