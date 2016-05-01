#pragma once

#include <tuple>
#include "std_hash_combine.h"

namespace std {

namespace {
// Recursive template code derived from Matthieu M.
template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
struct HashValueImpl {
	static void apply(size_t& seed, Tuple const& tuple) {
		HashValueImpl<Tuple, Index - 1>::apply(seed, tuple);
		hash_combine(seed, std::get<Index>(tuple));
	}
};

template <class Tuple>
struct HashValueImpl<Tuple, 0> {
	static void apply(size_t& seed, Tuple const& tuple){
		hash_combine(seed, std::get<0>(tuple));
	}
};
} // anonymous

template <typename ... TT>
struct hash<std::tuple<TT...>> {
	size_t operator()(std::tuple<TT...> const& tt) const {
		size_t seed = 0;
		HashValueImpl<std::tuple<TT...> >::apply(seed, tt);
		return seed;
	}
};

} // std
