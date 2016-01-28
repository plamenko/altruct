#pragma once

#include "structure/math/quadratic.h"

namespace altruct {
namespace math {

/**
 * Complex numbers
 */
template<typename T>
class complex : public quadratic<T, -1> {
public:
	complex(const complex &rhs) : quadratic(rhs) {}
	complex(T a = 0, T b = 0) : quadratic(a, b) {}

};

} // math
} // altruct
