#pragma once

#include "structure/math/quadratic.h"

namespace altruct {
namespace math {

/**
 * Complex numbers
 */
template<typename T>
using complex = quadratic<T, -1>;

} // math
} // altruct
