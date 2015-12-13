#include "structure/math/modulo.h"

namespace altruct {
namespace math {

// integral type specializations
// - to avoid negative numbers
// - to avoid overflow for big moduli
// - faster division

// modulo normalization
void modulo_normalize(long long& v, long long M) { modulo_normalize_int(v, M); }
void modulo_normalize(int& v, int M) { modulo_normalize_int(v, M); }

// modulo multiplication
long long modulo_multiply(long long x, long long y, long long M) {
	modulo_normalize_int(x, M);
	modulo_normalize_int(y, M);
	bool fit = (x < (1LL << 31) && y < (1LL << 31));
	return fit ? (x * y) % M : modulo_long_multiply_int(x, y, M);
}
int modulo_multiply(int x, int y, int M) {
	modulo_normalize_int(x, M);
	modulo_normalize_int(y, M);
	return ((long long)x * y) % M;
}

// modulo inversion
long long modulo_inverse(long long v, long long M) { return modulo_inverse_int(v, M); }
int modulo_inverse(int v, int M) { return modulo_inverse_int(v, M); }

// modulo division
long long modulo_divide(long long x, long long y, long long M) { return modulo_divide_int(x, y, M); }
int modulo_divide(int x, int y, int M) { return modulo_divide_int(x, y, M); }

} // math
} // altruct
