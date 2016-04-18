#include "algorithm/math/modulos.h"

namespace altruct {
namespace math {

int primitive_root(int m, prime_holder& prim) {
	int phi = euler_phi(prim.factor_integer(m));
	auto phi_factors = prime_factors(prim.factor_integer(phi));
	return primitive_root(m, phi, phi_factors);
}

std::set<int> kth_roots(int m, int k, prime_holder& prim) {
	int lam = carmichael_lambda(prim.factor_integer(m));
	int g = primitive_root(m, prim);
	return kth_roots(m, k, lam, g);
}

} // math
} // altruct
