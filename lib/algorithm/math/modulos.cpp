#include "algorithm/math/modulos.h"

namespace altruct {
namespace math {

int primitive_root(int m, prime_holder& prim) {
	int phi = euler_phi(prim.factor_integer(m));
	auto phi_factors = prime_factors(prim.factor_integer(phi));
	return primitive_root(m, phi, phi_factors);
}

std::vector<int> kth_roots(int m, int k, prime_holder& prim) {
	auto vf = prim.factor_integer(m);
	int lam = carmichael_lambda(vf);
	int phi = euler_phi(vf);
	auto phi_factors = prime_factors(prim.factor_integer(phi));
	return kth_roots(m, k, lam, phi, phi_factors);
}

} // math
} // altruct
