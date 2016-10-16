#include "algorithm/math/modulos.h"

namespace altruct {
namespace math {

int primitive_root(int m, prime_holder& prim) {
	int phi = euler_phi(prim.factor_integer(m));
	auto phi_factors = prime_factors(prim.factor_integer(phi));
	return primitive_root(m, phi, phi_factors);
}

int primitive_root_of_unity(int m, prime_holder& prim) {
	int lam = carmichael_lambda(prim.factor_integer(m));
	auto lam_factors = prime_factors(prim.factor_integer(lam));
	return primitive_root_of_unity(m, lam, lam_factors);
}

std::set<int> kth_roots_of_unity(int m, int k, prime_holder& prim) {
	int lam = carmichael_lambda(prim.factor_integer(m));
	auto lam_factors = prime_factors(prim.factor_integer(lam));
	int g = primitive_root_of_unity(m, lam, lam_factors);
	return kth_roots_of_unity(m, k, lam, g);
}

} // math
} // altruct
