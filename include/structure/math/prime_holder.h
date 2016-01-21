#pragma once

#include "algorithm/math/primes.h"

#include <vector>

namespace altruct {
namespace math {

class prime_holder {
private:
	typedef std::pair<int, int> fact_pair;

	int sz;                // upper bound (exclusive)
	int m;                 // number of primes up to sz
	std::vector<int> vp;   // primes
	std::vector<char> vq;  // prime flags
	std::vector<int> vpf;  // biggest prime factor
	std::vector<int> vpi;  // prime pi
	std::vector<int> vphi; // euler phi (totient)
	std::vector<int> vmu;  // moebius mu

	void ensure_pq();
	std::vector<int>& ensure(std::vector<int> &v, void(*f)(int*, int, int*, int));

public:
	prime_holder(int sz) : sz(sz) {}

	int size() { return sz; }
	int primes() { ensure_pq(); return m; }

	std::vector<int>& p() { ensure_pq(); return vp; }
	std::vector<char>& q() { ensure_pq(); return vq; }
	std::vector<int>& pf() { return ensure(vpf, altruct::math::factor); }
	std::vector<int>& pi() { return ensure(vpi, altruct::math::prime_pi); }
	std::vector<int>& phi() { return ensure(vphi, altruct::math::euler_phi); }
	std::vector<int>& mu() { return ensure(vmu, altruct::math::moebius_mu); }

	int p(int i) { return p().at(i); }
	int q(int i) { return q().at(i); }
	int pf(int i) { return pf().at(i); }
	int pi(int i) { return pi().at(i); }
	int phi(int i) { return phi().at(i); }
	int mu(int i) { return mu().at(i); }

	std::vector<fact_pair> factor_integer(int n);
	std::vector<fact_pair> factor_integer(std::vector<int> vn);
	std::vector<long long> divisors(const std::vector<fact_pair> &vf, long long maxd = LLONG_MAX);
	std::vector<long long> divisors(const std::vector<int> &vn, long long maxd = LLONG_MAX);
	std::vector<long long> divisors(int n, long long maxd = LLONG_MAX);
};

} // math
} // altruct
