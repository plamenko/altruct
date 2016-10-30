#include <ctime>
#include <iostream>
#include <iomanip>
#include <vector>

#include "algorithm/math/primes.h"
#include "algorithm/math/divisor_sums.h"

using namespace std;
using namespace altruct::math;

void dirichlet_sample() {
	cout << "=== dirichlet_sample ===" << endl;
	
	typedef long long ll;
	int n = 1000000;
	clock_t T;
	vector<int> p(n); T = clock(); int m = primes(p.data(), nullptr, n); printf("primes: %d  %d ms\n", m, clock() - T);
	vector<int> pf(n); T = clock(); factor(pf.data(), n, p.data(), m); printf("factor: %d  %d ms\n", n, clock() - T);
	vector<int> phi(n); T = clock(); euler_phi(phi.data(), n, p.data(), m); printf("phi: %d  %d ms\n", n, clock() - T);
	vector<int> mu(n); T = clock(); moebius_mu(mu.data(), n); printf("moebius: %d  %d ms\n", n, clock() - T);
	vector<int> s0(n); T = clock(); divisor_sigma0(s0.data(), n); printf("divisors: %d  %d ms\n", n, clock() - T);
	vector<ll> s1(n); T = clock(); divisor_sigma1(s1.data(), n); printf("sigma: %d  %d ms\n", n, clock() - T);
	auto _mu = [&](int n){ return mu[n]; };
	auto _amu = [&](int n){ return abs(mu[n]); };
	auto _s1 = [&](int n){ return (int)s1[n]; };
	auto _id = [&](int n){ return n; };
	auto _c1 = [&](int n){ return 1; };

	vector<int> mu_mult(n); // mu = 1^-1
	T = clock(); dirichlet_inverse_multiplicative(mu_mult, _c1, n, p.data(), m);
	printf("moebius mult:  %d  %d ms  %s\n", n, clock() - T, (mu_mult == mu) ? "OK" : "ERR");

	vector<int> ll_mult(n); // liouville_lambda = |mu|^-1
	T = clock(); dirichlet_inverse_multiplicative(ll_mult, _amu, n, p.data(), m);
	printf("liouville mult:  %d  %d ms\n", n, clock() - T);

	vector<int> phi_mult(n); // phi = Id * mu
	T = clock(); dirichlet_convolution_multiplicative(phi_mult, _id, _mu, n, p.data(), m);
	printf("totient mult:  %d  %d ms  %s\n", n, clock() - T, (phi_mult == phi) ? "OK" : "ERR");

	vector<int> s0_mult(n); // d = 1 * 1
	T = clock(); dirichlet_convolution_multiplicative(s0_mult, _c1, _c1, n, p.data(), m);
	printf("divisors mult:  %d  %d ms  %s\n", n, clock() - T, (s0_mult == s0) ? "OK" : "ERR");

	vector<ll> s1_mult(n); // s_k = Id_k * 1
	T = clock(); dirichlet_convolution_multiplicative(s1_mult, _id, _c1, n, p.data(), m);
	printf("sigma mult:  %d  %d ms  %s\n", n, clock() - T, (s1_mult == s1) ? "OK" : "ERR");

	vector<int> id(n);
	T = clock(); dirichlet_convolution(id, _s1, _mu, n);
	printf("convolution gen:  %d  %d ms\n", n, clock() - T);
	T = clock(); dirichlet_convolution_multiplicative(id, _s1, _mu, n, p.data(), m);
	printf("convolution mult: %d  %d ms\n", n, clock() - T);
	T = clock(); dirichlet_convolution_completely_multiplicative(id, _s1, _mu, n, pf.data());
	printf("convolution tot:  %d  %d ms\n", n, clock() - T);

	vector<int> i1(n);
	T = clock(); dirichlet_inverse(i1, _mu, n);
	printf("inverse gen:  %d  %d ms\n", n, clock() - T);
	T = clock(); dirichlet_inverse_multiplicative(i1, _mu, n, p.data(), m);
	printf("inverse mult: %d  %d ms\n", n, clock() - T);
	T = clock(); dirichlet_inverse_completely_multiplicative(i1, _mu, n, pf.data());
	printf("inverse tot:  %d  %d ms\n", n, clock() - T);

	vector<int> tr(n);
	T = clock(); moebius_transform(tr, _id, n);
	printf("moebius trans:  %d  %d ms  %s\n", n, clock() - T, (tr == phi) ? "OK" : "ERR");
	T = clock(); moebius_transform_multiplicative(tr, _id, n, p.data(), m);
	printf("moebius trans mult:  %d  %d ms  %s\n", n, clock() - T, (tr == phi) ? "OK" : "ERR");
	
	cout << endl;
}
