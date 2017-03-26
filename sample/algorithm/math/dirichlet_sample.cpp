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
	
    char* message_n = "O(n)";
    char* message_log = "O(n log n)";
    char* message_log_log = "O(n log log n)";
    char* message_log_log_vs = "O(n log log n) instead of O(n log n)";

    auto sec_since = [](clock_t T){ return double(clock() - T) / CLOCKS_PER_SEC; };

    typedef long long ll;
	int n = 10000000;
	clock_t T;
    vector<int> p(n); T = clock(); int m = primes(p.data(), nullptr, n);
    printf("primes:              %d  %0.3lf s  OK  %s\n", n, sec_since(T), message_log_log);
    vector<int> pf(n); T = clock(); factor(pf.data(), n, p.data(), m);
    printf("factor:              %d  %0.3lf s  OK  %s\n", n, sec_since(T), message_log_log);
    vector<int> phi(n); T = clock(); euler_phi(phi.data(), n, p.data(), m);
    printf("phi:                 %d  %0.3lf s  OK  %s\n", n, sec_since(T), message_log_log);
    vector<int> mu(n); T = clock(); moebius_mu(mu.data(), n);
    printf("moebius:             %d  %0.3lf s  OK  %s\n", n, sec_since(T), message_log_log);
    vector<int> s0(n); T = clock(); divisor_sigma0(s0.data(), n);
    printf("divisors:            %d  %0.3lf s  OK  %s\n", n, sec_since(T), message_log);
    vector<ll> s1(n); T = clock(); divisor_sigma1(s1.data(), n);
    printf("sigma:               %d  %0.3lf s  OK  %s\n", n, sec_since(T), message_log);
	auto _mu = [&](int n){ return mu[n]; };
	auto _amu = [&](int n){ return abs(mu[n]); };
	auto _s1 = [&](int n){ return (int)s1[n]; };
	auto _id = [&](int n){ return n; };
	auto _c1 = [&](int n){ return 1; };

	vector<int> mu_mult(n); // mu = 1^-1
	T = clock(); dirichlet_inverse_multiplicative(mu_mult, _c1, n, p.data(), m);
    printf("moebius mult:        %d  %0.3lf s  %s  %s\n", n, sec_since(T), (mu_mult == mu) ? "OK" : "ERR", message_log_log);

	vector<int> ll_mult(n); // liouville_lambda = |mu|^-1
	T = clock(); dirichlet_inverse_multiplicative(ll_mult, _amu, n, p.data(), m);
    printf("liouville mult:      %d  %0.3lf s  OK  %s\n", n, sec_since(T), message_log_log);

	vector<int> phi_mult(n); // phi = Id * mu
	T = clock(); dirichlet_convolution_multiplicative(phi_mult, _id, _mu, n, p.data(), m);
    printf("totient mult:        %d  %0.3lf s  %s  %s\n", n, sec_since(T), (phi_mult == phi) ? "OK" : "ERR", message_log_log);

	vector<int> s0_mult(n); // d = 1 * 1
	T = clock(); dirichlet_convolution_multiplicative(s0_mult, _c1, _c1, n, p.data(), m);
    printf("divisors mult:       %d  %0.3lf s  %s  %s\n", n, sec_since(T), (s0_mult == s0) ? "OK" : "ERR", message_log_log_vs);

	vector<ll> s1_mult(n); // s_k = Id_k * 1
	T = clock(); dirichlet_convolution_multiplicative(s1_mult, _id, _c1, n, p.data(), m);
    printf("sigma mult:          %d  %0.3lf s  %s  %s\n", n, sec_since(T), (s1_mult == s1) ? "OK" : "ERR", message_log_log_vs);

	vector<int> id(n);
	T = clock(); dirichlet_convolution(id, _s1, _mu, n);
    printf("convolution gen:     %d  %0.3lf s  OK  %s\n", n, sec_since(T), message_log);
	T = clock(); dirichlet_convolution_multiplicative(id, _s1, _mu, n, p.data(), m);
    printf("convolution mult:    %d  %0.3lf s  OK  %s\n", n, sec_since(T), message_log_log);
	T = clock(); dirichlet_convolution_completely_multiplicative(id, _s1, _mu, n, pf.data());
    printf("convolution tot:     %d  %0.3lf s  OK  %s\n", n, sec_since(T), message_n);

	vector<int> i1(n);
	T = clock(); dirichlet_inverse(i1, _mu, n);
    printf("inverse gen:         %d  %0.3lf s  OK  %s\n", n, sec_since(T), message_log);
	T = clock(); dirichlet_inverse_multiplicative(i1, _mu, n, p.data(), m);
    printf("inverse mult:        %d  %0.3lf s  OK  %s\n", n, sec_since(T), message_log_log);
	T = clock(); dirichlet_inverse_completely_multiplicative(i1, _mu, n, pf.data());
    printf("inverse tot:         %d  %0.3lf s  OK  %s\n", n, sec_since(T), message_n);

	vector<int> tr(n);
	T = clock(); moebius_transform(tr, _id, n);
    printf("moebius trans:       %d  %0.3lf s  %s  %s\n", n, sec_since(T), (tr == phi) ? "OK" : "ERR", message_log);
	T = clock(); moebius_transform_multiplicative(tr, _id, n, p.data(), m);
    printf("moebius trans mult:  %d  %0.3lf s  %s  %s\n", n, sec_since(T), (tr == phi) ? "OK" : "ERR", message_log_log);
	
	cout << endl;
}
