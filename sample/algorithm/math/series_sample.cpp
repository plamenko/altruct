#include <iostream>
#include <iomanip>

#include "algorithm/math/ranges.h"
#include "structure/math/modulo.h"
#include "structure/math/series.h"

using namespace std;
using namespace altruct::math;

void series_counting_sample() {
	cout << "=== series_ogf_sample ===" << endl;

	// A generating function for the number of submultisets of {inf a, inf b, inf c}
	// in which there are an odd number of `a`s, an even number of `b`s, and any number of `c`s.
	// g(x) = x / ((1 + x) ^ 2 (1 - x) ^ 3)
	typedef series<int, 100 + 1> ser;
	auto s1 = ser{ 0, 1 } / (powT(ser{ 1, 1 }, 2) * powT(ser{ 1, -1 }, 3));
	for (int i = 0; i <= 30; i++) {
		cout << s1[i] << " ";
	}
	cout << endl;

	// An exponential generating function for the number of
	// permutations with repetition of length `n` of the set `{a, b, c}`,
	// in which there are an odd number of `a`s, an even number of `b`s, and any number of `c`s.
	// g(x) = (e^3x - e^-x) / 4
	typedef series<double, 100 + 1> serd;
	auto s2 = ((serd::exp(3) - serd::exp(-1)) / 4).make_ordinary();
	for (int i = 0; i <= 15; i++) {
		cout << int(s2[i]) << " ";
	}
	cout << endl;

	cout << endl;
}

void series_egf_sample() {
	cout << "=== series_egf_sample ===" << endl;

	typedef modulo<int, 1000000006> mode;
	typedef modulo<int, 1000000007> mod;
	typedef series<mod, 1000 + 1> ser; // we are interested in `n` up to 1000

	int n = 1000, k = 2;
	
	// fact[n] = n! (mod M)
	vector<mod> fact(n + 1);
	factorials(fact.begin(), fact.end());

	// d[n] = 2^n (mod phi(M))
	vector<mode> d(n + 1);
	powers(d.begin(), d.end(), mode(2));
	
	// t[n] = 2^2^n (mod M)
	ser t; for (int i = 0; i <= n; i++) t[i] = powT(mod(2), d[i].v);
	// egf_t(x) = Sum[t[n] * x^n / n!, n]
	auto egf_t = t.make_exponential();
	
	// egf_q(x) = e^-x * t(x)
	auto egf_q = ser::exp(-1) * egf_t;
	
	// egf_f(x) = ln(egf_q(x)) + C, egf_f(0) = 0
	auto egf_f = egf_q.ln(0);
	
	// egf_r_k(x) = f(x)^k / k!
	auto egf_r_k = powT(egf_f, k) / fact[k];

	// egf_c_k(x) = e^x * r_k(x)
	ser egf_c_k = ser::exp(1) * egf_r_k;

	// c_k[n] calculated up to 1000, but display only up to 10
	auto c_k = egf_c_k.make_ordinary();
	for (int i = 0; i <= 10; i++) {
		cout << c_k[i].v << " ";
	}
	cout << endl;
	
	cout << endl;
}
