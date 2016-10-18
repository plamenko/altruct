#include <iostream>
#include <iomanip>

#include "structure/math/modulo.h"
#include "algorithm/math/recurrence.h"

using namespace std;
using namespace altruct::math;

void fibonacci_mod_m_sample() {
	cout << "=== fibonacci_mod_m_sample ===" << endl;
	cout << "f[0] = 0, f[1] = 1" << endl;
	cout << "f[n] = (f[n-1] + f[n-2]) % 1000" << endl;

	typedef modulo<int, 1000> mod;
	
	// altruct `fibonacci` method with `modulo` structure
	cout << "altruct: ";
	for (int n = 0; n <= 30; n++) {
		cout << fibonacci<mod>(n).v << " ";
	}
	cout << endl;

	// manual
	cout << "manual:  ";
	int fp = 0, fn = 1;
	cout << fp << " " << fn << " ";
	for (int n = 2; n <= 30; n++) {
		int f = (fp + fn) % 1000;
		cout << f << " ";
		fp = fn, fn = f;
	}
	cout << endl;
	
	cout << endl;
}

void linear_recurrence_sample() {
	cout << "=== linear_recurrence_sample ===" << endl;
	cout << "f[0] = 0.1, f[1] = -0.2, f[2] = 0.3" << endl;
	cout << "f[n] = 2.5 f[n-1] + 3.0 f[n-2] - 5.0 f[n-3]" << endl;

	// altruct `linear_recurrence` method with `double`
	cout << "altruct: ";
	for (int n = 0; n <= 15; n++) {
		cout << setprecision(3) << fixed << linear_recurrence<double, double>({ 2.5, 3.0, -5.0 }, { 0.1, -0.2, 0.3 }, n) << " ";
	}
	cout << endl;

	// manual
	cout << "manual:  ";
	double fn_3 = 0.1, fn_2 = -0.2, fn_1 = 0.3;
	cout << fn_3 << " " << fn_2 << " " << fn_1 << " ";
	for (int n = 3; n <= 15; n++) {
		double fn = 2.5 * fn_1 + 3 * fn_2 - 5 * fn_3;
		cout << setprecision(3) << fixed << fn << " ";
		fn_3 = fn_2, fn_2 = fn_1, fn_1 = fn;
	}
	cout << endl;
	cout << endl;
}

void linear_recurrence_sample2() {
	cout << "=== linear_recurrence_sample 2 ===" << endl;

	// manual calculation of linear recurrence element by element
	cout << "manual:      ";
	int a_n3 = 1, a_n2 = 4, a_n1 = 7; // a[n-3], a[n-2], a[n-1]
	cout << a_n3 << " " << a_n2 << " " << a_n1 << " ";
	for (int n = 3; n <= 15; n++) {
		int a_n = (2 * a_n1 - 1 * a_n2 + 1 * a_n3) % 1009; // a[n]
		cout << a_n << " ";
		a_n3 = a_n2, a_n2 = a_n1, a_n1 = a_n;
	}
	cout << endl;

	// using exponentiation modulo characteristic polynomial
	// to calculate the n-th element of a linear recurrence
	cout << "polynomial:  ";
	typedef modulo<int, 1009> mod;
	typedef polynom<mod> poly;
	typedef moduloX<poly> polymod;
	poly init = { 1, 4, 7 };
	poly p = { -1, +1, -2, 1 }; // 0 == a[n] - 2 a[n - 1] - a[n - 2] + a[n - 3]
	poly x = { 0, 1 };      // 0 + 1*x
	for (int n = 0; n <= 15; n++) {
		polymod xn = powT(polymod(x, p), n); // x^n % p(x)
		mod r = 0;
		for (int i = 0; i < p.size(); i++) {
			r += init[i] * xn.v[i];
		}
		cout << r.v << " ";
	}
	cout << endl;

	// using altruct built-in functions
	cout << "altruct:     ";
	for (int n = 0; n <= 15; n++) {
		cout << linear_recurrence<mod, mod>({ 2, -1, +1 }, { 1, 4, 7 }, n).v << " ";
	}
	cout << endl;

	cout << endl;
}
