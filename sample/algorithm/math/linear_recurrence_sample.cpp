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
