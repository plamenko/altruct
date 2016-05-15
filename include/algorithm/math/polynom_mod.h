#pragma once

#include "gmp_helpers.h"
#include "modulos.h"
#include "structure/math/modulo.h"
#include "structure/math/polynom.h"

namespace altruct {
namespace math {

/**
 * polynom<modulo<int>> specialization
 *
 * Important: not thread-safe as mod::M gets modified.
 */
template<int ID>
struct altruct::math::polynom_mul<modulo<int, ID, true>> {
	typedef modulo<int, ID, true> mod;
	static int threshold() { return 17000; /* karatsuba is better for smaller size */ }
	static void impl(polynom<mod> &pr, const polynom<mod> &p1, const polynom<mod> &p2, int lr = -1) {
		int l1 = p1.deg(), l2 = p2.deg(); if (lr < 0) lr = l1 + l2;
		// separate convolutions
		int M = mod::M; // store M
		const int P1 = 1012924417, R1 = 198, O1 = 1 << 21;
		const int P2 = 1004535809, R2 = 4172, O2 = 1 << 21;
		const int P3 = 985661441, R3 = 210, O3 = 1 << 22;
		mod::M = P1; vector<mod> r1 = convolution(p1.c.cbegin(), p1.c.cend(), p2.c.cbegin(), p2.c.cend(), mod(R1), O1);
		mod::M = P2; vector<mod> r2 = convolution(p1.c.cbegin(), p1.c.cend(), p2.c.cbegin(), p2.c.cend(), mod(R2), O2);
		mod::M = P3; vector<mod> r3 = convolution(p1.c.cbegin(), p1.c.cend(), p2.c.cbegin(), p2.c.cend(), mod(R3), O3);
		mod::M = M; // restore M
		// combine with CRT
		pr.c.resize(lr + 1);
		for (int i = 0; i <= lr; i++) {
			mpz zr = r1[i].v, zm = P1;
			chinese_remainder<mpz>(&zr, &zm, r2[i].v, P2);
			chinese_remainder<mpz>(&zr, &zm, r3[i].v, P3);
			zr %= mod::M;
			pr.c[i].v = zr.get_si();
		}
	}
};

} // math
} // altruct
