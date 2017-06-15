#include <chrono>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <vector>
#include <thread>

#include "altruct/chrono/chrono.h"
#include "altruct/structure/math/modulo.h"
#include "altruct/structure/math/prime_holder.h"
#include "altruct/io/iostream_overloads.h"

using namespace std;
using namespace std::chrono;
using namespace altruct::math;
using namespace altruct::chrono;

void test_sample() {
    typedef chrono::high_resolution_clock clk;

    typedef modulo<int, 1000000007, modulo_storage::CONSTANT> mod;
    mod r = 0;
    auto T0 = clk::now();
    for (int i = 1; i < 100000000; i++) {
        r *= i;
    }
    cout << r << " " << since(T0) << " s" << endl;

    typedef moduloX<int> modx;
    prime_holder prim(100);
    printf("primes: %d\n", prim.primes());

    for (int p : prim.p()) {
        modx a(100, p);
        modx r(0, p);
        printf("%d: ", p);
        for (int k = 0; k < p; k++) {
            r += powT(a, k);
            if (k > 0 && r.v == 1) break;
            printf("%d ", r.v);
        }
        printf("\n");
    }
    cout << flush;
}
