#include <cstdio>
#include <iostream>
#include <iomanip>
#include <vector>

#include "altruct/structure/math/modulo.h"
#include "altruct/structure/math/prime_holder.h"
#include "altruct/io/iostream_overloads.h"

using namespace std;
using namespace altruct::math;


void test_sample() {
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
}
