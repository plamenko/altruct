#include <iostream>
#include <fstream>

#include "altruct/algorithm/random/xorshift.h"

using namespace std;
using namespace altruct::random;

void xorshift_sample() {
    cout << "=== xorshift_sample ===" << endl;

    cout << "seeding with 1234 ..." << endl;
    xorshift_1024star rng(1234);
    cout << "random int [0, 2^64-1]: " << rng.next() << endl;
    cout << "random int [1000, 2000]: " << rng.next(1000, 2000) << endl;
    cout << "random int strongly uniform [1000, 2000]: " << rng.next_uniform(1000, 2000) << endl;
    cout << "random double [0, 1]: " << rng.next_0_1() << endl;

    // reseeding with the same seed produces the same result
    cout << "reseeding with 1234 ..." << endl;
    rng.seed(1234);
    cout << "random int [0, 2^64-1]: " << rng.next() << endl;
    cout << "random int [1000, 2000]: " << rng.next(1000, 2000) << endl;
    cout << "random int strongly uniform [1000, 2000]: " << rng.next_uniform(1000, 2000) << endl;
    cout << "random double [0, 1]: " << rng.next_0_1() << endl;

    cout << endl;
}
