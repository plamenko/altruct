#include "altruct/chrono/chrono.h"

#include "gtest/gtest.h"

#include <string>

using namespace std;
using namespace std::chrono;

using namespace altruct::chrono;

template <class clock>
void test_chrono_empty_loop() {
    typedef duration<long long, std::pico> picoseconds;
    typedef duration<double, typename clock::period> dduration;

    // do the loop
    const int N = 100000000;
    auto t0 = clock::now();
    for (volatile int j = 0; j < N; ++j);
    auto t1 = clock::now();

    // get the clock ticks per iteration
    auto ticks_per_iter = dduration(t1 - t0) / N;
    auto ticks_per_iter2 = since(t0) / N;
    std::cout << ticks_per_iter.count() << " clock ticks per iteration\n";

    // convert to real time units
    std::cout << duration_cast<seconds>(ticks_per_iter).count() << "s per iteration\n";
    std::cout << duration_cast<milliseconds>(ticks_per_iter).count() << "ms per iteration\n";
    std::cout << duration_cast<microseconds>(ticks_per_iter).count() << "us per iteration\n";
    std::cout << duration_cast<nanoseconds>(ticks_per_iter).count() << "ns per iteration\n";
    std::cout << duration_cast<picoseconds>(ticks_per_iter).count() << "ps per iteration\n";
    std::cout << setprecision(3) << ticks_per_iter2 << "s per iteration\n";
}

TEST(chono, chrono) {
    std::cout << "\nUsing rdtsc:\n";
    test_chrono_empty_loop<altruct::chrono::rdtsc_clock<2666666666>>();

    std::cout << "\nUsing std::chrono::high_resolution_clock:\n";
    test_chrono_empty_loop<std::chrono::high_resolution_clock>();

    std::cout << "\nUsing std::chrono::system_clock:\n";
    test_chrono_empty_loop<std::chrono::system_clock>();
}
