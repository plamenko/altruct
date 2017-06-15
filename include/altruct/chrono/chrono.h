#pragma once

#include <chrono>

#ifdef _MSC_VER
#include <intrin.h>
#pragma intrinsic(__rdtsc)
#elif __GNUC__
#include <x86intrin.h>
#endif

namespace altruct {
namespace chrono {

template<long long CPU_FREQUENCY = 2666666666>
struct rdtsc_clock {
    typedef long long rep;
    typedef std::ratio<1, CPU_FREQUENCY> period;
    typedef std::chrono::duration<rep, period> duration;
    typedef std::chrono::time_point<rdtsc_clock> time_point;
    static const bool is_steady = true;
    static time_point now() { return time_point(duration(_now())); }

#ifdef _MSC_VER
    static long long _now() { return __rdtsc(); }
#elif __GNUC__
    static long long _now() { return __rdtsc(); }
#endif
};

template<typename TIME_POINT>
double since(const TIME_POINT& T0) {
    return std::chrono::duration<double>(TIME_POINT::clock::now() - T0).count();
}

}
}
