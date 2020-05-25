#include <chrono>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <vector>

#include "altruct/chrono/chrono.h"
#include "altruct/algorithm/math/divisor_sums.h"
#include "altruct/algorithm/math/prime_counting.h"
#include "altruct/structure/math/modulo.h"
#include "altruct/structure/math/prime_holder.h"

using namespace std;
using namespace altruct::chrono;
using namespace altruct::math;

using clk = chrono::high_resolution_clock;

typedef modulo<int, 1000000007, modulo_storage::CONSTANT> mod;

void multiplicative_sum_sample() {
    cout << "=== multiplicative_sum_sample ===" << endl;
    for (int h = 6; h <= 6; h++) {
        auto T0 = clk::now();
        cout << "10^" << h << endl;

        const int64_t n = powT<int64_t>(10, h);
        int q = isqrt(n);

        auto T = clk::now();
        prime_holder prim(q + 1);
        prim.p();
        cout << "p: " << since(T) << endl;

        T = clk::now();
        auto pi = prime_power_sum_sqrt(0, n, int64_t(1));
        cout << "pi: " << since(T) << endl;

        T = clk::now();
        auto ps1 = prime_power_sum_sqrt(1, n, int64_t(1));
        cout << "ps1: " << since(T) << endl;

        T = clk::now();
        auto ps2 = prime_power_sum_sqrt(2, n, int64_t(1));
        cout << "ps2: " << since(T) << endl;

        {
            T = clk::now();
            // f(n) = ds0(n), f(p^e) = e + 1, s1(n) = pi(n) * 2
            auto f = [&](int64_t f_pe1, int p, int e) { return e + 1; };
            auto s1 = [&](int64_t n) { return pi[n] * 2; };
            auto ans = multiplicative_sum<int64_t>(s1, f, n, prim.p().data(), (int)prim.p().size());
            cout << "mul_sum ds0(k): " << ans << " " << since(T) << endl;
        }

        {
            T = clk::now();
            // f(n) = ds0(n^2), f(p^e) = 2 e + 1, s1(n) = pi(n) * 3
            auto f = [&](int64_t f_pe1, int p, int e) { return 2 * e + 1; };
            auto s1 = [&](int64_t n) { return pi[n] * 3; };
            auto ans = multiplicative_sum<int64_t>(s1, f, n, prim.p().data(), (int)prim.p().size());
            cout << "mul_sum ds0(k^2): " << ans << " " << since(T) << endl;
        }

        {
            T = clk::now();
            // f(n) = ds1(n), f(p^e) = (p^(e+1)-1)/(p-1), s1(n) = ps1(n) + pi(n)
            auto f = [&](int64_t f_pe1, int p, int e) { return f_pe1 * p + 1; };
            auto s1 = [&](int64_t n) { return ps1[n] + pi[n]; };
            auto ans = multiplicative_sum<int64_t>(s1, f, n, prim.p().data(), (int)prim.p().size());
            cout << "mul_sum ds1(k): " << ans << " " << since(T) << endl;
        }

        {
            T = clk::now();
            // f(n) = ds2(n), f(p^e) = (p^(2*e+2)-1)/(p^2-1), s1(n) = ps2(n) + pi(n)
            auto f = [&](mod f_pe1, int p, int e) { return f_pe1 * p * p + 1; };
            auto s1 = [&](int64_t n) { return castOf<mod>(ps2[n] + pi[n]); };
            auto ans = multiplicative_sum<mod>(s1, f, n, prim.p().data(), (int)prim.p().size());
            cout << "mul_sum ds2(k): " << ans.v << " " << since(T) << endl;
        }

        {
            T = clk::now();
            // f(n) = phi(n), f(p^e) = (p-1)p^(e-1), s1(n) = ps1(n) - pi(n)
            auto f = [&](int64_t f_pe1, int p, int e) { return (e == 1) ? p - 1 : f_pe1 * p; };
            auto s1 = [&](int64_t n) { return ps1[n] - pi[n]; };
            auto ans = multiplicative_sum<int64_t>(s1, f, n, prim.p().data(), (int)prim.p().size());
            cout << "mul_sum phi(k): " << ans << " " << since(T) << endl;
        }

        {
            T = clk::now();
            // f(n) = mu(n), f(p^e) = -[e == 1], s1(n) = -pi(n)
            auto f = [&](int64_t f_pe1, int p, int e) { return (e == 1) ? -1 : 0; };
            auto s1 = [&](int64_t n) { return -pi[n]; };
            auto ans = multiplicative_sum<int64_t>(s1, f, n, prim.p().data(), (int)prim.p().size());
            cout << "mul_sum mu(k): " << ans << " " << since(T) << endl;
        }

        {
            T = clk::now();
            // f(n) = rad(n), f(p^e) = p, s1(n) = ps1(n)
            auto f = [&](int64_t f_pe1, int p, int e) { return p; };
            auto s1 = [&](int64_t n) { return ps1[n]; };
            auto ans = multiplicative_sum<int64_t>(s1, f, n, prim.p().data(), (int)prim.p().size());
            cout << "mul_sum rad(k): " << ans << " " << since(T) << endl;
        }

        {
            T = clk::now();
            // f(n) = lambda(n), f(p^e) = (-1)^e, s1(n) = -pi(n)
            auto f = [&](int64_t f_pe1, int p, int e) { return (e % 2 == 0) ? +1 : -1; };
            auto s1 = [&](int64_t n) { return -pi[n]; };
            auto ans = multiplicative_sum<int64_t>(s1, f, n, prim.p().data(), (int)prim.p().size());
            cout << "mul_sum lambda(k): " << ans << " " << since(T) << endl;
        }

        {
            T = clk::now();
            // f(n) = 2^nu(n), f(p^e) = 2, s1(n) = pi(n) * 2
            auto f = [&](int64_t f_pe1, int p, int e) { return 2; };
            auto s1 = [&](int64_t n) { return pi[n] * 2; };
            auto ans = multiplicative_sum<int64_t>(s1, f, n, prim.p().data(), (int)prim.p().size());
            cout << "mul_sum 2^nu(k): " << ans << " " << since(T) << endl;
        }

        {
            T = clk::now();
            // f(n) = 2^omega(n), f(p^e) = 2^e, s1(n) = pi(n) * 2
            auto f = [&](int64_t f_pe1, int p, int e) { return 1 << e; };
            auto s1 = [&](int64_t n) { return pi[n] * 2; };
            auto ans = multiplicative_sum<int64_t>(s1, f, n, prim.p().data(), (int)prim.p().size());
            cout << "mul_sum 2^omega(k): " << ans << " " << since(T) << endl;
        }

        {
            T = clk::now();
            // f(n) = n^2, f(p^e) = p^(2*e), s1(n) = pi(n) * 2
            auto f = [&](mod f_pe1, int p, int e) { return f_pe1 * p * p; };
            auto s1 = [&](int64_t n) { return castOf<mod>(ps2[n]); };
            auto ans = multiplicative_sum<mod>(s1, f, n, prim.p().data(), (int)prim.p().size());
            cout << "mul_sum k^2: " << ans.v << " " << since(T) << endl;
        }

        cout << endl;
    }
}
