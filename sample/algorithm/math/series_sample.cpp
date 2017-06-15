#include <iostream>
#include <iomanip>

#include "altruct/algorithm/math/ranges.h"
#include "altruct/structure/math/modulo.h"
#include "altruct/structure/math/series.h"

using namespace std;
using namespace altruct::math;

void series_simple_counting_sample() {
    cout << "=== series_simple_counting_sample ===" << endl;

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

void series_combinatoric_sample() {
    cout << "=== series_combinatoric_sample ===" << endl;

    typedef modulo<int, 1000000006> mode;
    typedef modulo<int, 1000000007> mod;

    const int K = 5;
    const int N = 15; // we are interested in `n` up to 15
    typedef series<mod, N + 1> ser;

    // fact[n] = n! (mod M)
    vector<mod> fact(N + 1);
    factorials(fact.begin(), fact.end());

    // Sum[Binomial[n, k] * x^n / n!, n] = e^x x^k / k!
    cout << "Binomial[n, k]" << endl;
    for (int k = 1; k <= K; k++) {
        ser egf_bin_k = ser::exp(1) * powT(ser{ 0, 1 }, k) / fact[k];
        auto bin_k = egf_bin_k.make_ordinary();
        for (int n = 0; n <= N; n++) {
            cout << bin_k[n].v << " ";
        }
        cout << endl;
    }
    cout << endl;

    // Sum[StirlingS1[n, k] * x^n / n!, n] = (Log[1 + x])^k / k!
    cout << "StirlingS1[n, k]" << endl;
    for (int k = 1; k <= K; k++) {
        ser egf_s1_k = powT(ser{ 1, 1 }.ln(), k) / fact[k];
        auto s1_k = egf_s1_k.make_ordinary();
        for (int n = 0; n <= N; n++) {
            auto s1_n_k = ((n - k) % 2 == 1) ? -s1_k[n] : s1_k[n];
            cout << s1_n_k.v << " ";
        }
        cout << endl;
    }
    cout << endl;

    // Sum[StirlingS2[n, k] * x^n / n!, n] = (e^x - 1)^k / k!
    cout << "StirlingS2[n, k]" << endl;
    for (int k = 1; k <= K; k++) {
        ser egf_s2_k = powT(ser::exp(1) - ser{1}, k) / fact[k];
        auto s2_k = egf_s2_k.make_ordinary();
        for (int n = 0; n <= N; n++) {
            cout << s2_k[n].v << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "pe553[n, k]" << endl;
    // t[n] = 2^2^n (mod M)
    auto t = [](int n){ return powT(mod(2), powT(mode(2), n).v); };
    // egf_t(x) = Sum[t[n] * x^n / n!, n]
    // egf_q(x) = e^-x * egf_t(x)
    // egf_f(x) = ln(egf_q(x)) + C, egf_f(0) = 0
    auto egf_f = (ser::exp(-1) * ser::of(t).make_exponential()).ln(0);
    for (int k = 1; k <= K; k++) {
        // egf_r_k(x) = f(x)^k / k!
        auto egf_r_k = powT(egf_f, k) / fact[k];
        // egf_c_k(x) = e^x * r_k(x)
        ser egf_c_k = ser::exp(1) * egf_r_k;
        auto c_k = egf_c_k.make_ordinary();
        for (int n = 0; n <= N; n++) {
            cout << c_k[n].v << " ";
        }
        cout << endl;
    }
    cout << endl;
}
