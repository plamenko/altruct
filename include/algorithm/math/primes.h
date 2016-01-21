#pragma once

#include "base.h"
#include <vector>

namespace altruct {
namespace math {

/**
 * Prime sieve of Eratosthenes up to `n`
 *
 * Performs Eratosthenes prime sieving and stores the result in the provided
 * arrays. Array `p` contains all the primes up to `n`, whereas array `q`
 * contains flags whether the number at the given index is prime or not. If
 * null pointer is passed for one of those arrays, the corresponding output
 * won't be stored. At least one of `p` or `q` needs to be non-null.
 *
 * If only `p` is specified, `p` needs to be of size at least `n`.
 * If both `p` and `q` are specified, `p` needs to be of size at least `pi(n)`.
 * `q` if specified, needs to be of size at least `n`.
 *
 * Complexity: O(n log log n)
 *
 * @param p - array to store primes up to `n`, or null if not required
 * @param q - array to store prime flags up to `n`, or null if not required
 * @param n - performs sieving for numbers up to `n` (exclusive)
 * @return - number of primes up to `n`
 */
int primes(int *p, char *q, int n);

/**
 * Prime Pi (Number of primes) up to `n`
 *
 * Stores the number of primes up to `i` for each integer `i` less than `n`.
 *
 * Complexity: O(n)
 *
 * @param pi - array to store the number of primes
 * @param n - calculates number of primes up to `n` (exclusive)
 * @param p - array of prime numbers up to `n`
 * @param m - number of prime numbers up to `n`
 */
void prime_pi(int *pi, int n, int *p, int m);

/**
 * Euler's Phi (Number of coprimes; Totient) up to `n`
 *
 * Stores the number of numbers coprime to `i` up to `i` for each integer `i`
 * less than `n`.
 *
 * Complexity: O(n log log n)
 *
 * @param phi - array to store the number of coprimes
 * @param n - calculate number of coprimes up to `n` (exclusive)
 * @param p - array of prime numbers up to `n`
 * @param m - number of prime numbers up to `n`
 */
void euler_phi(int *phi, int n, int *p, int m);

/**
 * Moebius Mu (+/-1 if squarefree, 0 otherwise) up to `n`
 *
 * Moebius Mu for each integer `i` is defined as:
 *    0 if `i` is not a square-free positive integer
 *   +1 if `i` is a square-free positive integer with an even number of prime factors
 *   -1 if `i` is a square-free positive integer with an odd number of prime factors
 *
 * Complexity: O(n log log n)
 *
 * @param mu - array to store the result
 * @param n - calculate number of coprimes up to `n` (exclusive)
 * @param p - array of prime numbers up to `n`
 * @param m - number of prime numbers up to `n`
 */
void moebius_mu(int *mu, int n, int *p, int m);

/**
 * Segmented Euler's Phi (Number of coprimes; Totient) in range `[b, e)`
 *
 * Stores the number of numbers coprime to `i` up to `i` for each integer `i`
 * in range `[b, e)`.
 * Prime numbers up to `sqrt(e)` should provided. I.e. `p[m-1] >= sqrt(e-1)`.
 *
 * Complexity: O((e - b) log log e)
 *
 * @param phi - array to store the number of coprimes
 * @param tmp - temporary array
 * @param b, e - calculate number of coprimes in range `[b, e)`
 * @param p - array of prime numbers up to `sqrt(e)`
 * @param m - number of prime numbers up to `sqrt(e)`
 */
void segmented_phi(long long *phi, long long *tmp, long long b, long long e, int *p, int m);

/**
 * Segmented Moebius Mu (+/-1 if squarefree, 0 otherwise) in range `[b, e)`
 *
 * Moebius Mu for each integer `i` is defined as:
 *    0 if `i` is not a square-free positive integer
 *   +1 if `i` is a square-free positive integer with an even number of prime factors
 *   -1 if `i` is a square-free positive integer with an odd number of prime factors
 * Prime numbers up to `sqrt(e)` should provided. I.e. `p[m-1] >= sqrt(e-1)`.
 *
 * Complexity: O((e - b) log log e)
 *
 * @param mu - array to store the result
 * @param b, e - calculate number of coprimes in range `[b, e)`
 * @param p - array of prime numbers up to `sqrt(e)`
 * @param m - number of prime numbers up to `sqrt(e)`
 */
void segmented_mu(long long *mu, long long b, long long e, int *p, int m);

/**
 * Divisor Sigma 0 (Number of divisors) up to `n`
 *
 * Complexity: O(n log n)
 *
 * @param ds0 - array to store the result
 * @param n - calculate number of divisors up to `n` (exclusive)
 */
void divisor_sigma0(int *ds0, int n);

/**
 * Divisor Sigma 1 (Sum of divisors) up to `n`
 *
 * Complexity: O(n log n)
 *
 * @param ds1 - array to store the result
 * @param n - calculate sum of divisors up to `n` (exclusive)
 */
void divisor_sigma1(long long *ds1, int n);

/**
 * Biggest prime factor for integers up to `n`
 *
 * Sotres the biggest prime factor for integers up to n to the array `bpf`.
 * Output array `bpf` needs to be of size `n`.
 *
 * Complexity: O(n log log n)
 *
 * @param bpf - array to store the biggest prime factors for integers up to `n`
 * @param n - factors numbers up to `n` (exclusive)
 * @param p - array of prime numbers up to `n`
 * @param m - number of prime numbers up to `n`
 */
void factor(int *bpf, int n, int *p, int m);

/**
 * Prime factorization of integer `n`
 *
 * Stores the prime factors and their exponents to the vector `vf`. This
 * function requires a prime factors array for integers up to `n` inclusive
 * to be provided. That array can be precalculated with the `factor` function.
 *
 * Complexity: O(log n / log log n)
 *
 * @param vf - vector to store prime factors and the corresponding exponents
 * @param n - integer to factor
 * @param pf - array of prime factors for integers up to `n` inclusive
 */
void factor_integer(std::vector<std::pair<int, int>> &vf, int n, int *pf);

/**
 * Prime factorization of the product of integers `vn`
 *
 * Stores the prime factors and their exponents to the vector `vf`. This
 * function requires a prime factors array for integers up to `n` inclusive
 * to be provided, where `n` is the largest element in `vn`. That array can
 * be precalculated with the `factor` function.
 *
 * Complexity: O(k log n / log log n)
 *
 * @param vf - vector to store prime factors and the corresponding exponents
 * @param vn - integers whose product to factor
 * @param pf - array of prime factors for integers up to `n` inclusive
 */
void factor_integer(std::vector<std::pair<int, int>> &vf, std::vector<int> &vn, int *pf);

/**
 * Divisors of the number whose factorization is `vf`
 *
 * Stores the divisors to the vector `vd`. Only divisors up to `maxd` are stored.
 *
 * @param vd - vector to store divisors
 * @param vf - prime factorization of the original number
 * @param maxd - maximum divisor to store
 * @param d - reserved
 * @param i - reserved
 */
void divisors(std::vector<long long> &vd, const std::vector<std::pair<int, int>> &vf, long long maxd = LLONG_MAX, long long d = 1, int i = 0);

} // math
} // altruct
