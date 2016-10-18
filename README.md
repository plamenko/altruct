Altruct is a collection of algorithms and data structures, hence the name: **Al**go + S**truct**.
It is written in C++ and predominately uses templates. This is both a curse and a blessing.
Templates allow powerful and fast abstractions at a cost of higher compilation times and harder debugging.


## Motivation
* Reusability. It is tedious to implement the same things over and over again.
 When you read a problem on a programming contest you may immediately have an idea of how to solve it.
 But before that, you need to implement several tricky algorithms and/or structures. Instead of focusing on
 a bigger picture and how to combine those building blocks, you need to spend time writing blocks themselves.
* Reliability. Related to the above point, it is a big relief to know that your building blocks are correct.
 If your solution is not working correctly, chances are much higher that the bug is in the top-level code
 and/or logical reasoning than an edge case in the implementation of the building blocks used.
* Apart from the practical reasons above, I also wanted to deepen my knowledge in computer science and mathematics.
* I also wanted to deepen my knowledge in software engineering by designing a library.
 There are many challenges involved in properly organizing the code, properly testing it,
 making sure it is robust and portable, and finally making it open source.

 
## Unit test coverage and samples
Altruct has extensive unit tests coverage powered by `googletest`.
Those unit tests are also its best usage examples.

### Experimental
Nearly every bit of the code has a unit test that covers it in multiple scenarios.
The code that doesn't yet have a proper unit test is included in the `experinemntal` directory.
Much like the rest of the library, tode in `experimental` has been successfully used
in multiple programming contest problems and chances of having a bug are not big.
However, until the proper unit tests get added, it will remain experimental.


## Supported plaftorms
* Altruct has been compiled on Windows, Linux and Mac OSX.
* Altruct has been compiled with MSVS, GCC and Clang.
* Altruct comes with a MSVS solution and a BUCK file.
* Code is written in endian independent way.

I occasionally check that Altruct still compiles on all of the above,
however I am almost exclusively using MSVS on Windows so there is a
posibility of something being temporarily broken for other compilers.
In the future I may set up some continuous integration tests running.


## Example

Here is a naive implementation for calculating the n-th term of a linear recurrence of degree L (in modulo M arithmetic) that works in O(L * n).

```
	int a_k3 = 1, a_k2 = 4, a_k1 = 7; // a[k-3], a[k-2], a[k-1]
	for (int k = 3; k <= n; k++) {
		int a_k = (2 * a_k1 - 1 * a_k2 + 1 * a_k3) % 1009;
		a_k3 = a_k2, a_k2 = a_k1, a_k1 = a_k;
	}
	cout << a_k1 << endl;
```

Here is how to do this in O(L^2 log n) by composing the mathematical structures Altruct provides.
The algorithm calculates the fast exponentiation of polynom `x` to the `n-th` power in O(log n) steps,
where each operation is performed modulo characteristic polynomial in O(L^2). There is a similar
algorithm that uses matrix eponentiation but that would cost O(L^3) instead of O(L^2) per step.

```
	typedef modulo<int, 1009> mod;
	typedef polynom<mod> poly;
	typedef moduloX<poly> polymod;
	// the initial values
	poly init = { 1, 4, 7 };
	// the characteristic polynomial of our recurrence
	// 0 == a[n] - 2 a[n - 1] - a[n - 2] + a[n - 3]
	poly p = { -1, +1, -2, 1 };
	poly x = { 0, 1 }; // 0 + 1*x
	polymod xn = powT(polymod(x, p), n); // x^n % p(x)
	mod r = 0;
	for (int i = 0; i < p.size(); i++) {
		r += init[i] * xn.v[i];
	}
	cout << r.v << endl;
```

But of course, the aforementioned algorithm is already implemented as `linear_recurrence`,
so for that you can just use the Altruct function and not have to implement it by yourself.

```
    typedef modulo<int, 1009> mod;
	cout << linear_recurrence<mod, mod>({ 2, -1, +1 }, { 1, 4, 7 }, n).v << endl;
```

Observe how in the above excample, the actual type used is `modulo<polynom<modulo<int>>>`.
All of the math structures Altruct provides can be composed in such a way where it makes sense.

The same implementation of Fast Fourier Transform and related convolutions can be used both
on complex numbers and on modulo numbers. Or on any other field provided that you know how to
calculate the n-th root of unity in it.

See below for a comprehensive listing of the library content.


## Library structure

There are several levels of source code branching:
* source - This is the main body of the library.
* experimental - Same as `source`, but without proper unit test coverage.
* sample - A limited set of examples as `test` already shows the usage.
* test - Unit tests for the rest of the library.
* test_util - Utility code that comes handy in testing, but is not the main body of the library.

Further subdivision is about headers files and source files. All of the above use one or both of:
* include - A list of header files containing templates and declarations.
* src - A list of C++ source files.
* Further subdivision is about algorithms and data structures and the respective subdomains.

Finally there are build files, source control files and other repository metadata.


## Implemented algorithms and structures
* Algorithm:
  * Graph:
    * Dinic blocking flow (max flow)
    * Heavy-light Decomposition of a tree
    * Iterative DFS
    * Lowest Common Ancestor
  * Hash:
    * Polynomial hash (useful for string related problems)
    * std::tuple hash
  * Math:
    * Base:
      * IdentityT, ZeroT,
      * absT, minT, maxT,
      * powT, sqT, sqrtT, cbT, cbrtT
      * gcd, gcd_ex, gcd_max, lcm
      * div_floor, div_ceil, div_round, multiple
    * Bits:
      * ilog2, lzc, tzc, bit_cnt1, bit_reverse,
      * or_down, xor_down, gray_to_bin, bin_to_gray
      * hi_bit, lo_bit, is_pow2, next_pow2
    * Combinatorial primes:
      * factorial/binomial/multinomial prime exponent
    * Combinatorics:
      * next_partition, next_combination, nth_permutation
    * Convolutions:
      * Slow and/or/xor/max/cyclic convolution (O(n^2) implementation)
      * Fast and/or/xor/max convolution (FFT-like O(n log n) implementation)
      * Fast Walsh-Hadamard transform, fast arithmetic transform
      * Fast Fourier transform
      * Fast ordinary/cyclic convolution
    * Fractions:
      * Farey sequence
    * Modulos:
      * Chinese Remainder Theorem, Garner algorithm
      * Jacobi symbol
      * Cipolla algorithm (modular square root)
      * Hensel lift (square root modulo prime power)
      * Primitive root, primitive root of unity
      * k-th root, k-th roots of unity
      * Factorial/binomial modulo prime power
    * Polynoms:
      * Monotonic search
      * Zeros
      * Discrete integral
    * Primes:
      * Precompute for a range (1 to n, or segmented):
        * Primes, Prime-Pi, Euler-Phi (Totient), Moebius-Mu, Divisor-Sigma, Prime-Factor
        * Moebius transform, factorization 
      * Compute from a given factorization:
        * Divisors, Euler-Phi, CarmichaelLambda
      * Pollard-Rho factorization, Miller-Rabin primality test
      * Integer digits for a base
    * Ranges:
      * Arithmetic progression
      * Powers, Factorials, Inverse
      * Negation, Alternate sign
      * Accumulate, Differences
    * Recurrences:
      * Linear recurrence
      * Fibonacci, Lucas-L, Lucas-U, Lucas-V
      * Bernoulli numbers
      * Berlekamp-Massey algorithm (finds the characteristic polynomial of a linear recurrence)
    * Reduce:
      * Sum, Product
      * Min, Max, Mex (Grundy minimal excludant)
    * Sums:
      * Sum of f
      * Sum of powers
      * Sum of floors in O(sqrt n)
      * Sum of multiplicative f
      * Mertens function (Sum of Moebius-Mu)
      * Sum of primes
  * Random:
    * Mersene Twister
    * XorShift
  * Search:
    * Binary search
    * Knuth–Morris–Pratt search (string search)
* I/O:
  * Fast I/O
  * Reader: File / Stream / String / Buffered / Simple()
  * Writer: File / Stream / String / Buffered / Simple
	
* Structure:
  * Container:
    * Aho-Corasick trie
    * Palindrome tree
    * Siffix array
    * Bit vector
    * Interval tree (range update/get with lazy propagation)
    * Segment tree (single element update, range get)
    * Low-High map, Sqrt map (efficient storage for a certain sparse data)
    * Range Minimum Query structure (O(n) build, O(1) query)
  * Graph:
    * Disjoint set
  * Math:
    * complex
    * double_int
    * fraction
    * galois_field_2
    * matrix
    * modulo
    * nimber
    * permuation (cycle/transposition/array notations)
    * polynom
    * quadratic
    * series (a formal power series)
    * vector2d, vector3d, vectorNd
    * fenwick_tree (a.k.a. "logaritamska struktura" in Croatia)
    * prime_holder (keeps a range of primes and related functions)
    * root_wrapper (wraps the root powers used in FFT)

