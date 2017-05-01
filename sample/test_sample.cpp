#include <cstdio>
#include <iostream>
#include <iomanip>
#include <vector>

#include "structure/math/modulo.h"
#include "structure/math/prime_holder.h"
#include "io/iostream_overloads.h"

using namespace std;
using namespace altruct::math;

//       0
//   1       2
// 3   4   5   6

template<typename T>
class binary_heap {
public:
  std::vector<T> v;

  binary_heap() {}
  binary_heap(size_t sz) : v(sz) {}
  binary_heap(std::vector<T>& v) : v(v) { rebuild(); }
  template<typename It>
  binary_heap(It begin, It end) : v(begin, end) { rebuild(); }

  size_t size() const {
    return v.size();
  }

  bool cmp(const T& lhs, const T& rhs) {
    return lhs < rhs;
  }

  void insert(const T& val) {
    v.push_back(val);
    heapify_up(v.size() - 1);
  }

  const T& front() const {
    return v.front();
  }

  void pop_front() {
    std::swap(v.front(), v.back());
    heapify_down(0, v.size());
  }

  void heapify_up(size_t index) {
    while (true) {
      size_t parent = (index - 1) / 2;
      if (!cmp(v[index], v[parent])) break;
      std::swap(v[index], v[parent]);
      index = parent;
    }
  }

  void heapify_down(size_t index, size_t sz) {
    size_t best = index;
    while (true) {
      size_t left = index * 2 + 1, right = index * 2 + 2;
      if (left < sz && cmp(v[left], v[best])) best = left;
      if (right < sz && cmp(v[right], v[best])) best = right;
      if (best == index) break;
      std::swap(v[index], v[best]);
      index = best;
    }
  }

  void rebuild() {
    for (size_t index = v.size(); index > 0; index--) {
      heapify_down(index - 1, v.size());
    }
  }

  /** converts heap to sorted array (heap order is destroyed) */
  void sort() {
    for (size_t index = size(); index > 0; index--) {
      std::swap(v[0], v[index - 1]);
      heapify_down(0, index);
    }
    reverse(v.begin(), v.end());
  }
};

void test_sample() {
  vector<int> va(10);
  for (int& a : va) a = rand() % 10;
  
  binary_heap<int> bh(va);
  bh.sort();
  
  sort(va.begin(), va.end());

  if (bh.v != va) {
    printf("ERROR\n");
    cout << va << endl;
    cout << bh.v << endl;
  }

  // typedef moduloX<int> modx;
  // prime_holder prim(100);
  // printf("primes: %d\n", prim.primes());

  // for (int p : prim.p()) {
  //   modx a(100, p);
  //   modx r(0, p);
  //   printf("%d: ", p);
  //   for (int k = 0; k < p; k++) {
  //     r += powT(a, k);
  //     if (k > 0 && r.v == 1) break;
  //     printf("%d ", r.v);
  //   }
  //   printf("\n");
  // }

}
