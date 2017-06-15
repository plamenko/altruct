#pragma once

#include "altruct/algorithm/math/base.h"

#include <array>

namespace altruct {
namespace math {

/**
 * A vector of type T and fixed size N.
 */
template<typename T, int N>
class vectorNd {
public:
    std::array<T, N> a;

    vectorNd() { a.fill(0); }
    vectorNd(const std::array<T, N>& rhs) : a(rhs) { }
    vectorNd(const vectorNd& rhs) : a(rhs.a) { }
    vectorNd(const T& a0) { a.fill(a0); }
    // construct from int, but only if T is not integral to avoid constructor clashing
    template <typename I = T, typename = std::enable_if_t<!std::is_integral<I>::value>>
    vectorNd(int a0) { a.fill(a0); } // to allow constructing from 0 and 1
    vectorNd(std::initializer_list<T> list) { assign(list.begin()); }

    int size() const { return N; }
    const T& operator [] (int index) const { return a[index]; }
    T& operator [] (int index) { return a[index]; }

    void assign(const T* rhs) {
        for (int i = 0; i < N; i++) {
            a[i] = rhs[i];
        }
    }

    int cmp(const vectorNd& rhs) const {
        for (int i = 0; i < N; i++) {
            if (a[i] < rhs.a[i]) return -1;
            if (a[i] > rhs.a[i]) return +1;
        }
        return 0;
    }

    bool operator == (const vectorNd& rhs) const { return cmp(rhs) == 0; }
    bool operator != (const vectorNd& rhs) const { return cmp(rhs) != 0; }
    bool operator <  (const vectorNd& rhs) const { return cmp(rhs) < 0; }
    bool operator >  (const vectorNd& rhs) const { return cmp(rhs) > 0; }
    bool operator <= (const vectorNd& rhs) const { return cmp(rhs) <= 0; }
    bool operator >= (const vectorNd& rhs) const { return cmp(rhs) >= 0; }

    vectorNd  operator +  (const vectorNd& rhs) const { vectorNd t(*this); t += rhs; return t; }
    vectorNd  operator -  (const vectorNd& rhs) const { vectorNd t(*this); t -= rhs; return t; }
    vectorNd  operator -  ()                    const { vectorNd t; for (int i = 0; i < N; i++) t.a[i] = -a[i]; return t; }
    vectorNd  operator *  (const vectorNd& rhs) const { vectorNd t(*this); t *= rhs; return t; }
    vectorNd  operator /  (const vectorNd& rhs) const { vectorNd t(*this); t /= rhs; return t; }
    vectorNd  operator %  (const vectorNd& rhs) const { vectorNd t(*this); t %= rhs; return t; }

    vectorNd  operator +  (const T& rhs) const { vectorNd t(*this); t += rhs; return t; }
    vectorNd  operator -  (const T& rhs) const { vectorNd t(*this); t -= rhs; return t; }
    vectorNd  operator *  (const T& rhs) const { vectorNd t(*this); t *= rhs; return t; }
    vectorNd  operator /  (const T& rhs) const { vectorNd t(*this); t /= rhs; return t; }
    vectorNd  operator %  (const T& rhs) const { vectorNd t(*this); t %= rhs; return t; }

    vectorNd& operator += (const vectorNd& rhs) { for (int i = 0; i < N; i++) a[i] += rhs.a[i]; return *this; }
    vectorNd& operator -= (const vectorNd& rhs) { for (int i = 0; i < N; i++) a[i] -= rhs.a[i]; return *this; }
    vectorNd& operator *= (const vectorNd& rhs) { for (int i = 0; i < N; i++) a[i] *= rhs.a[i]; return *this; }
    vectorNd& operator /= (const vectorNd& rhs) { for (int i = 0; i < N; i++) a[i] /= rhs.a[i]; return *this; }
    vectorNd& operator %= (const vectorNd& rhs) { for (int i = 0; i < N; i++) a[i] %= rhs.a[i]; return *this; }

    vectorNd& operator += (const T& rhs) { for (int i = 0; i < N; i++) a[i] += rhs; return *this; }
    vectorNd& operator -= (const T& rhs) { for (int i = 0; i < N; i++) a[i] -= rhs; return *this; }
    vectorNd& operator *= (const T& rhs) { for (int i = 0; i < N; i++) a[i] *= rhs; return *this; }
    vectorNd& operator /= (const T& rhs) { for (int i = 0; i < N; i++) a[i] /= rhs; return *this; }
    vectorNd& operator %= (const T& rhs) { for (int i = 0; i < N; i++) a[i] %= rhs; return *this; }

    T abs2() const { T r = 0; for (int i = 0; i < N; i++) r += a[i] * a[i]; return r; }
};

template<typename T, int N>
struct identityT<vectorNd<T, N>> {
    static vectorNd<T, N> of(const vectorNd<T, N>& x) {
        vectorNd<T, N> e1;
        for (int i = 0; i < N; i++) {
            e1[i] = identityT<T>::of(x[i]);
        }
        return e1;
    }
};

template<typename T, int N>
struct zeroT<vectorNd<T, N>> {
    static vectorNd<T, N> of(const vectorNd<T, N>& x) {
        vectorNd<T, N> e0;
        for (int i = 0; i < N; i++) {
            e0[i] = zeroT<T>::of(x[i]);
        }
        return e0;
    }
};

} // math
} // altruct
