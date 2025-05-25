#pragma once

#include "altruct/algorithm/math/ranges.h"

#include <algorithm>
#include <vector>

namespace altruct {
namespace math {

template<typename T>
class factorial_holder {
private:
    int sz;                 // upper bound (exclusive)
    std::vector<T> _fact;   // factorials
    std::vector<T> _ifact;  // inverse factorials
    std::vector<T> _inv;    // inverses 

public:
    factorial_holder(int sz, T id = T(1)) : sz(sz) {
        _fact = make_factorials<T>(sz, id);
        _ifact = make_inv_factorials<T>(sz, _fact[sz - 1], sz - 1);
        _inv = _ifact; inverses_from_ifact(_inv.begin(), _inv.end(), id);
    }

    int size() { return sz; }

    std::vector<T>& fact() { return _fact; }
    std::vector<T>& ifact() { return _ifact; }
    std::vector<T>& inv() { return _inv; }

    T fact(size_t k) { return fact().at(k); }
    T ifact(size_t k) { return ifact().at(k); }
    T inv(size_t k) { return inv().at(k); }
    T bin(size_t n, size_t k) {
        if (k < 0 || k > n) return zeroOf(_fact[0]);
        return _fact[n] * _ifact[n - k] * _ifact[k];
    }
};

} // math
} // altruct
