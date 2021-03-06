#pragma once

#include "altruct/structure/math/polynom.h"

namespace altruct {
namespace math {

// series storage type

namespace series_storage {
    enum type { INSTANCE, STATIC, CONSTANT };
}

template<typename T, int ID, int STORAGE_TYPE>
struct series_members;

template<typename T, int ID>
struct series_members<T, ID, series_storage::INSTANCE> {
    int _N;
    series_members(int _N = 0) : _N(_N) {}
    int N() const { return _N; }
    int& refN() { return _N; }
};

template<typename T, int ID>
struct series_members<T, ID, series_storage::STATIC> {
    static int _N;
    series_members(int _N = 0) {}
    static int N() { return _N; }
    static int& refN() { return _N; }
};

template<typename T, int ID>
struct series_members<T, ID, series_storage::CONSTANT> {
    series_members(int _N = 0) {}
    static int N() { return ID; }
    static int& refN() { static int _N = 0; return _N; }
};

/**
 * A formal power series
 *
 * `s(x) = p(x) + O(x^N)`, where p(x) is a polynomial of degree N-1
 *
 * @param T - the underlying type
 * @param ID - ID of the series type (useful with series_storage::CONSTANT)
 * @param STORAGE_TYPE - whether N is a constant, static or instance member
 *   See `STORAGE_TYPE` in the `modulo` class for details.
 */
template<typename T, int ID, int STORAGE_TYPE = series_storage::CONSTANT>
class series : public series_members<T, ID, STORAGE_TYPE> {
    typedef series_members<T, ID, STORAGE_TYPE> my_series_members;
public:
    // s(x) = p(x) + O(x^N)
    polynom<T> p;

    series(const T& c0 = T(0)) : my_series_members(1), p(c0) { p.reserve(this->N()); }
    // construct from int, but only if T is not integral to avoid constructor clashing
    template <typename I = T, typename = std::enable_if_t<!std::is_integral<I>::value>>
    series(int c0) : my_series_members(1) , p(c0) { p.reserve(this->N()); } // to allow constructing from 0 and 1
    series(series&& rhs) : my_series_members(rhs.N()), p(std::move(rhs.p)) { p.reserve(this->N()); }
    series(const series& rhs) : my_series_members(rhs.N()), p(rhs.p) { p.reserve(this->N()); }
    series(polynom<T>&& rhs, int _N) : my_series_members(_N), p(std::move(rhs)) { p.reserve(this->N()); }
    series(const polynom<T>& rhs, int _N) : my_series_members(_N), p(rhs) { p.reserve(this->N()); }
    series(polynom<T>&& rhs) : my_series_members(rhs.size()), p(std::move(rhs)) { p.reserve(this->N()); }
    series(const polynom<T>& rhs) : my_series_members(rhs.size()), p(rhs) { p.reserve(this->N()); }
    series(std::vector<T>&& c) : my_series_members((int)c.size()), p(std::move(c)) { p.reserve(this->N()); }
    series(const std::vector<T>& c) : my_series_members((int)c.size()), p(c) { p.reserve(this->N()); }
    template<typename It> series(It begin, It end) : my_series_members((int)std::distance(begin, end)), p(begin, end) { p.reserve(this->N()); }
    series(std::initializer_list<T> list) : my_series_members((int)list.size()), p(list) { p.reserve(this->N()); }
    series& operator = (series&& rhs) { this->refN() = rhs.N(), this->p = std::move(rhs.p); return *this; }
    series& operator = (const series& rhs) { this->refN() = rhs.N(), this->p = rhs.p; return *this; }

    series& swap(series& rhs) { p.swap(rhs.p); std::swap(this->refN(), rhs.refN()); return *this; }

    series& resize(int _N) { this->refN() = _N; p.resize(_N); return *this; }
    int size() const { return this->N(); }
    const T& at(int index) const { return p.at(index); }
    const T& operator [] (int index) const { return p[index]; }
    T& operator [] (int index) { return p[index]; }

    bool operator == (const series& rhs) const { return p == rhs.p; }
    bool operator != (const series& rhs) const { return p != rhs.p; }
    bool operator <  (const series& rhs) const { return p <  rhs.p; }
    bool operator >  (const series& rhs) const { return p >  rhs.p; }
    bool operator <= (const series& rhs) const { return p <= rhs.p; }
    bool operator >= (const series& rhs) const { return p >= rhs.p; }

    series  operator +  (const series &rhs) const { series t(*this); t += rhs; return t; }
    series  operator -  (const series &rhs) const { series t(*this); t -= rhs; return t; }
    series  operator -  ()                  const { return series(-p, this->N()); }
    series  operator *  (const series& rhs) const { series t(*this); t *= rhs; return t; }
    series  operator /  (const series& rhs) const { series t(*this); t /= rhs; return t; }

    series  operator *  (const T &val) const { series t(*this); t *= val; return t; }
    series  operator /  (const T &val) const { series t(*this); t /= val; return t; }

    series& operator += (const series &rhs) { p += rhs.p; p.reserve(this->N()); return *this; }
    series& operator -= (const series &rhs) { p -= rhs.p; p.reserve(this->N()); return *this; }
    series& operator *= (const series& rhs) { polynom<T>::mul(p, p, rhs.p, this->N() - 1); return *this; }
    series& operator /= (const series& rhs) { return *this *= rhs.inverse(); }

    series& operator *= (const T &val) { p *= val; p.reserve(this->N()); return *this; }
    series& operator /= (const T &val) { p /= val; p.reserve(this->N()); return *this; }

    series derivative() const { return series(p.derivative(), this->N()); }
    series integral() const { return integral(p.ZERO_COEFF); }
    series integral(const T& c0) const { auto pi = p.integral(c0); pi.resize(this->N()); return series(std::move(pi), this->N()); }

    // s(x)*x^l - shifts coefficients of s(x) by l places
    series shift(int l) const {
        if (l < 0) {
            polynom<T> t{ p.c.begin() - l, p.c.end() };
            t.c.insert(t.c.end(), -l, p.ZERO_COEFF);
            return series(t, this->N());
        } else {
            polynom<T> t{ p.c.begin(), p.c.end() - l };
            t.c.insert(t.c.begin(), l, p.ZERO_COEFF);
            return series(t, this->N());
        }
    }

    // s(x*a)
    series sub_mul(const T& a) const {
        series s = *this;
        T a_n = id_coeff();
        for (int i = 1; i < this->N(); i++) {
            a_n *= a;
            s[i] *= a_n;
        }
        return s;
    }

    // s(x^a)
    series sub_pow(int a) const {
        series s(polynom<T>{ p.ZERO_COEFF }, this->N());
        if (a > 0) {
            for (int i = 0, j = 0; i < this->N(); i += a, j++) {
                s[i] = p[j];
            }
        }
        return s;
    }

    // t(x) so that s(x) * t(x) == 1 + O(x^N); O(M(N))
    series inverse() const {
        // ensure that p[0] is 1 before inverting
        if (p[0] == zeroT<T>::of(p[0])) return series(polynom<T>{ p.ZERO_COEFF }, this->N());
        if (p[0] != identityT<T>::of(p[0])) return (*this / p[0]).inverse() / p[0];
        polynom<T> r{id_coeff()}, t;
        for (int l = 1; l < this->N() * 2; l *= 2) {
            int m = std::min(this->N() - 1, l), k = l / 2 + 1;
            t.c.assign(p.c.begin(), p.c.begin() + m + 1);
            polynom<T>::mul(t, t, r, l + 1);
            t.c.erase(t.c.begin(), t.c.begin() + k);
            polynom<T>::mul(t, t, r, l - k);
            for (int i = m; i >= k; i--) {
                r[i] = -t[i - k];
            }
        }
        return series(std::move(r), this->N());
    }

    // exp(s(x)) - series expansion of exponential of s(x)
    // the following should hold: s(0) == 0
    series exp() const {
        // See R.P.Brent & H.T.Kung - Fast Algorithms for Manipulating Formal Power Series
        typedef series<T, 0, series_storage::INSTANCE> serx;
        polynom<T> r{ id_coeff() }, t;
        for (int l = 1; l < this->N(); l *= 2) {
            int m = std::min(this->N(), l * 2);
            t.c.assign(p.c.begin(), p.c.begin() + m);
            t -= serx(r, m).ln().p;
            t[0] += id_coeff();
            polynom<T>::mul(r, r, t, m - 1);
        }
        return series(std::move(r), this->N());
    }

    // ln(s(x)) - series expansion of natural logarithm of s(x)
    // the following should hold: s(0) == 1
    series ln() const { return ln(p.ZERO_COEFF); }
    series ln(const T& c0) const {
        return (derivative() / *this).integral(c0);
    }

    // s(x)^a - a-th power of s(x)
    series pow(int64_t a, int64_t threshold = 200) const {
        if (a < threshold) {
            return powT(*this, a);
        } else if (p[0] == p.ZERO_COEFF) {
            int l = p.lowest();
            return shift(-p.lowest()).pow(a).shift(p.lowest() * int(a));
        } else if (p[0] == id_coeff()) {
            return (ln() * castOf(p[0], a)).exp();
        } else {
            return (*this / p[0]).pow(a) * powT(p[0], a);
        }
    }

    // series expansion of exp(a*x) = Sum[a^n * x^n / n!, n]
    static series exp(const T& a, int _N = 0) {
        auto id = identityT<T>::of(a);
        series s(polynom<T>{ id }, _N);
        for (int i = 1; i < s.size(); i++) {
            s[i] = s[i - 1] * a;
        }
        return s.make_exponential();
    }

    // converts the ordinary generating function to the exponential
    // by dividing coefficient of each x^n by n!
    series make_exponential() const {
        series s = *this;
        T fact = id_coeff();
        for (int i = 1; i < this->N(); i++) {
            fact *= i;
        }
        T ifact = id_coeff() / fact;
        for (int i = this->N() - 1; i > 0; i--) {
            s[i] *= ifact;
            ifact *= i;
        }
        return s;
    }

    // converts the exponential generating function to the ordinary
    // by multiplying coefficient of each x^n by n!
    series make_ordinary() const {
        series s = *this;
        T fact = id_coeff();
        for (int i = 1; i < this->N(); i++) {
            fact *= i;
            s[i] *= fact;
        }
        return s;
    }

    // Sum[f(n) * x^n, n]
    template<typename F>
    static series of(F f, int _N = 0) {
        series s(polynom<T>{ f(0) }, _N);
        for (int n = 1; n < s.size(); n++) {
            s[n] = f(n);
        }
        return s;
    }

    // identity coefficient
    T id_coeff() const {
        return identityT<T>::of(p.ZERO_COEFF);
    }
};

template<typename T>
using seriesX = series<T, 0, series_storage::INSTANCE>;

template<typename T, int ID>
int series_members<T, ID, series_storage::STATIC>::_N = ID;

template<typename T, int ID, int STORAGE_TYPE, typename I>
struct castT<series<T, ID, STORAGE_TYPE>, I> {
    typedef series<T, ID, STORAGE_TYPE> ser;
    static ser of(const I& x) {
        return ser(castOf<polynom<T>>(x));
    }
    static ser of(const ser& ref, const I& x) {
        return ser(castOf(ref.p, x), ref.N());
    }
};
template<typename T, int ID, int STORAGE_TYPE>
struct castT<series<T, ID, STORAGE_TYPE>, series<T, ID, STORAGE_TYPE>> : nopCastT<series<T, ID, STORAGE_TYPE>>{};

template<typename T, int ID, int STORAGE_TYPE>
struct identityT<series<T, ID, STORAGE_TYPE>> {
    typedef series<T, ID, STORAGE_TYPE> ser;
    static ser of(const ser& s) {
        return ser(identityOf(s.p), s.N());
    }
};

template<typename T, int ID, int STORAGE_TYPE>
struct zeroT<series<T, ID, STORAGE_TYPE>> {
    typedef series<T, ID, STORAGE_TYPE> ser;
    static ser of(const ser& s) {
        return ser(zeroOf(s.p), s.N());
    }
};

} // math
} // altruct
