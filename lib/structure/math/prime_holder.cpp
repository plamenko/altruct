#include "structure/math/prime_holder.h"

#include <algorithm>

using namespace std;

namespace altruct {
namespace math {

typedef long long ll;
typedef std::pair<int, int> fact_pair;

void prime_holder::ensure_pq() {
	if (vp.empty() || vq.empty()) {
		m = int(1.25 * sz / log(sz)) + 1; // upper bound on pi(sz)
		vp.resize(m);
		vq.resize(sz);
		m = altruct::math::primes(&vp[0], &vq[0], sz);
		vp.resize(m);
	}
}

vector<int>& prime_holder::ensure(vector<int> &v, void(*f)(int*, int, int*, int)) {
	if (v.empty()) {
		v.resize(sz);
		f(&v[0], sz, &p()[0], primes());
	}
	return v;
}

vector<fact_pair> prime_holder::factor_integer(int n) {
	vector<fact_pair> vf;
	altruct::math::factor_integer(vf, n, &pf()[0]);
	sort(vf.begin(), vf.end());
	return vf;
}

vector<fact_pair> prime_holder::factor_integer(vector<int> vn) {
	vector<fact_pair> vf;
	altruct::math::factor_integer(vf, vn, &pf()[0]);
	sort(vf.begin(), vf.end());
	return vf;
}

vector<ll> prime_holder::divisors(const vector<fact_pair> &vf, ll maxd) {
	vector<ll> vd;
	altruct::math::divisors(vd, vf, maxd);
	sort(vd.begin(), vd.end());
	return vd;
}

vector<ll> prime_holder::divisors(int n, ll maxd) {
	return divisors(factor_integer(n), maxd);
}

vector<ll> prime_holder::divisors(const vector<int> &vn, ll maxd) {
	return divisors(factor_integer(vn), maxd);
}

} // math
} // altruct
