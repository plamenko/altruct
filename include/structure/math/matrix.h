#pragma once

#include "algorithm/math/base.h"

#include <vector>

namespace altruct {
namespace math {

template<typename T>
class matrix {
public:
	typedef std::vector<T> row_type;

	std::vector<row_type> a;

	matrix() : matrix(0, 0) {
	}

	matrix(int n, int m) {
		if (!m) m = n;
		a.resize(n);
		for (auto &row : a) {
			row.resize(m);
		}
	}

	matrix(std::initializer_list<row_type> list) : a(list) {}

	matrix(const matrix& rhs, int i, int j, int n = 0, int m = 0) {
		if (!n) n = rhs.rows() - i;
		if (!m) m = rhs.cols() - j;
		a.resize(n);
		for (int k = 0; k < n; k++) {
			auto it = rhs[k + i].begin() + j;
			a[k].assign(it, it + m);
		}
	}

	int rows() const {
		return (int)a.size();
	}

	int cols() const {
		return a.empty() ? 0 : (int)a.back().size();
	}

	matrix& swap(matrix &rhs) {
		std::swap(a, rhs.a);
		return (*this);
	}

	row_type& operator [] (int i) { return a[i]; }
	const row_type& operator [] (int i) const { return a[i]; }

	bool operator == (const matrix &rhs) const { return a == rhs.a; }
	bool operator != (const matrix &rhs) const { return a != rhs.a; }
	bool operator <  (const matrix &rhs) const { return a <  rhs.a; }
	bool operator >  (const matrix &rhs) const { return a >  rhs.a; }
	bool operator <= (const matrix &rhs) const { return a <= rhs.a; }
	bool operator >= (const matrix &rhs) const { return a >= rhs.a; }

	// matrices must be of same dimensions
	matrix& operator += (const matrix &rhs) {
		int n = rows(), m = cols();
		for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			a[i][j] += rhs[i][j];
		return *this;
	}
	matrix operator + (const matrix &rhs) const {
		return matrix(*this) += rhs;
	}

	// matrices must be of same dimensions
	matrix& operator -= (const matrix &rhs) {
		int n = rows(), m = cols();
		for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			a[i][j] -= rhs[i][j];
		return *this;
	}
	matrix operator - (const matrix &rhs) const {
		return matrix(*this) -= rhs;
	}

	// lhs.cols() must be equal to rhs.rows()
	matrix& operator *= (const matrix &rhs) {
		int n = rows(), m = cols(), p = rhs.cols();
		matrix t(n, p);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < p; j++) {
				t[i][j] = 0;
				for (int k = 0; k < m; k++)
					t[i][j] += a[i][k] * rhs[k][j];
			}
		}
		return swap(t);
	}
	matrix operator * (const matrix &rhs) const {
		return matrix(*this) *= rhs;
	}

	matrix& operator *= (const T& s) {
		int n = rows(), m = cols();
		for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			a[i][j] *= s;
		return *this;
	}
	matrix operator * (const T& s) const {
		return matrix(*this) *= s;
	}

	matrix& operator /= (const matrix &rhs) {
		return (*this) *= rhs.inverse();
	}
	matrix operator / (const matrix &rhs) const {
		return matrix(*this) *= rhs.inverse();
	}

	matrix& operator /= (const T& s) {
		int n = rows(), m = cols();
		for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			a[i][j] /= s;
		return *this;
	}
	matrix operator / (const T& s) const {
		return matrix(*this) /= s;
	}

	// matrix must be a square matrix
	matrix pow(long long p) const {
		if (p < 0) {
			return powT(inverse(), -p);
		} else {
			return powT(*this, p);
		}
	}

	// matrix must be a square matrix
	static bool gauss(matrix &mat, matrix &inv, T &det) {
		T e0 = zeroT<T>::of(mat.a[0][0]), e1 = identityT<T>::of(mat.a[0][0]);
		int i, j, k, n = mat.rows();
		inv = identity(n);
		det = e1;
		for (j = 0; j < n; j++) {
			for (i = j; i < n && mat.a[i][j] == e0; i++);
			// matrix singular?
			if (i >= n) {
				inv = matrix(n, n);
				det = e0;
				return false;
			}
			// rows transposition
			if (i != j) {
				std::swap(mat.a[i], mat.a[j]);
				std::swap(inv.a[i], inv.a[j]);
				det = -det;
			}
			// normalize row
			T pi = e1 / mat.a[j][j];
			det *= mat.a[j][j];
			for (k = 0; k < n; k++) {
				mat.a[j][k] *= pi;
				inv.a[j][k] *= pi;
			}
			// eliminate below
			for (i = j + 1; i < n; i++) {
				T p = mat.a[i][j]; if (p == e0) continue;
				for (k = 0; k < n; k++) {
					mat.a[i][k] -= mat.a[j][k] * p;
					inv.a[i][k] -= inv.a[j][k] * p;
				}
			}
		}
		for (j = n - 1; j >= 0; j--) {
			// eliminate above
			for (i = j - 1; i >= 0; i--) {
				T p = mat.a[i][j]; if (p == e0) continue;
				for (k = 0; k < n; k++) {
					mat.a[i][k] -= mat.a[j][k] * p;
					inv.a[i][k] -= inv.a[j][k] * p;
				}
			}
		}
		return true;
	}

	// matrix must be a square matrix
	matrix inverse() const {
		matrix mat(*this), inv; T det;
		gauss(mat, inv, det);
		return inv;
	}

	// matrix must be a square matrix
	T det() const {
		matrix mat(*this), inv; T det;
		gauss(mat, inv, det);
		return det;
	}

	matrix transpose() const {
		int n = rows(), m = cols();
		matrix t(m, n);
		for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			t[j][i] = a[i][j];
		return t;
	}

	static matrix identity(int n) {
		matrix t(n, n);
		for (int i = 0; i < n; i++) {
			t.a[i][i] = 1;
		}
		return t;
	}
};

template<typename T>
struct identityT<matrix<T>> {
	static matrix<T> of(const matrix<T>& x) {
		return matrix<T>::identity(x.rows());
	}
};

template<typename T>
struct zeroT<matrix<T>> {
	static matrix<T> of(const matrix<T>& x) {
		return matrix<T>(x.rows(), x.cols());
	}
};

} // math
} // altruct
