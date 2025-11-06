#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#include <vector>
#include <stdexcept>
#include <gsl/gsl_linalg.h>


using Matrix = std::vector<std::vector<double>>;
using Vec = std::vector<double>;


struct SPDMatrix {
// Holds Cholesky factor L of SPD matrix A: A = L L^T
std::size_t n{};
Matrix L; // lower-triangular


static SPDMatrix cholesky(const Matrix& A) {
const std::size_t n = A.size();
if (n == 0) throw std::invalid_argument("Empty matrix");
for (const auto& r : A) if (r.size() != n) throw std::invalid_argument("Non-square matrix");


gsl_matrix* M = gsl_matrix_alloc(n, n);
for (std::size_t i=0;i<n;++i)
for (std::size_t j=0;j<n;++j)
gsl_matrix_set(M, i, j, A[i][j]);


int s = gsl_linalg_cholesky_decomp(M); // in-place, lower stored in M (strict lower + diag)
if (s) { gsl_matrix_free(M); throw std::runtime_error("Cholesky failed: matrix not SPD"); }


SPDMatrix out; out.n = n; out.L.assign(n, Vec(n, 0.0));
for (std::size_t i=0;i<n;++i)
for (std::size_t j=0;j<=i;++j)
out.L[i][j] = gsl_matrix_get(M, i, j);


gsl_matrix_free(M);
return out;
}


// Solve A x = b via L L^T x = b
Vec solve(const Vec& b) const {
if (b.size() != n) throw std::invalid_argument("Dimension mismatch");
Vec y(n,0.0), x(n,0.0);
// forward: L y = b
for (std::size_t i=0;i<n;++i) {
double acc = b[i];
for (std::size_t k=0;k<i;++k) acc -= L[i][k]*y[k];
y[i] = acc / L[i][i];
}
// backward: L^T x = y
for (std::size_t i=n; i-- > 0;) {
double acc = y[i];
for (std::size_t k=i+1;k<n;++k) acc -= L[k][i]*x[k];
x[i] = acc / L[i][i];
}
return x;
}


// Quadratic form b^T A^{-1} b
double quad_inv(const Vec& b) const {
Vec x = solve(b);
double v = 0.0; for (std::size_t i=0;i<n;++i) v += b[i]*x[i];
return v;
}
};

#endif