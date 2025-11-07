#include "RandomVectorGenerator.h"



void CorrelationMatrixValidator::validate(const Matrix& R) const {
    const size_t n = R.size();
    if (n == 0) throw std::invalid_argument("Invalid matrix.");
    for (const auto& row : R) {
        if (row.size() != n) throw std::invalid_argument("Non squared matrix.");
    }

    //->check sym
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            if (std::fabs(R[i][j] - R[j][i]) > EPS_SYM) {
                throw std::invalid_argument("Non symmetric matrix.");
            }
        }
    }
    //->check diag
    // for (size_t i = 0; i < n; ++i) {
    //     if (std::fabs(R[i][i] - 1.0) > EPS_DIAG) {
    //         throw std::invalid_argument("Diagonal elements of the matrix needs to be ones");
    //     }
    // }
    // positive-dev check with Cholesky
}


Matrix CholeskyDecomposition::factorize(const Matrix& R) {
    const size_t n = R.size();
    Matrix L(n, std::vector<double>(n, 0.0));

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            double sum = R[i][j];
            for (size_t k = 0; k < j; ++k) {
                sum -= L[i][k] * L[j][k];
            }
            if (i == j) {
                if (sum <= 0.0) {
                    throw std::invalid_argument(
                        "Matrix is not positive-définie (for Cholesky).");
                }
                L[i][j] = std::sqrt(sum);
            } else {
                L[i][j] = sum / L[j][j];
            }
        }
    }
    return L;
}

RandomVectorGenerator::RandomVectorGenerator(std::unique_ptr<IDistribution> dist,
                        std::unique_ptr<IDecomposition> decomp)
    : dist_(std::move(dist)), decomp_(std::move(decomp)) {}

// y = L * z, z ~ i.i.d. (E=0, Var=1). Cov(y) = L L^T = R.
Vector RandomVectorGenerator::generate(const Matrix& correlation) const {
    CorrelationMatrixValidator().validate(correlation);
    Matrix L = decomp_->factorize(correlation);
    const std::size_t n = L.size();

    Vector z = dist_->sample(n);
    Vector y(n, 0.0);

    for (std::size_t i = 0; i < n; ++i) {
        double acc = 0.0;
        for (std::size_t k = 0; k <= i; ++k) {
            acc += L[i][k] * z[k];
        }
        y[i] = acc;
    }
    return y;
}


