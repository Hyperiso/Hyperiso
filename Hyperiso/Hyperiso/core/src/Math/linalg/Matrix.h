#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <vector>
#include <array>
#include <memory>
#include <stdexcept>
#include <algorithm>
#include "Math.h"

// Minimal RAII for GSL pointers

using gsl_matrix_ptr = std::unique_ptr<gsl_matrix, decltype(&gsl_matrix_free)>;
using gsl_vector_ptr = std::unique_ptr<gsl_vector, decltype(&gsl_vector_free)>;
using gsl_permutation_ptr = std::unique_ptr<gsl_permutation, decltype(&gsl_permutation_free)>;
using gsl_eigen_workspace_ptr = std::unique_ptr<gsl_eigen_symmv_workspace, decltype(&gsl_eigen_symmv_free)>;

inline gsl_matrix_ptr make_gsl_matrix(size_t rows, size_t cols) {
    return gsl_matrix_ptr(gsl_matrix_alloc(rows, cols), &gsl_matrix_free);
}

inline gsl_vector_ptr make_gsl_vector(size_t size) {
    return gsl_vector_ptr(gsl_vector_alloc(size), &gsl_vector_free);
}

inline gsl_permutation_ptr make_gsl_permutation(size_t size) {
    return gsl_permutation_ptr(gsl_permutation_alloc(size), &gsl_permutation_free);
}

inline gsl_eigen_workspace_ptr make_eigen_workspace(size_t n) {
    return { gsl_eigen_symmv_alloc(n), &gsl_eigen_symmv_free };
}


// -------------------------------------------------------

struct EigenSystem {
    RealMatrix D;
    RealMatrix P;
};

struct SignedLogDet {
    double logdet;
    int sign;
};

class RealMatrix {
public:
    RealMatrix() = default;
    RealMatrix(std::vector<double> data, std::size_t rows, std::size_t cols);
    RealMatrix(std::vector<std::vector<double>> data);
    RealMatrix(std::size_t rows, std::size_t cols);

    // io

    double& at(size_t i, size_t j);
    double& unchecked_at(size_t i, size_t j);
    const double& at(size_t i, size_t j) const;
    const double& unchecked_at(size_t i, size_t j) const;

    std::size_t rows() const;
    std::size_t cols() const;

    // GSL support

    static RealMatrix from_gsl_copy(const gsl_matrix* A);
    gsl_matrix_ptr to_gsl_matrix() const;

    // Linear algebra

    bool is_symmetric() const;
    RealMatrix transpose() const;
    EigenSystem eig() const;
    SignedLogDet slogdet() const;
    RealMatrix inv() const;

    // Arithmetics

    RealMatrix operator-() const;
    RealMatrix& operator+=(const RealMatrix& rhs);
    RealMatrix& operator-=(const RealMatrix& rhs);
    RealMatrix& operator*=(const RealMatrix& rhs);
    RealMatrix& operator*=(double scalar);
    RealMatrix& operator/=(double scalar);

    friend RealMatrix operator+(RealMatrix lhs, const RealMatrix& rhs) {
        lhs += rhs;
        return lhs;
    }

    friend RealMatrix operator-(RealMatrix lhs, const RealMatrix& rhs) {
        lhs -= rhs;
        return lhs;
    }

    friend RealMatrix operator*(RealMatrix lhs, const RealMatrix& rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend RealMatrix operator*(RealMatrix lhs, double scalar) {
        lhs *= scalar;
        return lhs;
    }

    friend RealMatrix operator*(double scalar, RealMatrix rhs) {
        rhs *= scalar;
        return rhs;
    }

    friend RealMatrix operator/(RealMatrix lhs, double scalar) {
        lhs /= scalar;
        return lhs;
    }
    
private:
    std::vector<double> data;
    std::size_t rows_, cols_;
};

RealMatrix eye(std::size_t n);
RealMatrix diag(gsl_vector* X);
RealMatrix nearest_psd(RealMatrix R, double thr = 1e-12);
RealMatrix cholesky_L(RealMatrix R);

#endif // __MATRIX_H__