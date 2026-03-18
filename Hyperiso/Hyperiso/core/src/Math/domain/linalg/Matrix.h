#ifndef MATRIX_H
#define MATRIX_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <vector>
#include <array>
#include <memory>
#include <stdexcept>
#include <algorithm>

#include "../special/special_generic.h"

/**
 * @file Matrix.h
 * @brief Lightweight dense real-matrix utilities built on top of STL storage and GSL backends.
 *
 * This file provides:
 * - a small RAII layer for owning GSL pointers safely,
 * - the @ref RealMatrix class for dense real-valued matrices,
 * - the @ref EigenSystem and @ref SignedLogDet helper structs,
 * - several utility constructors and decompositions:
 *   - identity matrix,
 *   - diagonal matrix from a GSL vector,
 *   - nearest positive semi-definite projection,
 *   - Cholesky factor extraction,
 *   - block-diagonal concatenation.
 *
 * Design goals:
 * - keep ownership simple using `std::vector<double>` row-major storage,
 * - expose a minimal, readable matrix API,
 * - delegate numerically sensitive operations (LU, eigensystem, Cholesky)
 *   to GSL,
 * - provide enough functionality for the statistics / covariance layer.
 *
 * Storage convention:
 * - matrices are stored in row-major order,
 * - element `(i,j)` is stored at `data[i * cols + j]`.
 *
 * @note This is a real-matrix utility only. Complex-valued matrix algebra is
 *       intentionally out of scope here.
 */

// -----------------------------------------------------------------------------
// RAII wrappers around GSL resources
// -----------------------------------------------------------------------------

/// Convenient alias for a dense vector of real values.
using Vector = std::vector<double>;

/// Owning smart pointer for a GSL matrix.
using gsl_matrix_sptr = std::unique_ptr<gsl_matrix, decltype(&gsl_matrix_free)>;

/// Owning smart pointer for a GSL vector.
using gsl_vector_sptr = std::unique_ptr<gsl_vector, decltype(&gsl_vector_free)>;

/// Owning smart pointer for a GSL permutation.
using gsl_permutation_sptr = std::unique_ptr<gsl_permutation, decltype(&gsl_permutation_free)>;

/// Owning smart pointer for a GSL symmetric-eigensystem workspace.
using gsl_eigen_workspace_sptr = std::unique_ptr<gsl_eigen_symmv_workspace, decltype(&gsl_eigen_symmv_free)>;

/**
 * @brief Allocates a GSL matrix with automatic RAII destruction.
 * @param rows Number of rows.
 * @param cols Number of columns.
 * @return Unique pointer owning the allocated matrix.
 */
inline gsl_matrix_sptr make_gsl_matrix(size_t rows, size_t cols) {
    return gsl_matrix_sptr(gsl_matrix_alloc(rows, cols), &gsl_matrix_free);
}

/**
 * @brief Allocates a GSL vector with automatic RAII destruction.
 * @param size Vector size.
 * @return Unique pointer owning the allocated vector.
 */
inline gsl_vector_sptr make_gsl_vector(size_t size) {
    return gsl_vector_sptr(gsl_vector_alloc(size), &gsl_vector_free);
}

/**
 * @brief Allocates a GSL permutation with automatic RAII destruction.
 * @param size Permutation size.
 * @return Unique pointer owning the allocated permutation.
 */
inline gsl_permutation_sptr make_gsl_permutation(size_t size) {
    return gsl_permutation_sptr(gsl_permutation_alloc(size), &gsl_permutation_free);
}

/**
 * @brief Allocates a GSL symmetric eigensolver workspace with RAII destruction.
 * @param n Matrix dimension.
 * @return Unique pointer owning the workspace.
 */
inline gsl_eigen_workspace_sptr make_eigen_workspace(size_t n) {
    return { gsl_eigen_symmv_alloc(n), &gsl_eigen_symmv_free };
}


// -------------------------------------------------------

struct EigenSystem; 

/**
 * @struct SignedLogDet
 * @brief Signed logarithmic determinant representation.
 *
 * This is the standard robust decomposition of a determinant:
 * - `logdet = log(|det(A)|)`
 * - `sign   = sign(det(A))`
 *
 * It is useful when determinants may span many orders of magnitude and direct
 * evaluation would underflow or overflow.
 */
struct SignedLogDet {
    double logdet;  /// Natural logarithm of the absolute determinant.
    int sign;       /// Sign of the determinant (`-1`, `0`, or `+1` depending on decomposition result).
};

/**
 * @class RealMatrix
 * @brief Dense real-valued matrix with STL storage and GSL-backed linear algebra.
 *
 * @ref RealMatrix is a small utility class representing a dense matrix of
 * double-precision real numbers. It stores data in a single contiguous
 * `std::vector<double>` in row-major order.
 *
 * Features include:
 * - checked and unchecked element access,
 * - shape queries,
 * - row / column removal,
 * - conversion to/from GSL matrices,
 * - symmetry checks,
 * - transpose,
 * - eigendecomposition for symmetric matrices,
 * - LU-based inverse and signed log-determinant,
 * - basic arithmetic operations.
 *
 * The class is intentionally lightweight and does not attempt to be a full
 * linear algebra framework.
 */
class RealMatrix {
public:
    /**
     * @brief Constructs an empty 0×0 matrix.
     */
    RealMatrix() : rows_(0), cols_(0) {}

    /**
     * @brief Constructs a matrix from flat row-major data.
     *
     * @param data_ Flat row-major data buffer.
     * @param rows  Number of rows.
     * @param cols  Number of columns.
     *
     * @throws std::invalid_argument if `data_.size() != rows * cols`.
     */
    RealMatrix(std::vector<double> data_, std::size_t rows, std::size_t cols);

    /**
     * @brief Constructs a matrix from nested row vectors.
     *
     * @param data_ Matrix data as `data_[row][col]`.
     *
     * @throws std::invalid_argument if the outer vector is empty or if rows do
     *         not all have the same length.
     */
    RealMatrix(std::vector<std::vector<double>> data_);

    /**
     * @brief Constructs a zero-initialized matrix of given shape.
     *
     * @param rows Number of rows.
     * @param cols Number of columns.
     */
    RealMatrix(std::size_t rows, std::size_t cols);

     // -------------------------------------------------------------------------
    // Element access / shape
    // -------------------------------------------------------------------------

    /**
     * @brief Returns a mutable reference to element `(i,j)` with bounds checking.
     *
     * @param i Row index.
     * @param j Column index.
     * @return Mutable reference to the requested entry.
     *
     * @throws std::out_of_range if `(i,j)` lies outside the matrix shape.
     */
    double& at(size_t i, size_t j);

    /**
     * @brief Returns a mutable reference to element `(i,j)` without bounds checking.
     *
     * @param i Row index.
     * @param j Column index.
     * @return Mutable reference to the requested entry.
     *
     * @warning No bounds check is performed.
     */
    double& unchecked_at(size_t i, size_t j);

    /**
     * @brief Returns a const reference to element `(i,j)` with bounds checking.
     *
     * @param i Row index.
     * @param j Column index.
     * @return Const reference to the requested entry.
     *
     * @throws std::out_of_range if `(i,j)` lies outside the matrix shape.
     */
    const double& at(size_t i, size_t j) const;

    /**
     * @brief Returns a const reference to element `(i,j)` without bounds checking.
     *
     * @param i Row index.
     * @param j Column index.
     * @return Const reference to the requested entry.
     *
     * @warning No bounds check is performed.
     */
    const double& unchecked_at(size_t i, size_t j) const;

    /**
     * @brief Returns the number of rows.
     */
    std::size_t rows() const;

    /**
     * @brief Returns the number of columns.
     */
    std::size_t cols() const;

    // -------------------------------------------------------------------------
    // Row / column manipulation
    // -------------------------------------------------------------------------

    /**
     * @brief Removes one row from the matrix.
     *
     * All rows below @p row_idx are shifted up by one.
     *
     * @param row_idx Index of the row to remove.
     *
     * @warning The implementation assumes `row_idx < rows()`.
     *          It does not currently throw on invalid indices.
     */
    void remove_row(std::size_t row_idx);

    /**
     * @brief Removes one column from the matrix.
     *
     * All columns to the right of @p col_idx are shifted left by one.
     *
     * @param col_idx Index of the column to remove.
     *
     * @warning The implementation assumes `col_idx < cols()`.
     *          It does not currently throw on invalid indices.
     */
    void remove_column(std::size_t col_idx);

    /**
     * @brief Removes one row and one column with the same index.
     *
     * This is useful for principal-submatrix extraction, e.g. when removing one
     * nuisance direction from a covariance / correlation matrix.
     *
     * @param dim_idx Index of the row and column to remove.
     *
     * @warning The implementation assumes a valid index.
     * @note Implementation detail: the rewritten buffer must use the new column
     *       count (`cols()-1`) and the shifted column index.
     */
    void remove_row_and_column(std::size_t dim_idx);

    // -------------------------------------------------------------------------
    // GSL support
    // -------------------------------------------------------------------------

    /**
     * @brief Creates a RealMatrix by copying data from a GSL matrix.
     *
     * @param A Source GSL matrix.
     * @return Copied matrix in row-major STL storage.
     */
    static RealMatrix from_gsl_copy(const gsl_matrix* A);

    /**
     * @brief Converts this matrix to a newly allocated GSL matrix.
     *
     * @return RAII-owned GSL matrix containing a full copy of the data.
     */
    gsl_matrix_sptr to_gsl_matrix() const;

    // -------------------------------------------------------------------------
    // Linear algebra
    // -------------------------------------------------------------------------

    /**
     * @brief Checks whether the matrix is symmetric.
     *
     * The check is exact up to the project helper `fpeq(...)`.
     *
     * @return True if the matrix is square and symmetric, false otherwise.
     */
    bool is_symmetric() const;

    /**
     * @brief Returns the transpose of the matrix.
     *
     * @return Transposed matrix.
     */
    RealMatrix transpose() const;

    /**
     * @brief Computes the eigensystem of a symmetric matrix.
     *
     * Uses GSL symmetric eigendecomposition and sorts eigenvalues in descending order.
     *
     * @return A struct containing:
     * - `D`: diagonal matrix of eigenvalues,
     * - `P`: matrix of eigenvectors.
     *
     * @throws std::invalid_argument if the matrix is not square.
     * @throws std::runtime_error if the matrix is not symmetric.
     */
    EigenSystem eig() const;

    /**
     * @brief Computes the signed logarithmic determinant via LU decomposition.
     *
     * @return Signed-logdet representation of the determinant.
     *
     * @throws std::invalid_argument if the matrix is not square.
     */
    SignedLogDet slogdet() const;

    /**
     * @brief Computes the inverse of the matrix via LU decomposition.
     *
     * @return Inverse matrix.
     *
     * @throws std::invalid_argument if the matrix is not square.
     * @throws std::runtime_error if the matrix is singular.
     */
    RealMatrix inv() const;

    // -------------------------------------------------------------------------
    // Arithmetic
    // -------------------------------------------------------------------------

    /**
     * @brief Unary minus.
     * @return Matrix with all coefficients negated.
     */
    RealMatrix operator-() const;

    /**
     * @brief In-place matrix addition.
     *
     * @param rhs Right-hand side matrix.
     * @return Reference to `*this`.
     *
     * @throws std::invalid_argument if shapes differ.
     */
    RealMatrix& operator+=(const RealMatrix& rhs);

    /**
     * @brief In-place matrix subtraction.
     *
     * @param rhs Right-hand side matrix.
     * @return Reference to `*this`.
     *
     * @throws std::invalid_argument if shapes differ.
     */
    RealMatrix& operator-=(const RealMatrix& rhs);

    /**
     * @brief In-place matrix multiplication.
     *
     * Performs dense matrix multiplication using a blocked algorithm.
     *
     * @param rhs Right-hand side matrix.
     * @return Reference to `*this`.
     *
     * @throws std::invalid_argument if inner dimensions do not match.
     */
    RealMatrix& operator*=(const RealMatrix& rhs);

    /**
     * @brief In-place scalar multiplication.
     *
     * @param scalar Multiplicative factor.
     * @return Reference to `*this`.
     */
    RealMatrix& operator*=(double scalar);

    /**
     * @brief In-place scalar division.
     *
     * @param scalar Divisor.
     * @return Reference to `*this`.
     *
     * @throws std::invalid_argument if scalar is numerically zero.
     */
    RealMatrix& operator/=(double scalar);

    /**
     * @brief Matrix addition.
     */
    friend RealMatrix operator+(RealMatrix lhs, const RealMatrix& rhs) {
        lhs += rhs;
        return lhs;
    }

    /**
     * @brief Matrix subtraction.
     */
    friend RealMatrix operator-(RealMatrix lhs, const RealMatrix& rhs) {
        lhs -= rhs;
        return lhs;
    }

    /**
     * @brief Matrix product.
     */
    friend RealMatrix operator*(RealMatrix lhs, const RealMatrix& rhs) {
        lhs *= rhs;
        return lhs;
    }

    /**
     * @brief Right scalar multiplication.
     */
    friend RealMatrix operator*(RealMatrix lhs, double scalar) {
        lhs *= scalar;
        return lhs;
    }

    /**
     * @brief Left scalar multiplication.
     */
    friend RealMatrix operator*(double scalar, RealMatrix rhs) {
        rhs *= scalar;
        return rhs;
    }

    /**
     * @brief Scalar division.
     */
    friend RealMatrix operator/(RealMatrix lhs, double scalar) {
        lhs /= scalar;
        return lhs;
    }

    /**
     * @brief Stream output helper.
     *
     * Prints the matrix in a nested bracket form.
     *
     * @param os Output stream.
     * @param A  Matrix to print.
     * @return Output stream.
     */
    friend std::ostream& operator<<(std::ostream& os, RealMatrix A);
    
private:
    std::vector<double> data;   /// Flat row-major data buffer.
    std::size_t rows_, cols_;   /// Number of rows and columns.
};

/**
 * @struct EigenSystem
 * @brief Container for an eigendecomposition.
 *
 * Convention used here:
 * - `D` is the diagonal eigenvalue matrix,
 * - `P` is the eigenvector matrix.
 *
 * For a symmetric matrix `A`, the decomposition is intended to satisfy:
 * @code
 *   A = P * D * P.transpose()
 * @endcode
 */
struct EigenSystem {
    RealMatrix D;   /// Diagonal matrix of eigenvalues.
    RealMatrix P;   /// Matrix of eigenvectors.
};

/**
 * @brief Returns the identity matrix of size `n`.
 *
 * @param n Matrix dimension.
 * @return `n × n` identity matrix.
 */
RealMatrix eye(std::size_t n);

/**
 * @brief Builds a diagonal matrix from a GSL vector.
 *
 * @param X Input vector.
 * @return Square diagonal matrix whose diagonal entries are the values of @p X.
 */
RealMatrix diag(const gsl_vector* X);

/**
 * @brief Projects a matrix to the nearest positive semi-definite correlation-like matrix.
 *
 * Algorithm:
 * - symmetrize the input,
 * - diagonalize it,
 * - clip eigenvalues below @p thr,
 * - reconstruct the matrix,
 * - renormalize the diagonal to 1,
 * - symmetrize again.
 *
 * This is especially useful for repairing numerically inconsistent correlation matrices.
 *
 * @param R   Input matrix.
 * @param thr Lower bound applied to eigenvalues during clipping.
 * @return Positive semi-definite, unit-diagonal repaired matrix.
 *
 * @throws std::invalid_argument if the matrix is not square.
 */
RealMatrix nearest_psd(RealMatrix R, double thr = 1e-12);

/**
 * @brief Returns the lower-triangular Cholesky factor of a PSD matrix.
 *
 * If `R = L L^T`, this function returns `L`.
 *
 * @param R Input matrix.
 * @return Lower-triangular Cholesky factor.
 *
 * @throws std::invalid_argument if the matrix is not square.
 * @throws std::runtime_error if Cholesky decomposition fails.
 */
RealMatrix cholesky_L(RealMatrix R);

/**
 * @brief Builds a block-diagonal matrix from several square blocks.
 *
 * @param blocks List of square matrices.
 * @return Block-diagonal concatenation of all input matrices.
 *
 * @throws std::invalid_argument if one input block is not square.
 */
RealMatrix block_diag(const std::vector<RealMatrix>& blocks);

/**
 * @brief Variadic convenience overload for @ref block_diag.
 *
 * @tparam Ms Matrix-like arguments convertible to @ref RealMatrix.
 * @param ms  Input blocks.
 * @return Block-diagonal concatenation.
 */
template<typename... Ms>
RealMatrix block_diag(const Ms&... ms) {
    return block_diag(std::vector<RealMatrix>{ms...});
}

#endif // MATRIX_H