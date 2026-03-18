#include "Matrix.h"

std::ostream &operator<<(std::ostream &os, RealMatrix A) {
    os << "[";
    for (size_t i = 0; i < A.rows(); i++) {
        os << (i == 0 ? "[ " : " [ ");
        for (size_t j = 0; j < A.cols(); j++) {
            os << A.at(i, j) << (j == A.cols() - 1 ? " " : ", ");
        }
        if (i != A.rows() - 1)
            os << "]" << std::endl;
    }
    os << "]]" << std::endl;
    return os;
}

RealMatrix eye(std::size_t n)
{
    std::vector<double> new_data(n * n, 0.0);

    for (size_t i = 0; i < n; ++i) {
        new_data[i * (n + 1)] = 1.0;
    }

    return RealMatrix(new_data, n, n);
}

RealMatrix diag(const gsl_vector *X) {
    const std::size_t n = X->size;
    std::vector<double> new_data(n * n, 0.0);

    for (size_t i = 0; i < n; ++i) {
        new_data[i * (n + 1)] = gsl_vector_get(X, i);
    }

    return RealMatrix(new_data, n, n);
}

RealMatrix nearest_psd(RealMatrix R, double thr) {
    if (R.rows() != R.cols())
        throw std::invalid_argument("Matrix should be square");

    size_t n = R.rows();

    // Symmetrize
    R = 0.5 * (R + R.transpose());

    // Compute eigensystem, clip eigenvalues and recombine
    EigenSystem e = R.eig();
    for (size_t i = 0; i < n; i++) {
        e.D.at(i, i) = std::max(e.D.at(i, i), thr);
    }
    RealMatrix R_psd = e.P * e.D * e.P.transpose();

    // Normalize, re-symmetrize and enforce unit correlations on the diagonal
    std::vector<double> inv_sqrt_diag(n);
    for (size_t i = 0; i < n; ++i)
        inv_sqrt_diag[i] = 1.0 / std::sqrt(R_psd.at(i, i));

    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            R_psd.at(i, j) *= inv_sqrt_diag[i] * inv_sqrt_diag[j];

    R_psd = 0.5 * (R_psd + R_psd.transpose());

    for (size_t i = 0; i < n; i++)
        R_psd.at(i, i) = 1.0;

    return R_psd;
}

RealMatrix cholesky_L(RealMatrix R) {
    if (R.rows() != R.cols())
        throw std::invalid_argument("Matrix should be square");

    gsl_matrix_sptr R_gsl = R.to_gsl_matrix();
    
    if (gsl_linalg_cholesky_decomp1(R_gsl.get()) != GSL_SUCCESS)
        throw std::runtime_error("Cholesky decomposition failed (matrix not PSD)");

    const size_t n = R.rows();
    RealMatrix L(n, n);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            L.at(i, j) = gsl_matrix_get(R_gsl.get(), i, j);
        }
    }

    return L;
}

RealMatrix block_diag(const std::vector<RealMatrix> &blocks) {
    if (blocks.empty())
        return RealMatrix{};

    std::size_t total_dim = 0;

    for (const auto& B : blocks) {
        if (B.rows() != B.cols())
            throw std::invalid_argument("block_diag: matrices must be square");
        total_dim += B.rows();
    }

    RealMatrix M(total_dim, total_dim);

    std::size_t offset = 0;
    for (const auto& B : blocks) {
        const std::size_t n = B.rows();
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                M.unchecked_at(offset + i, offset + j) = B.unchecked_at(i, j);
        offset += n;
    }

    return M;
}

RealMatrix::RealMatrix(std::vector<double> data_, std::size_t rows, std::size_t cols) : data(data_), rows_(rows), cols_(cols) {
    if (data.size() != rows * cols)
        throw std::invalid_argument("Data size does not match matrix shape.");
}

RealMatrix::RealMatrix(std::vector<std::vector<double>> data) {
    if (data.size() == 0)
        throw std::invalid_argument("Matrix should not be sizeless.");

    this->rows_ = data.size();

    std::size_t size_2 = data[0].size();
    for (auto& row : data) {
        if (row.size() != size_2)
            throw std::invalid_argument("All rows should have the same number of columns.");
    }

    this->cols_ = size_2;

    std::vector<double> new_data (rows_ * cols_, 0.0);
    for (size_t i = 0; i < rows_; i++) {
        for (size_t j = 0; j < cols_; j++) {
            new_data[i * cols_ + j] = data[i][j];
        }
    }

    this->data = std::move(new_data);
}

RealMatrix::RealMatrix(std::size_t rows, std::size_t cols) : data(rows * cols, 0.0), rows_(rows), cols_(cols) 
{}

double &RealMatrix::at(size_t i, size_t j) {
    if (i >= rows_ || j >= cols_) 
        throw std::out_of_range("Indices outside matrix shape");

    return this->data[i * cols_ + j];
}

const double &RealMatrix::at(size_t i, size_t j) const {
    if (i >= rows_ || j >= cols_) 
        throw std::out_of_range("Indices outside matrix shape");

    return this->data[i * cols_ + j];
}

std::size_t RealMatrix::rows() const {
    return this->rows_;
}

std::size_t RealMatrix::cols() const {
    return this->cols_;
}

void RealMatrix::remove_row(std::size_t row_idx) {
    std::vector<double> new_data = std::vector<double>((this->rows_ - 1) * this->cols_, 0.0);

    for (size_t i = 0; i < this->rows_; i++) {
        if (i == row_idx) continue;
        std::size_t new_i = i < row_idx ? i : i - 1;
        for (size_t j = 0; j < this->cols_; j++) {
            new_data[new_i * this->cols_ + j] = this->at(i, j);
        }
    }

    this->data = std::move(new_data);
    this->rows_--;
}

void RealMatrix::remove_column(std::size_t col_idx) {
    std::vector<double> new_data = std::vector<double>(this->rows_ * (this->cols_ - 1), 0.0);

    for (size_t j = 0; j < this->cols_; j++) {
        if (j == col_idx) continue;
        std::size_t new_j = j < col_idx ? j : j - 1;
        for (size_t i = 0; i < this->rows_; i++) {
            new_data[i * (this->cols_ - 1) + new_j] = this->at(i, j);
        }
    }

    this->data = std::move(new_data);
    this->cols_--;
}

void RealMatrix::remove_row_and_column(std::size_t dim_idx) {
    std::vector<double> new_data((this->rows_ - 1) * (this->cols_ - 1), 0.0);

    for (size_t i = 0; i < this->rows_; i++) {
        if (i == dim_idx) continue;
        std::size_t new_i = i < dim_idx ? i : i - 1;

        for (size_t j = 0; j < this->cols_; j++) {
            if (j == dim_idx) continue;
            std::size_t new_j = j < dim_idx ? j : j - 1;
            new_data[new_i * (this->cols_ - 1) + new_j] = this->at(i, j);
        }
    }

    this->data = std::move(new_data);
    this->rows_--;
    this->cols_--;
}

RealMatrix RealMatrix::from_gsl_copy(const gsl_matrix* A) {
    const size_t n = A->size1;
    const size_t m = A->size2;
    std::vector<double> data(n * m, 0.0);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            data[i * m + j] = gsl_matrix_get(A, i, j);
        }
    }

    return RealMatrix(std::move(data), n, m);
}

gsl_matrix_sptr RealMatrix::to_gsl_matrix() const {
    gsl_matrix_sptr A = make_gsl_matrix(rows_, cols_);

    for (size_t i = 0; i < rows_; ++i) {
        for (size_t j = 0; j < cols_; ++j) {
            gsl_matrix_set(A.get(), i, j, this->data[i * cols_ + j]);
        }
    }

    return A;
}

bool RealMatrix::is_symmetric() const {
    if (rows_ != cols_)
        return false;

    for (std::size_t i = 0; i < rows_; i++) {
        for (std::size_t j = i + 1; j < cols_; j++) {
            if (!fpeq(data[i * cols_ + j], data[j * cols_ + i])) 
                return false;
        }
    }

    return true;
}

RealMatrix RealMatrix::transpose() const {
    auto new_data = std::vector(this->cols_ * this->rows_, 0.0);

    for (size_t i = 0; i < this->rows_; ++i) {
        for (size_t j = 0; j < this->cols_; ++j) {
            new_data[j * rows_ + i] = this->unchecked_at(i, j);
        }
    }

    return RealMatrix(std::move(new_data), cols_, rows_);
}

EigenSystem RealMatrix::eig() const {
    if (rows_ != cols_)
        throw std::invalid_argument("Eigen decomposition requires a square matrix.");

    if (!is_symmetric())
        throw std::runtime_error("Matrix is not symmetric.");

    gsl_matrix_sptr M = this->to_gsl_matrix();
    gsl_vector_sptr eval = make_gsl_vector (rows_);
    gsl_matrix_sptr evec = make_gsl_matrix (rows_, rows_);
    gsl_eigen_workspace_sptr w = make_eigen_workspace(rows_);
    gsl_eigen_symmv (M.get(), eval.get(), evec.get(), w.get());
    gsl_eigen_symmv_sort(eval.get(), evec.get(), GSL_EIGEN_SORT_VAL_DESC);

    EigenSystem e;
    e.P = from_gsl_copy(evec.get());
    e.D = diag(eval.get());

    return e;
}

SignedLogDet RealMatrix::slogdet() const {
    if (rows_ != cols_)
        throw std::invalid_argument("LU decomposition requires a square matrix.");

    gsl_matrix_sptr M = this->to_gsl_matrix();
    gsl_permutation_sptr p = make_gsl_permutation(rows_);
    int p_sign;

    if (gsl_linalg_LU_decomp(M.get(), p.get(), &p_sign) != GSL_SUCCESS) {
        return {0.0, 1};
    }
        
    SignedLogDet sld;
    sld.logdet = gsl_linalg_LU_lndet(M.get());
    sld.sign = gsl_linalg_LU_sgndet(M.get(), p_sign);

    return sld;
}

RealMatrix RealMatrix::inv() const {
    if (rows_ != cols_)
        throw std::invalid_argument("inversion requires a square matrix.");

    gsl_matrix_sptr M = this->to_gsl_matrix();
    gsl_matrix_sptr M_inv = make_gsl_matrix(rows_, cols_);
    gsl_permutation_sptr p = make_gsl_permutation(rows_);
    int p_sign;

    if (gsl_linalg_LU_decomp(M.get(), p.get(), &p_sign) != GSL_SUCCESS) {
        throw std::runtime_error("Matrix is singular");
    }
        
    gsl_linalg_LU_invert(M.get(), p.get(), M_inv.get());

    return from_gsl_copy(M_inv.get());
}

RealMatrix RealMatrix::operator-() const {
    auto new_data = std::vector(this->rows_ * this->cols_, 0.0);

    for (size_t i = 0; i < this->rows_; ++i) {
        for (size_t j = 0; j < this->cols_; ++j) {
            new_data[i * cols_ + j] = -this->unchecked_at(i, j);
        }
    }

    return RealMatrix(std::move(new_data), rows_, cols_);
}

RealMatrix &RealMatrix::operator+=(const RealMatrix &rhs) {
    if (this->rows_ != rhs.rows() || this->cols_ != rhs.cols())
        throw std::invalid_argument("Matrices should have the same shape.");

    for (size_t i = 0; i < this->rows_; ++i) {
        for (size_t j = 0; j < this->cols_; ++j) {
            this->unchecked_at(i, j) += rhs.unchecked_at(i, j);
        }
    }

    return *this;
}

RealMatrix &RealMatrix::operator-=(const RealMatrix &rhs) {
    if (this->rows_ != rhs.rows() || this->cols_ != rhs.cols())
        throw std::invalid_argument("Matrices should have the same shape.");

    for (size_t i = 0; i < this->rows_; ++i) {
        for (size_t j = 0; j < this->cols_; ++j) {
            this->unchecked_at(i, j) -= rhs.unchecked_at(i, j);
        }
    }

    return *this;
}

RealMatrix &RealMatrix::operator*=(const RealMatrix &rhs) {
    if (this->cols_ != rhs.rows()) 
        throw std::invalid_argument("Matrices don't have the right shape to be multiplied.");

    const size_t M = rows_;
    const size_t N = rhs.cols_;
    const size_t K = cols_;

    std::vector<double> result(M * N, 0.0);

    constexpr size_t BLOCK = 32;

    for (size_t ii = 0; ii < M; ii += BLOCK) {
        for (size_t kk = 0; kk < K; kk += BLOCK) {
            for (size_t jj = 0; jj < N; jj += BLOCK) {

                size_t i_max = std::min(ii + BLOCK, M);
                size_t k_max = std::min(kk + BLOCK, K);
                size_t j_max = std::min(jj + BLOCK, N);

                for (size_t i = ii; i < i_max; ++i) {
                    for (size_t k = kk; k < k_max; ++k) {
                        double aik = unchecked_at(i, k);
                        const double* rhs_row = &rhs.data[k * rhs.cols_];

                        double* res_row = &result[i * N];

                        for (size_t j = jj; j < j_max; ++j) {
                            res_row[j] += aik * rhs_row[j];
                        }
                    }
                }
            }
        }
    }

    data = std::move(result);
    cols_ = N;
    return *this;
}

RealMatrix &RealMatrix::operator*=(double scalar) {
    for (double& x: data)
        x *= scalar;

    return *this;
}

RealMatrix &RealMatrix::operator/=(double scalar) {
    if (fpeq(scalar, 0.0))
        throw std::invalid_argument("Division by zero.");

    for (double& x: data)
        x /= scalar;

    return *this;
}

double &RealMatrix::unchecked_at(size_t i, size_t j) {
    return data[i * cols_ + j];
}

const double &RealMatrix::unchecked_at(size_t i, size_t j) const {
    return data[i * cols_ + j];
}
