#ifndef GAUSSIANCOPULA_H
#define GAUSSIANCOPULA_H

#include "ICopula.h"
#include "GenericCopula.h"
#include "Matrix.h"
#include "AbstractConfig.h"

/**
 * @file GaussianCopula.h
 * @brief Gaussian copula implementation.
 *
 * This header defines:
 * - @ref GaussianCopulaConfig, a lightweight config object holding the
 *   correlation matrix,
 * - @ref GaussianCopula, a Gaussian copula based on a correlation matrix.
 *
 * The Gaussian copula is built from a multivariate normal latent vector:
 * - sample \f$Z \sim \mathcal{N}(0, R)\f$,
 * - transform component-wise via the standard normal CDF,
 * - obtain \f$U_i = \Phi(Z_i)\f$.
 *
 * @see ICopula
 * @see GenericCopula
 * @see StudentTCopula
 */

/**
 * @struct GaussianCopulaConfig
 * @brief Configuration object for a Gaussian copula.
 *
 * The matrix @ref R is expected to represent a correlation matrix
 * (square, symmetric, unit diagonal, positive semidefinite ideally).
 *
 * In practice, the implementation projects it to the nearest PSD matrix
 * before use.
 */
struct GaussianCopulaConfig : public AbstractConfig {
    RealMatrix R;   /// Correlation matrix driving the Gaussian dependence structure.
};

/**
 * @class GaussianCopula
 * @brief Gaussian copula based on a correlation matrix.
 *
 * Internally, the copula stores:
 * - a correlation matrix @ref R,
 * - its inverse,
 * - its Cholesky factor,
 * - and its log-determinant.
 *
 * Sampling:
 * - draw a standard Gaussian vector,
 * - correlate it with the Cholesky factor,
 * - map each component through the Gaussian CDF.
 *
 * Density evaluation:
 * - transform uniforms back to Gaussian latent variables,
 * - evaluate the standard Gaussian copula log-density.
 *
 * @note The input matrix is regularized through @ref nearest_psd before use.
 */
class GaussianCopula final : public GenericCopula {
public:
    /**
     * @brief Constructs a Gaussian copula from a correlation matrix.
     *
     * The matrix is projected to a nearest positive semidefinite correlation
     * matrix before decomposition.
     *
     * @param seed Seed for the RNG.
     * @param R    Correlation matrix.
     */
    explicit GaussianCopula(unsigned int seed, RealMatrix R);

    /**
     * @copydoc ICopula::sample_u(std::size_t)
     */
    std::vector<Vector> sample_u(std::size_t n) override;

    /**
     * @copydoc ICopula::sample_u()
     */
    Vector sample_u() override;

    /**
     * @copydoc ICopula::log_density(Vector)
     *
     * This evaluates the Gaussian copula density:
     * \f[
     *   \log c(u)
     *   = -\frac12 \log\det(R)
     *     - \frac12 z^\top (R^{-1} - I) z
     * \f]
     * where \f$z_i = \Phi^{-1}(u_i)\f$.
     */
    double log_density(Vector u) override;
    RealMatrix dlog_density(std::vector<double> u) override;
    RealMatrix ddlog_density(std::vector<double> u) override;
    LogDensityDiff log_c_dc_ddc(std::vector<double> u) override;


private:
    
    RealMatrix R;       /// Correlation matrix used by the copula.
    RealMatrix R_inv;   /// Inverse correlation matrix.
    RealMatrix L;       /// Lower Cholesky factor such that \f$R = LL^\top\f$.
    double logdet;      /// Log-determinant of the correlation matrix.
};

#endif // GAUSSIANCOPULA_H
