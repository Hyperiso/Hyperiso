#ifndef STUDENTTCOPULA_H
#define STUDENTTCOPULA_H

#include <gsl/gsl_sf.h>

#include "GenericCopula.h"
#include "Math.h"
#include "AbstractConfig.h"

/**
 * @file StudentTCopula.h
 * @brief Student-t copula implementation.
 *
 * This header defines:
 * - @ref StudentTCopulaConfig, a lightweight config object holding the
 *   correlation matrix and the number of degrees of freedom,
 * - @ref StudentTCopula, a Student-t copula based on a correlation matrix.
 *
 * The Student-t copula is obtained from a multivariate Student latent vector:
 * - sample \f$Z \sim t_\nu(0, R)\f$,
 * - transform component-wise via the univariate Student-t CDF,
 * - obtain dependent uniforms \f$U_i\f$.
 *
 * Compared with the Gaussian copula, the Student-t copula can model
 * stronger tail dependence.
 *
 * @see ICopula
 * @see GenericCopula
 * @see GaussianCopula
 */

/**
 * @struct StudentTCopulaConfig
 * @brief Configuration object for a Student-t copula.
 */
struct StudentTCopulaConfig : public AbstractConfig {
    RealMatrix R;   /// Correlation matrix of the latent Student-t vector.
    int nu;         /// Degrees of freedom of the Student-t law.
};

/**
 * @class StudentTCopula
 * @brief Student-t copula with correlation matrix and degrees of freedom.
 *
 * Internally, the copula stores:
 * - a PSD-corrected correlation matrix,
 * - its inverse,
 * - its Cholesky factor,
 * - its log-determinant,
 * - and the number of degrees of freedom @ref nu.
 *
 * Sampling:
 * - draw a multivariate Gaussian latent vector,
 * - correlate it,
 * - rescale by a chi-square variable to obtain a Student-t latent vector,
 * - apply the univariate Student-t CDF component-wise.
 *
 * Density evaluation uses the usual Student-t copula formula based on:
 * - the multivariate Student-t density,
 * - divided by the product of marginal Student-t densities.
 */
class StudentTCopula : public GenericCopula {
public:
    /**
     * @brief Constructs a Student-t copula.
     *
     * The correlation matrix is projected to the nearest PSD correlation matrix.
     *
     * @param seed Seed for the RNG.
     * @param R    Correlation matrix.
     * @param nu   Degrees of freedom. Must be at least 2.
     *
     * @throws std::invalid_argument if @p nu < 2.
     */
    explicit StudentTCopula(unsigned int seed, RealMatrix R, int nu);

    /**
     * @copydoc ICopula::sample_u(std::size_t)
     */
    std::vector<std::vector<double>> sample_u(std::size_t n) override;

    /**
     * @copydoc ICopula::sample_u()
     */
    std::vector<double> sample_u() override;

    /**
     * @copydoc ICopula::log_density(std::vector<double>)
     *
     * The implementation:
     * - transforms \f$u_i\f$ to Student-t quantiles,
     * - evaluates the multivariate Student-t log-density,
     * - subtracts the sum of marginal log-densities,
     * yielding the copula log-density.
     */
    double log_density(std::vector<double> u) override;

private:
    int nu;                 /// Degrees of freedom of the latent Student-t vector.
    RealMatrix R;           /// Correlation matrix used by the copula.
    RealMatrix R_inv;       /// Inverse correlation matrix.
    RealMatrix L;           /// Lower Cholesky factor of @ref R.
    double logdet;          /// Log-determinant of the correlation matrix.
};

#endif // STUDENTTCOPULA_H
