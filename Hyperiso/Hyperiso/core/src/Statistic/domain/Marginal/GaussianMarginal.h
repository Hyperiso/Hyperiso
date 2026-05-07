#ifndef STANDARD_NORMAL_H
#define STANDARD_NORMAL_H

#include <random>

#include "IMarginalDistribution.h"
#include "AbstractConfig.h"
#include "Math.h"

/**
 * @file StandardNormal.h
 * @brief Gaussian marginal distribution with arbitrary mean and standard deviation.
 *
 * This header defines:
 * - @ref GaussianMarginalCfg: configuration object for Gaussian marginals,
 * - @ref GaussianMarginal: normal distribution implementation backed by GSL.
 *
 * The distribution is:
 * \f[
 *   X \sim \mathcal{N}(\mu,\sigma^2).
 * \f]
 *
 * @see IMarginalDistribution
 */

using gsl_rng_sptr = std::unique_ptr<gsl_rng, decltype(&gsl_rng_free)>;

/**
 * @struct GaussianMarginalCfg
 * @brief Configuration object for @ref GaussianMarginal.
 */
struct GaussianMarginalCfg : public AbstractConfig {
    double mu {0.};     /// Mean of the Gaussian distribution.
    double sigma {1.};  /// Standard deviation of the Gaussian distribution.

    /// Default constructor.
    GaussianMarginalCfg() = default;

    /**
     * @brief Constructs a Gaussian configuration.
     * @param mu Mean.
     * @param sigma Standard deviation.
     */
    GaussianMarginalCfg(double mu, double sigma) : mu(mu), sigma(sigma) {}
};

/**
 * @class GaussianMarginal
 * @brief One-dimensional Gaussian marginal distribution.
 *
 * This class models the normal distribution:
 * \f[
 *   f(x)=\frac{1}{\sqrt{2\pi}\sigma}\exp\left(-\frac{(x-\mu)^2}{2\sigma^2}\right).
 * \f]
 *
 * It provides:
 * - iid sampling,
 * - exact log-density,
 * - CDF,
 * - inverse CDF,
 * - analytical mean and standard deviation.
 *
 * Sampling and distribution utilities are provided through GSL.
 */
class GaussianMarginal final : public IMarginalDistribution {
public:
    /**
     * @brief Constructs a Gaussian marginal.
     *
     * @param mu Mean of the distribution.
     * @param sigma Standard deviation of the distribution.
     * @param seed Seed for the internal RNG.
     */
    explicit GaussianMarginal(double mu, double sigma, unsigned int seed = std::random_device{}());

    /// \copydoc IMarginalDistribution::rvs
    std::vector<double> rvs(std::size_t n) override;

    /// \copydoc IMarginalDistribution::logpdf
    double logpdf(double x) override;

    PDFDiff f_df_ddf(double x) override;

    /// \copydoc IMarginalDistribution::cdf
    double cdf(double x) override;

    /// \copydoc IMarginalDistribution::ppf
    double ppf(double p) override;

    /// \copydoc IMarginalDistribution::mean
    double mean() override;

    /// \copydoc IMarginalDistribution::mean
    double std() override;

private:
    double mu, sigma;                               /// Mean and standard deviation of the Gaussian.
    const gsl_rng_type* rng_tp {gsl_rng_mt19937};   /// GSL RNG type used for sampling.
    gsl_rng_sptr eng_{gsl_rng_alloc(rng_tp), &gsl_rng_free};      /// Owned GSL random-number generator.
};

#endif