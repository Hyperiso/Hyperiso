#ifndef SPLITGAUSSIANMARGINAL_H
#define SPLITGAUSSIANMARGINAL_H

#include <random>

#include "IMarginalDistribution.h"
#include "Math.h"
#include "AbstractConfig.h"

/**
 * @file SplitGaussianMarginal.h
 * @brief Asymmetric split-Gaussian marginal distribution.
 *
 * This header defines:
 * - @ref SplitGaussianMarginalCfg: configuration object,
 * - @ref SplitGaussianMarginal: asymmetric Gaussian distribution with
 *   different widths on each side of the mode.
 *
 * The split Gaussian is useful to model asymmetric uncertainties, with:
 * - left width  @f$\sigma_-\f$,
 * - right width @f$\sigma_+\f$,
 * around a central location @f$\mu\f$.
 *
 * @see IMarginalDistribution
 */

using gsl_rng_sptr = std::unique_ptr<gsl_rng, decltype(&gsl_rng_free)>;

/**
 * @struct SplitGaussianMarginalCfg
 * @brief Configuration object for @ref SplitGaussianMarginal.
 */
struct SplitGaussianMarginalCfg : public AbstractConfig {
    double mu {0.0};        /// Central value / mode.
    double sigma_p {1.0};   /// Standard deviation on the right side (x > mu).
    double sigma_m {1.0};   /// Standard deviation on the left side (x <= mu).
};

/**
 * @class SplitGaussianMarginal
 * @brief Piecewise Gaussian marginal with asymmetric widths.
 *
 * This class represents a split normal / two-piece Gaussian distribution:
 * - for x <= mu, the width is sigma_m,
 * - for x >  mu, the width is sigma_p.
 *
 * The distribution is normalized through:
 * - @ref N : normalization constant,
 * - @ref w : left/right probability split used in the CDF/PPF.
 *
 * Sampling is currently implemented by inverse transform:
 * a uniform random number is drawn and passed through @ref ppf.
 */
class SplitGaussianMarginal final : public IMarginalDistribution {
public:
    /**
     * @brief Constructs an asymmetric split-Gaussian marginal.
     *
     * @param mu Central value / mode.
     * @param sigma_p Right-side width.
     * @param sigma_m Left-side width.
     * @param seed Seed for the internal RNG.
     */
    explicit SplitGaussianMarginal(double mu, double sigma_p, double sigma_m, unsigned int seed = std::random_device{}());

    /// \copydoc IMarginalDistribution::rvs
    Vector rvs(std::size_t n) override;

    /// \copydoc IMarginalDistribution::logpdf
    double logpdf(double x) override;

    /// \copydoc IMarginalDistribution::cdf
    double cdf(double x) override;

    /// \copydoc IMarginalDistribution::ppf
    double ppf(double p) override;

    /// \copydoc IMarginalDistribution::mean
    double mean() override;

    /// \copydoc IMarginalDistribution::std
    double std() override;

private:
    double mu;                                      /// Central value / mode.
    double sigma_p;                                 /// Right-side standard deviation.
    double sigma_m;                                 /// Left-side standard deviation.
    double N;                                       /// Overall normalization factor.
    double w;                                       /// Weight controlling the left/right CDF split.
    const gsl_rng_type* rng_tp {gsl_rng_mt19937};   /// GSL RNG type used for sampling.
    gsl_rng_sptr eng_{gsl_rng_alloc(rng_tp), &gsl_rng_free};      /// Owned GSL random-number generator.
};

#endif // SPLITGAUSSIANMARGINAL_H
