#ifndef STANDARD_FLAT_H
#define STANDARD_FLAT_H

#include <random>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "IMarginalDistribution.h"
#include "Include.h"
#include "AbstractConfig.h"

/**
 * @file StandardFlat.h
 * @brief Uniform (flat) marginal distribution on a bounded interval.
 *
 * This header defines:
 * - @ref FlatMarginalCfg: lightweight configuration object for a flat marginal,
 * - @ref FlatMarginal: implementation of a continuous uniform distribution
 *   on the interval [a, b].
 *
 * Sampling and distribution functions are implemented using GSL.
 *
 * @see IMarginalDistribution
 */

 using gsl_rng_sptr = std::unique_ptr<gsl_rng, decltype(&gsl_rng_free)>;

/**
 * @struct FlatMarginalCfg
 * @brief Configuration object for @ref FlatMarginal.
 *
 * The interval bounds define the support of the uniform law:
 * \f[
 *   X \sim \mathcal{U}(a,b).
 * \f]
 */
struct FlatMarginalCfg : public AbstractConfig {
    double a;   /// Lower bound of the support.
    double b;   /// Upper bound of the support.

    /// Default constructor.
    FlatMarginalCfg() = default;

    /**
     * @brief Constructs a flat marginal configuration.
     * @param a Lower bound.
     * @param b Upper bound.
     */
    FlatMarginalCfg(double a, double b) : a(a), b(b) {}
};

/**
 * @class FlatMarginal
 * @brief Continuous uniform marginal distribution on [a,b].
 *
 * This class models the flat probability density:
 * \f[
 *   f(x) = \frac{1}{b-a}, \qquad x \in [a,b].
 * \f]
 *
 * It provides:
 * - random variate generation,
 * - log-density,
 * - CDF,
 * - inverse CDF,
 * - analytical mean and standard deviation.
 *
 * A GSL random-number generator is owned internally.
 *
 * The constructor rejects non-finite or reversed bounds. The density and
 * log-density use the closed support convention at the two endpoints; this
 * choice has no effect on probabilities for a continuous distribution.
 */
class FlatMarginal final : public IMarginalDistribution {
public:
    /**
     * @brief Constructs a uniform marginal on [a,b].
     *
     * @param a Lower bound of the support.
     * @param b Upper bound of the support.
     * @param seed Seed for the internal RNG.
     */
    explicit FlatMarginal(double a, double b, unsigned int seed = std::random_device{}());

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

    /// \copydoc IMarginalDistribution::std
    double std() override;

private:
    double a, b;                                    /// Lower and upper bound of the support.
    const gsl_rng_type* rng_tp {gsl_rng_mt19937};   /// Lower bound of the support.
    gsl_rng_sptr eng_{gsl_rng_alloc(rng_tp), &gsl_rng_free};      /// Owned GSL random-number generator.
};

#endif