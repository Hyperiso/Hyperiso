#ifndef LIKELIHOOD_DISCRETE_H
#define LIKELIHOOD_DISCRETE_H

#include <random>
#include <vector>
#include <numeric>
#include <stdexcept>
#include <cmath>

#include "IMarginalDistribution.h"
#include "AbstractConfig.h"

/**
 * @file LikelihoodDiscrete.h
 * @brief Discrete marginal distribution built from weighted support points.
 *
 * This header defines:
 * - @ref LikelihoodMarginalCfg: configuration object carrying support values
 *   and associated weights,
 * - @ref LikelihoodMarginal: discrete sampler based on Vose's alias method.
 *
 * The distribution is defined by a finite support:
 * \f[
 *   X \in \{x_1,\dots,x_n\}, \qquad P(X=x_i)\propto w_i.
 * \f]
 *
 * This implementation is primarily intended for efficient sampling from
 * empirical or likelihood-derived discrete supports.
 */

/**
 * @struct LikelihoodMarginalCfg
 * @brief Configuration object for @ref LikelihoodMarginal.
 *
 * The vectors @ref values and @ref weights must have the same size, and
 * weights are interpreted up to an overall normalization factor.
 */
struct LikelihoodMarginalCfg : public AbstractConfig {
    Vector values;  /// Support points of the discrete distribution.
    Vector weights; /// Non-negative weights associated with each support point.
};

/**
 * @class LikelihoodMarginal
 * @brief Discrete marginal sampler using Vose's alias method.
 *
 * This class implements efficient O(1) sampling from a finite weighted support.
 *
 * The constructor:
 * - checks consistency of values/weights,
 * - optionally computes weighted mean and standard deviation,
 * - builds alias tables for fast random generation.
 *
 * If @p standardize is enabled, samples returned by @ref rvs are transformed as:
 * \f[
 *   z = \frac{x-\mu}{\sigma}.
 * \f]
 *
 * @note In the current implementation:
 * - @ref rvs is fully implemented,
 * - @ref logpdf, @ref cdf, @ref ppf, @ref mean and @ref std return placeholder
 *   values rather than the exact discrete-distribution quantities.
 *   This matches the current code but should be documented clearly for users.
 */
class LikelihoodMarginal final : public IMarginalDistribution {
public:
    /**
     * @brief Constructs a discrete likelihood marginal.
     *
     * @param values Support values.
     * @param weights Non-negative weights associated with support values.
     * @param seed Seed for the internal random engine.
     * @param standardize If true, sampled values are centered and rescaled
     *        using the weighted mean and standard deviation.
     *
     * @throws std::invalid_argument if:
     * - values is empty,
     * - values and weights sizes differ,
     * - a weight is invalid,
     * - the total weight is non-positive.
     */
    LikelihoodMarginal(Vector values,
                       Vector weights,
                       unsigned int seed = std::random_device{}(),
                       bool standardize = false);
    
    /// \copydoc IMarginalDistribution::rvs
    Vector rvs(std::size_t n) override;

    /**
     * @brief Placeholder log-density implementation.
     *
     * @note The current implementation always returns 0.0.
     *       This class is currently intended mainly as a sampler.
     */
    double logpdf(double x) override { return 0.0; }

    /**
     * @brief Placeholder CDF implementation.
     *
     * @note The current implementation always returns 0.0.
     */

    double cdf(double x) override { return 0.0; }

    /**
     * @brief Placeholder inverse-CDF implementation.
     *
     * @note The current implementation always returns 0.0.
     */
    double ppf(double p) override { return 0.0; }

    /**
     * @brief Placeholder mean implementation.
     *
     * @note The current implementation always returns 0.0, even when the
     *       internal weighted mean has been computed for standardization.
     */
    double mean() override { return 0.0; }

    /**
     * @brief Placeholder standard deviation implementation.
     *
     * @note The current implementation always returns 0.0, even when the
     *       internal weighted standard deviation has been computed.
     */
    double std() override { return 0.0; }

private:
    /**
     * @brief Builds alias tables for O(1) discrete sampling.
     *
     * @param weights Raw non-negative weights.
     */
    void build_alias_tables(std::vector<double> weights);

    std::mt19937 eng_;                              /// Underlying pseudo-random engine.
    std::uniform_real_distribution<double> u01_;    /// Uniform real distribution on [0,1).
    std::vector<double> values_;                    /// Discrete support values.

    std::vector<double> prob_;                      /// Alias probability table.
    std::vector<std::size_t> alias_;                /// Alias indirection table.

    bool standardize_;                              /// Whether returned samples should be standardized.
    double mean_ = 0.0;                             /// Weighted mean used when standardization is enabled.
    double std_  = 1.0;                             /// Weighted standard deviation used when standardization is enabled.
};

#endif
