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
    std::vector<double> values;  /// Support points of the discrete distribution.
    std::vector<double> weights; /// Non-negative weights associated with each support point.
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
 * @note Random sampling is implemented. Analytical distribution methods are
 * deliberately unavailable until their discrete semantics are finalized; they
 * throw `std::logic_error` instead of returning a plausible placeholder value.
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
    LikelihoodMarginal(std::vector<double> values,
                       std::vector<double> weights,
                       unsigned int seed = std::random_device{}(),
                       bool standardize = false);
    
    /// \copydoc IMarginalDistribution::rvs
    std::vector<double> rvs(std::size_t n) override;

    /** @brief Not implemented; use this class only for sampling. */
    double logpdf(double) override {
        throw std::logic_error("LikelihoodMarginal::logpdf is not implemented");
    }

    /** @brief Not implemented; use this class only for sampling. */
    PDFDiff f_df_ddf(double) override {
        throw std::logic_error("LikelihoodMarginal::f_df_ddf is not implemented");
    }

    /** @brief Not implemented; use this class only for sampling. */

    double cdf(double) override {
        throw std::logic_error("LikelihoodMarginal::cdf is not implemented");
    }

    /** @brief Not implemented; use this class only for sampling. */
    double ppf(double) override {
        throw std::logic_error("LikelihoodMarginal::ppf is not implemented");
    }

    /** @brief Not implemented; use this class only for sampling. */
    double mean() override {
        throw std::logic_error("LikelihoodMarginal::mean is not implemented");
    }

    /** @brief Not implemented; use this class only for sampling. */
    double std() override {
        throw std::logic_error("LikelihoodMarginal::std is not implemented");
    }

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
