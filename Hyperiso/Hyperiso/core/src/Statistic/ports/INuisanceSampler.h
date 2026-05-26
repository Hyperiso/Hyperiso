#pragma once
#include <vector>
#include <random>
#include "Matrix.h"
#include "Include.h"


using Vec = std::vector<double>;

/**
 * @file INuisanceSampler.h
 * @brief Interface for nuisance-parameter random samplers.
 *
 * A nuisance sampler draws parameter maps used to propagate theoretical or
 * external-input uncertainties through model predictions.
 */

/**
 * @class INuisanceSampler
 * @brief Abstract source of nuisance-parameter samples.
 *
 * Implementations define the actual sampling law, for example a joint
 * distribution with marginal laws and a copula. Samples are returned as maps
 * keyed by @ref ParamId so they can be passed directly to @ref IModel.
 */
class INuisanceSampler {
public:
    /// Default virtual destructor.
    virtual ~INuisanceSampler() = default;

    /**
     * @brief Draws one nuisance-parameter sample.
     *
     * @return Map from nuisance parameter identifiers to sampled values.
     */
    virtual std::map<ParamId, double> sample() const = 0;

    /**
     * @brief Draws several nuisance-parameter samples.
     *
     * @param n Number of samples to generate.
     *
     * @return Vector of sampled nuisance maps.
     */
    virtual std::vector<std::map<ParamId, double>> sample(std::size_t n) const = 0;
};