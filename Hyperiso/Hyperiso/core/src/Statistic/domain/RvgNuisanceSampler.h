#ifndef RVG_NUISANCE_SAMPLER_H
#define RVG_NUISANCE_SAMPLER_H

#include <vector>
#include <random>
#include <stdexcept>
#include <cmath>

#include "INuisanceSampler.h"
#include "JointDistribution.h"
#include "MarginalType.h"
#include "MarginalFactory.h"
#include "CopulaFactory.h"
#include "Indexing.h"

/**
 * @file RvgNuisanceSampler.h
 * @brief Nuisance sampler backed by a joint random-vector generator.
 *
 * @see INuisanceSampler
 * @see JointDistribution
 */

/**
 * @class RvgNuisanceSampler
 * @brief Draws nuisance-parameter maps from a @ref JointDistribution.
 *
 * The sampler draws vectors from the provided joint distribution and associates
 * each component with the corresponding @ref ParamId stored in construction
 * order.
 */
class RvgNuisanceSampler final : public INuisanceSampler {
public:
    /**
     * @brief Constructs a sampler from parameter identifiers and a joint law.
     *
     * @param ids Parameter identifiers defining the order of the sampled vector.
     * @param rvg Joint distribution used to generate nuisance values.
     */
    explicit RvgNuisanceSampler(const std::vector<ParamId>& ids, std::unique_ptr<JointDistribution> rvg);

    /**
     * @copydoc INuisanceSampler::sample()
     */
    std::map<ParamId, double> sample() const override;

    /**
     * @copydoc INuisanceSampler::sample(std::size_t)
     */
    std::vector<std::map<ParamId, double>> sample(std::size_t n) const override;

private:
    std::unique_ptr<JointDistribution> rvg_;    ///< Joint distribution used for random draws.
    std::vector<ParamId> ids_;                  ///< Parameter identifiers matching the distribution coordinates.
};

#endif