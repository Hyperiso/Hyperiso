#ifndef DISTRIBUTION_FACTORY_H
#define DISTRIBUTION_FACTORY_H

#include <memory>
#include <random>
#include <algorithm>
#include <variant>

#include "MarginalType.h"
#include "IMarginalDistribution.h"
#include "GaussianMarginal.h"
#include "FlatMarginal.h"
#include "LikelihoodMarginal.h"
#include "SplitGaussianMarginal.h"

/**
 * @file DistributionFactory.h
 * @brief Factory for instantiating concrete marginal distributions.
 *
 * This header defines:
 * - @ref MarginalConfig, the variant of supported distribution configs,
 * - @ref DistributionFactory, which materializes a concrete
 *   @ref IMarginalDistribution from a @ref MarginalType and a compatible
 *   configuration object.
 *
 * This is the second stage after @ref MarginalConfigFactory:
 * - first, build the configuration from physics/statistical metadata,
 * - then, build the concrete random distribution object.
 *
 * @see IMarginalDistribution
 * @see MarginalConfigFactory
 */

/**
 * @typedef MarginalConfig
 * @brief Variant of all supported marginal configuration objects.
 *
 * The contained type must be compatible with the requested @ref MarginalType.
 */
using MarginalConfig = std::variant<FlatMarginalCfg, GaussianMarginalCfg, SplitGaussianMarginalCfg, LikelihoodMarginalCfg>;

/**
 * @class DistributionFactory
 * @brief Static factory creating scalar marginal-distribution objects.
 *
 * This factory converts:
 * - a distribution tag (@ref MarginalType),
 * - a typed configuration object (@ref MarginalConfig),
 * - and an RNG seed,
 * into a concrete @ref IMarginalDistribution implementation.
 *
 * Internally, construction is delegated through a `std::visit` over the config
 * variant, so the actual type stored in @ref MarginalConfig matters.
 *
 * @note The current implementation accepts the following effective mappings:
 * - GAUSSIAN      -> usually @ref GaussianMarginal
 * - HALF_GAUSSIAN -> currently relies on @ref SplitGaussianMarginalCfg
 *                    / @ref SplitGaussianMarginal when provided
 * - FLAT          -> @ref FlatMarginal
 *
 * @warning In the current implementation, the LIKELIHOOD branch throws
 * before reaching the generic variant-based constructor, so it is effectively
 * disabled through this API entry point.
 */
class DistributionFactory {
public:
    /**
     * @brief Creates a concrete marginal-distribution object.
     *
     * The returned distribution is heap-allocated and owned through
     * `std::unique_ptr<IMarginalDistribution>`.
     *
     * @param name Requested marginal-distribution family.
     * @param cfg  Variant containing the concrete configuration data.
     * @param seed Seed passed to the underlying distribution constructor.
     * @return Unique pointer to a newly created marginal distribution.
     *
     * @throws std::invalid_argument if:
     * - the marginal type is unknown,
     * - the current implementation rejects the selected type
     *   (notably LIKELIHOOD through this entry point).
     */
    static std::unique_ptr<IMarginalDistribution> create(MarginalType name, MarginalConfig cfg, unsigned int seed = std::random_device{}());
};

#endif