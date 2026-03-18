#ifndef DISTRIBUTION_TYPE_H
#define DISTRIBUTION_TYPE_H

/**
 * @file MarginalType.h
 * @brief Enumeration of supported one-dimensional marginal distribution types.
 *
 * This header defines @ref MarginalType, used by the statistical layer to
 * select which scalar probability law should be attached to a parameter
 * or observable.
 *
 * Typical workflow:
 * - a @ref MarginalType is chosen by the user or inference setup,
 * - @ref MarginalConfigFactory builds the appropriate configuration object
 *   from the current parameter/observable metadata,
 * - @ref DistributionFactory materializes the concrete
 *   @ref IMarginalDistribution implementation.
 *
 * @see MarginalConfigFactory
 * @see DistributionFactory
 * @see IMarginalDistribution
 */

/**
 * @enum MarginalType
 * @brief Supported marginal-distribution families.
 *
 * These values identify the probability law to use for a scalar nuisance
 * parameter or observable fluctuation.
 */
enum class MarginalType {
    GAUSSIAN,       ///< Symmetric Gaussian marginal.
    HALF_GAUSSIAN,  ///< Asymmetric / half-Gaussian-like marginal (currently mapped to split Gaussian logic).
    FLAT,           ///< Uniform (flat) marginal on a finite interval.
    LIKELIHOOD      ///< Discrete likelihood-based marginal built from weighted support points.
};

#endif