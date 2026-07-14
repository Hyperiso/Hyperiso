#ifndef MARGINALCONFIGFACTORY_H
#define MARGINALCONFIGFACTORY_H

#include <variant>

#include "Include.h"
#include "MarginalType.h"
#include "StatParameterProxy.h"
#include "GaussianMarginal.h"
#include "FlatMarginal.h"
#include "SplitGaussianMarginal.h"
#include "LikelihoodMarginal.h"
#include "ExperimentObs.h"
#include "NuisanceSpec.h"

/**
 * @file MarginalConfigFactory.h
 * @brief Factory for building marginal-distribution configuration objects.
 *
 * This header defines:
 * - @ref MarginalConfig, a tagged union of all supported marginal
 *   configuration structs,
 * - @ref MarginalConfigFactory, which derives the appropriate configuration
 *   from an existing parameter or observable and a requested
 *   @ref MarginalType.
 *
 * The goal of this factory is to separate:
 * - extraction of numerical inputs (central values, uncertainties, etc.),
 * - from actual construction of concrete distributions
 *   (handled later by @ref DistributionFactory).
 *
 * Typical workflow:
 * @code
 *   MarginalConfigFactory cfg_factory;
 *   MarginalConfig cfg = cfg_factory.create(pid, MarginalType::GAUSSIAN);
 * @endcode
 *
 * @see DistributionFactory
 * @see StatParameterProxy
 */

/**
 * @typedef MarginalConfig
 * @brief Variant holding any supported marginal configuration object.
 *
 * The active alternative depends on the requested @ref MarginalType and on
 * the information available for the targeted parameter / observable.
 */
using MarginalConfig = std::variant<FlatMarginalCfg, GaussianMarginalCfg, SplitGaussianMarginalCfg, LikelihoodMarginalCfg>;

/**
 * @class MarginalConfigFactory
 * @brief Builds distribution-configuration objects from parameters or observables.
 *
 * This factory converts a statistical object description into a
 * configuration suitable for a concrete marginal distribution.
 *
 * Current implemented policies:
 * - For a @ref ParamId:
 *   - GAUSSIAN  -> centered at parameter value, width = combined uncertainty
 *   - FLAT      -> interval chosen so that the flat std matches the combined uncertainty
 * - For a @ref BinnedObservableId:
 *   - GAUSSIAN  -> mean fixed to 0, width = combined observable uncertainty
 *   - FLAT      -> symmetric interval around 0 with matching std
 *
 * Current non-implemented cases:
 * - HALF_GAUSSIAN
 * - LIKELIHOOD
 *
 * Those currently raise exceptions in the implementation.
 */
class MarginalConfigFactory {
public:
    /**
     * @brief Builds a marginal configuration from a parameter identifier.
     *
     * The parameter is queried through an internal @ref StatParameterProxy.
     *
     * Current behavior:
     * - @ref MarginalType::GAUSSIAN returns @ref GaussianMarginalCfg{mu, sigma}
     * - @ref MarginalType::FLAT returns @ref FlatMarginalCfg with bounds
     *   \f$\mu \pm \sigma\sqrt{3}\f$, so that the flat distribution has
     *   standard deviation \f$\sigma\f$
     *
     * @param pid      Parameter identifier.
     * @param marginal Requested marginal family.
     * @return A @ref MarginalConfig variant containing the matching config type.
     *
     * @throws std::runtime_error for currently unsupported marginals.
     * @throws std::invalid_argument if the marginal type is unknown.
     */
    MarginalConfig create(ParamId pid, MarginalType marginal);

    /**
     * @brief Builds a parameter marginal while honoring an explicit nuisance specification.
     *
     * For a flat nuisance, the configured support is taken directly from
     * @p spec.bounds instead of being reconstructed from the parameter standard
     * deviation. Other marginal families retain the standard parameter-derived
     * configuration.
     *
     * @param pid Parameter identifier.
     * @param marginal Effective marginal family after default/user resolution.
     * @param spec Effective nuisance specification associated with @p pid.
     * @return A marginal configuration compatible with @ref MarginalFactory.
     *
     * @throws std::invalid_argument if the specification block/code does not match @p pid,
     *         or if flat bounds are non-finite or not strictly ordered.
     */
    MarginalConfig create(ParamId pid, MarginalType marginal, const NuisanceSpec& spec);

    /**
     * @brief Builds a vector of marginal configurations from a binned observable identifier.
     *
     * Current convention for observables:
     * - the marginal is centered at 0,
     * - the spread is derived from the observable combined uncertainty.
     *
     * This is useful when the marginal describes a fluctuation / nuisance-like
     * variable rather than an absolute physical central value.
     *
     * @param pid      ExperimentObs observable identifier.
     * @param marginal Requested marginal family.
     * @return A @ref MarginalConfig variant containing the matching config type.
     *
     * @throws std::runtime_error for currently unsupported marginals.
     * @throws std::invalid_argument if the marginal type is unknown.
     */
    MarginalConfig create(ExperimentObs pid, MarginalType marginal);

private:
    /// Statistical parameter proxy used to fetch central values and uncertainties.
    const StatParameterProxy p {};
};

#endif // MARGINALCONFIGFACTORY_H
