#ifndef ILIKELIHOOD_H
#define ILIKELIHOOD_H

#include "Math.h"

/**
 * @file ILikelihood.h
 * @brief Abstract interface for likelihood functions used by the fitting layer.
 *
 * This header defines the minimal contract required from any likelihood object:
 * evaluating a negative log-likelihood, exposing its parameter definitions and
 * reporting the total dimension of the optimization vector.
 *
 * @see IProfileableLikelihood
 * @see BaseLikelihood
 */

/**
 * @class ILikelihood
 * @brief Base interface for negative log-likelihood evaluators.
 *
 * Implementations receive a parameter vector \f$\theta\f$ and return the
 * corresponding negative log-likelihood value. The interface deliberately keeps
 * the parameter layout abstract; concrete implementations are responsible for
 * defining how \f$\theta\f$ is split and interpreted.
 */
class ILikelihood {
public:
    /**
     * @brief Virtual destructor for safe polymorphic deletion.
     */
    virtual ~ILikelihood() = default;

    /**
     * @brief Evaluates the negative log-likelihood at a given parameter point.
     *
     * @param theta Parameter vector evaluated by the likelihood implementation.
     * @return Negative log-likelihood value associated with @p theta.
     */
    virtual double nll(const std::vector<double>& theta) const = 0;

    /**
     * @brief Returns the metadata describing the likelihood parameters.
     *
     * The returned definitions are expected to follow the same ordering as the
     * parameter vector accepted by @ref nll.
     *
     * @return Ordered list of parameter definitions.
     */
    virtual std::vector<fit_app::ParameterDefinition> get_param_defs() const = 0;

    /**
     * @brief Returns the total dimension of the parameter vector.
     *
     * @return Number of scalar parameters expected by @ref nll.
     */
    virtual std::size_t dim() const = 0;
};

#endif // ILIKELIHOOD_H
