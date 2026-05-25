#ifndef WITHGAUSSIANCONSTRAINTS_H
#define WITHGAUSSIANCONSTRAINTS_H

#include "ILikelihood.h"
#include "JointDistribution.h"

/**
 * @file WithGaussianConstraints.h
 * @brief Likelihood decorator adding an external Gaussian constraint term.
 *
 * This header defines @ref WithGaussianConstraints, a lightweight wrapper around
 * an existing likelihood. The wrapper evaluates the base likelihood and adds the
 * negative log-density of a joint constraint distribution applied to a selected
 * subset of the parameter vector.
 *
 * @see ILikelihood
 * @see JointDistribution
 */

/**
 * @class WithGaussianConstraints
 * @brief Decorates a likelihood with Gaussian constraints on selected parameters.
 *
 * The class is useful when a profile or projection should keep some parameters
 * weakly constrained around a fitted central value. For a full parameter vector
 * \f$\theta\f$, the objective evaluated by this wrapper is
 * \f[
 *   \mathrm{NLL}_{\mathrm{total}}(\theta)
 *   = \mathrm{NLL}_{\mathrm{base}}(\theta)
 *     - \log p_{\mathrm{constraint}}(\theta_{I}),
 * \f]
 * where \f$\theta_I\f$ is the sub-vector selected by @ref constrained_params.
 *
 * @note The wrapper does not own or modify parameter definitions. It forwards
 *       metadata and dimensionality queries to the wrapped likelihood.
 */
class WithGaussianConstraints : public ILikelihood {
public:
    /**
     * @brief Constructs a constrained likelihood wrapper.
     *
     * @param base Base likelihood to evaluate first.
     * @param constraints_dist Joint distribution defining the Gaussian prior-like
     *                         constraint on the selected parameters.
     * @param constrained_params Indices of the full parameter-vector components passed to
     *                           \p constraints_dist, in distribution order.
     */
    WithGaussianConstraints(std::shared_ptr<ILikelihood> base, std::shared_ptr<JointDistribution> constraints_dist, std::vector<std::size_t> constrained_params);

    /**
     * @brief Evaluates the constrained negative log-likelihood.
     *
     * The returned value is the sum of the base NLL and the constraint penalty
     * \f$-\log p_{\mathrm{constraint}}\f$ evaluated on the constrained subset.
     *
     * @param theta Full parameter vector passed to the base likelihood.
     * @return Constrained negative log-likelihood value.
     *
     * @throws std::out_of_range if a constrained parameter index is outside
     *         \p theta.
     * @throws std::exception forwarded from the base likelihood or constraint
     *         distribution.
     */
    double nll(const std::vector<double>& theta) const override;

    /**
     * @copydoc ILikelihood::get_param_defs()
     */
    std::vector<fit_app::ParameterDefinition> get_param_defs() const override;

    /**
     * @copydoc ILikelihood::dim()
     */
    std::size_t dim() const override;

private:
    std::shared_ptr<ILikelihood> base;                      ///< Wrapped likelihood.
    std::shared_ptr<JointDistribution> constraints_dist;    ///< Constraint distribution for selected parameters.
    std::vector<std::size_t> constrained_params;            ///< Indices constrained by @ref constraints_dist.
};

#endif // WITHGAUSSIANCONSTRAINTS_H
