#ifndef CHISQUAREDLIKELIHOOD_H
#define CHISQUAREDLIKELIHOOD_H

#include <cmath>
#include <stdexcept>
#include <vector>

#include "BaseLikelihood.h"


/**
 * @file ChiSquaredLikelihood.h
 * @brief Chi-square likelihood implementation for models without explicit nuisance parameters.
 *
 * This file declares @ref ChiSquaredLikelihood, a lightweight likelihood model
 * based on a fixed inverse covariance matrix. It is intended for profiling or
 * minimization problems where the optimization vector contains only the
 * parameters of interest.
 *
 * @see BaseLikelihood
 * @see ILikelihood
 */

/**
 * @class ChiSquaredLikelihood
 * @brief Global chi-square likelihood with no explicit nuisance parameters.
 *
 * The optimization vector is interpreted as the parameter vector \f$p\f$ and
 * the nuisance vector \f$\eta\f$ is always empty. The negative log-likelihood
 * is evaluated as
 *
 * \f[
 *   \mathrm{NLL}(p)
 *   = \frac{1}{2}
 *     \left(f(p, \emptyset) - O_{\mathrm{exp}}\right)^\top
 *     C^{-1}
 *     \left(f(p, \emptyset) - O_{\mathrm{exp}}\right),
 * \f]
 *
 * where \f$f(p, \emptyset)\f$ is the model prediction, \f$O_{\mathrm{exp}}\f$
 * is the vector of experimental observations and \f$C^{-1}\f$ is the supplied
 * inverse covariance matrix.
 *
 * @note This class derives from @ref BaseLikelihood for API compatibility, but
 *       overrides the likelihood evaluation and reports a zero-dimensional
 *       nuisance space.
 */
class ChiSquaredLikelihood final : public BaseLikelihood {
public:
    /**
     * @brief Constructs a chi-square likelihood.
     *
     * @param model Model function mapping \f$(p, \eta)\f$ to observable predictions.
     *              For this likelihood, @p eta is expected to be empty.
     * @param ctx Shared likelihood context containing observations and parameter
     *            definitions.
     * @param p_dim Number of parameters of interest.
     * @param covariance_inv Inverse covariance matrix used in the quadratic form.
     */
    ChiSquaredLikelihood(
        const ModelFn& model,
        std::shared_ptr<LikelihoodContext> ctx,
        std::size_t p_dim,
        RealMatrix covariance_inv
    );

    /**
     * @brief Evaluates the chi-square negative log-likelihood.
     *
     * @param theta Parameter vector. Its size must be exactly @ref p_dim.
     *
     * @return Finite NLL value when the model and covariance are valid; a large
     *         penalty value is returned when non-finite numerical values are
     *         encountered.
     *
     * @throws std::invalid_argument if @p theta does not have dimension @ref p_dim.
     */
    double nll(const std::vector<double>& theta) const override;

    /**
     * @brief Returns the total optimization dimension.
     *
     * For this class, the dimension is equal to the number of parameters of
     * interest because no nuisance parameters are explicitly represented.
     *
     * @return Parameter-space dimension.
     */
    std::size_t dim() const override;

    /**
     * @brief Returns the observable-space curvature matrix.
     *
     * @param r Residual vector. The argument is accepted for interface
     *          compatibility and is not used because the curvature is constant.
     *
     * @return The inverse covariance matrix supplied at construction.
     */
    RealMatrix observable_curvature(
        const std::vector<double>& r
    ) const override;

    /**
     * @brief Returns the nuisance-space curvature matrix.
     *
     * @param eta Nuisance vector. Must be empty for this likelihood.
     *
     * @return An empty \f$0 \times 0\f$ matrix.
     *
     * @throws std::invalid_argument if @p eta is not empty.
     */
    RealMatrix nuisance_curvature(
        const std::vector<double>& eta
    ) const override;

private:
    RealMatrix covariance_inv_;     ///< Inverse covariance matrix used in the chi-square quadratic form.
};

#endif // CHISQUAREDLIKELIHOOD_H
