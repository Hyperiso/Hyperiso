#ifndef IPROFILEABLE_LIKELIHOOD_H
#define IPROFILEABLE_LIKELIHOOD_H

#include "ILikelihood.h"

/**
 * @file IProfileableLikelihood.h
 * @brief Interface for likelihoods that separate fitted and nuisance parameters.
 *
 * A profileable likelihood exposes the parameter vector as two blocks:
 * - physics or fitted parameters \f$p\f$,
 * - nuisance parameters \f$\eta\f$.
 *
 * This separation is useful for profile-likelihood scans, uncertainty
 * propagation, curvature-based approximations and diagnostic tooling.
 *
 * @see ILikelihood
 * @see BaseLikelihood
 */

/**
 * @class IProfileableLikelihood
 * @brief Extension of @ref ILikelihood with explicit parameter-block access.
 *
 * Implementations must provide utilities to evaluate the likelihood from split
 * parameter blocks, compute model predictions and residuals, and expose local
 * curvature matrices for the observable and nuisance terms.
 */
class IProfileableLikelihood : public ILikelihood {
public:
    /**
     * @brief Returns the dimension of the fitted-parameter block.
     *
     * @return Number of scalar parameters in \f$p\f$.
     */
    virtual std::size_t p_dimension() const = 0;

    /**
     * @brief Returns the dimension of the nuisance-parameter block.
     *
     * @return Number of scalar parameters in \f$\eta\f$.
     */
    virtual std::size_t eta_dimension() const = 0;

    /**
     * @brief Returns the central values of the fitted parameters.
     *
     * @return Vector containing the nominal values of \f$p\f$.
     */
    virtual std::vector<double> central_p() const = 0;

    /**
     * @brief Returns the central values of the nuisance parameters.
     *
     * @return Vector containing the nominal values of \f$\eta\f$.
     */
    virtual std::vector<double> central_eta() const = 0;

    /**
     * @brief Evaluates the model prediction for split parameters.
     *
     * @param p   Fitted-parameter block.
     * @param eta Nuisance-parameter block.
     * @return Model prediction vector associated with @p p and @p eta.
     */
    virtual std::vector<double> predict(
        const std::vector<double>& p,
        const std::vector<double>& eta
    ) const = 0;

    /**
     * @brief Computes observable residuals for split parameters.
     *
     * Residuals are generally defined as model predictions minus the experimental
     * observable central values.
     *
     * @param p   Fitted-parameter block.
     * @param eta Nuisance-parameter block.
     * @return Residual vector in observable space.
     */
    virtual std::vector<double> residuals(
        const std::vector<double>& p,
        const std::vector<double>& eta
    ) const = 0;

    /**
     * @brief Evaluates the negative log-likelihood from split parameters.
     *
     * @param p   Fitted-parameter block.
     * @param eta Nuisance-parameter block.
     * @return Negative log-likelihood value at \f$(p, \eta)\f$.
     */
    virtual double nll_from_split(
        const std::vector<double>& p,
        const std::vector<double>& eta
    ) const = 0;

    /**
     * @brief Computes the observable-term curvature matrix.
     *
     * @param residuals Residual vector in observable space.
     * @return Curvature matrix associated with the observable likelihood term.
     */
    virtual RealMatrix observable_curvature(
        const std::vector<double>& residuals
    ) const = 0;

    /**
     * @brief Computes the nuisance-term curvature matrix.
     *
     * @param eta Nuisance-parameter vector.
     * @return Curvature matrix associated with the nuisance likelihood term.
     */
    virtual RealMatrix nuisance_curvature(
        const std::vector<double>& eta
    ) const = 0;
};

#endif