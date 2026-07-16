#ifndef BASELIKELIHOOD_H
#define BASELIKELIHOOD_H

#include <memory>
#include <vector>

#include "IProfileableLikelihood.h"
#include "JointDistribution.h"
#include "Math.h"

/**
 * @file BaseLikelihood.h
 * @brief Concrete profileable likelihood built from a model and joint distributions.
 *
 * The likelihood combines:
 * - a user-provided model function mapping fitted and nuisance parameters to
 *   observable predictions,
 * - a joint distribution for the experimental observables,
 * - a joint distribution for the nuisance parameters.
 *
 * The optimization vector is ordered as \f$\theta = (p, \eta)\f$, where
 * \f$p\f$ denotes fitted parameters and \f$\eta\f$ denotes nuisance parameters.
 *
 * @see IProfileableLikelihood
 * @see JointDistribution
 */

/**
 * @brief Model function signature used by @ref BaseLikelihood.
 *
 * The function receives the fitted-parameter block \f$p\f$ and the nuisance
 * block \f$\eta\f$, and returns the corresponding model prediction in
 * observable space.
 */
using ModelFn = std::function<std::vector<double>(const std::vector<double>& p, const std::vector<double>& eta)>;

/**
 * @struct LikelihoodContext
 * @brief Shared immutable-like data required to evaluate a likelihood.
 *
 * The context owns the probability distributions used by the likelihood and
 * stores the nominal observable and parameter definitions. It is shared through
 * @c std::shared_ptr by @ref BaseLikelihood so that several likelihood-related
 * objects can reuse the same statistical description.
 */
struct LikelihoodContext {
    std::unique_ptr<JointDistribution> nuisance_dist;       ///< Joint distribution of nuisance parameters.
    std::unique_ptr<JointDistribution> exp_obs_dist;        ///< Joint distribution of observable residuals.
    std::vector<double>exp_obs_values;                      ///< Nominal experimental observable values.
    std::vector<fit_app::ParameterDefinition> nuis_defs;    ///< Definitions of nuisance parameters.
    std::vector<fit_app::ParameterDefinition> fp_defs;      ///< Definitions of fitted parameters.
};

/**
 * @class BaseLikelihood
 * @brief Default implementation of a profileable negative log-likelihood.
 *
 * @ref BaseLikelihood evaluates
 * \f[
 *   -\log L(p, \eta)
 *   = -\left[\log f_{obs}(m(p,\eta) - y_{obs})
 *       + \log f_{nuis}(\eta)\right],
 * \f]
 * where \f$m(p,\eta)\f$ is the model prediction and \f$y_{obs}\f$ are the
 * nominal experimental observable values stored in @ref LikelihoodContext.
 *
 * The class also exposes diagnostic helpers for profiling workflows, including
 * central values, split-parameter evaluation, residual computation and local
 * curvature matrices.
 *
 * @note If model evaluation or density evaluation produces a non-finite value,
 *       the implementation returns a large finite penalty value instead of
 *       propagating NaNs through the optimizer.
 */
class BaseLikelihood : public IProfileableLikelihood {
public:
    /**
     * @brief Constructs a likelihood from a model and statistical context.
     *
     * @param model Model function used to compute observable predictions.
     * @param ctx   Shared likelihood context containing distributions and
     *              parameter definitions.
     * @param p_dim Dimension of the fitted-parameter block \f$p\f$.
     */
    BaseLikelihood(const ModelFn& model, std::shared_ptr<LikelihoodContext> ctx, size_t p_dim);

    /**
     * @brief Evaluates the negative log-likelihood for a full parameter vector.
     *
     * The input vector is split as \f$\theta = (p, \eta)\f$ using the fitted
     * parameter dimension provided at construction.
     *
     * @param theta Full parameter vector ordered as fitted parameters followed
     *              by nuisance parameters.
     * @return Negative log-likelihood value, or a large finite penalty if the
     *         evaluation is numerically invalid.
     */
    double nll(const std::vector<double>& theta) const override;

    /**
     * @brief Returns all parameter definitions in likelihood-vector order.
     *
     * The returned list contains fitted-parameter definitions first, followed by
     * nuisance-parameter definitions.
     *
     * @return Ordered list of fitted and nuisance parameter definitions.
     */
    std::vector<fit_app::ParameterDefinition> get_param_defs() const override {
        auto p_defs = ctx->fp_defs;
        p_defs.insert(p_defs.end(), ctx->nuis_defs.begin(), ctx->nuis_defs.end());
        return p_defs;
    }

    /**
     * @brief Returns the total likelihood dimension.
     *
     * @return Sum of fitted-parameter and nuisance-parameter dimensions.
     */
    std::size_t dim() const override;

    /**
     * @brief Enables debug tracing for likelihood evaluations.
     *
     * When enabled, the first @p max_evals evaluations print diagnostic
     * information such as the negative log-likelihood value, observable and
     * nuisance log-density terms, and maximum shifts relative to the first
     * traced point.
     *
     * @param max_evals Maximum number of evaluations to print.
     */
    void enable_debug_trace(std::size_t max_evals = 25) {
        debug_trace_enabled_ = true;
        debug_trace_max_evals_ = max_evals;
        debug_eval_count_ = 0;
        debug_have_ref_theta_ = false;
        debug_have_ref_res_ = false;
        debug_ref_theta_.clear();
        debug_ref_res_.clear();
    }

    /**
     * @brief Disables debug tracing of likelihood evaluations.
     */
    void disable_debug_trace() {
        debug_trace_enabled_ = false;
    }

    /** @copydoc IProfileableLikelihood::p_dimension() */
    std::size_t p_dimension() const override;

    /** @copydoc IProfileableLikelihood::eta_dimension() */
    std::size_t eta_dimension() const override;

    /** @copydoc IProfileableLikelihood::central_p() */
    std::vector<double> central_p() const override;

    /** @copydoc IProfileableLikelihood::central_eta() */
    std::vector<double> central_eta() const override;

    /** @copydoc IProfileableLikelihood::predict(const std::vector<double>&, const std::vector<double>&) */
    std::vector<double> predict(
        const std::vector<double>& p,
        const std::vector<double>& eta
    ) const override;

    /** @copydoc IProfileableLikelihood::residuals(const std::vector<double>&, const std::vector<double>&) */
    std::vector<double> residuals(
        const std::vector<double>& p,
        const std::vector<double>& eta
    ) const override;

    /** @copydoc IProfileableLikelihood::nll_from_split(const std::vector<double>&, const std::vector<double>&) */
    double nll_from_split(
        const std::vector<double>& p,
        const std::vector<double>& eta
    ) const override;

    /** @copydoc IProfileableLikelihood::observable_curvature(const std::vector<double>&) */
    RealMatrix observable_curvature(
        const std::vector<double>& r
    ) const override;

    /** @copydoc IProfileableLikelihood::nuisance_curvature(const std::vector<double>&) */
    RealMatrix nuisance_curvature(
        const std::vector<double>& eta
    ) const override;

protected:
    std::shared_ptr<LikelihoodContext> ctx;             ///< Shared statistical context used by the likelihood.
    ModelFn model;                                      ///< Model function evaluated by the likelihood.
    std::size_t p_dim;                                  ///< Dimension of the fitted-parameter block.

    mutable bool debug_trace_enabled_ = false;          ///< Whether evaluation debug tracing is enabled.
    mutable std::size_t debug_trace_max_evals_ = 0;     ///< Maximum number of traced evaluations.
    mutable std::size_t debug_eval_count_ = 0;          ///< Number of evaluations already traced.

    mutable bool debug_have_ref_theta_ = false;         ///< Whether the debug reference parameter vector is set.
    mutable bool debug_have_ref_res_ = false;           ///< Whether the debug reference residual vector is set.
    mutable std::vector<double> debug_ref_theta_;       ///< First traced parameter vector used as debug reference.
    mutable std::vector<double> debug_ref_res_;         ///< First traced residual vector used as debug reference.
    
private:
    /**
     * @brief Computes the maximum absolute component-wise difference.
     *
     * @param a First vector.
     * @param b Second vector.
     * @return Maximum absolute difference, or NaN when vector sizes differ.
     */
    static double max_abs_diff(const std::vector<double>& a,
                               const std::vector<double>& b) {
        if (a.size() != b.size()) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        double out = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i) {
            out = std::max(out, std::abs(a[i] - b[i]));
        }
        return out;
    }

    /**
     * @brief Prints one debug trace line for a likelihood evaluation when enabled.
     *
     * @param theta     Full parameter vector.
     * @param p         Fitted-parameter block.
     * @param eta       Nuisance-parameter block.
     * @param res       Observable residual vector.
     * @param ell_obs   Observable log-density contribution.
     * @param ell_nuis  Nuisance log-density contribution.
     * @param nll_value Final negative log-likelihood value.
     */
    void maybe_log_debug_eval(const std::vector<double>& theta,
                              const std::vector<double>& p,
                              const std::vector<double>& eta,
                              const std::vector<double>& res,
                              double ell_obs,
                              double ell_nuis,
                              double nll_value) const;
};

#endif // BASELIKELIHOOD_H
