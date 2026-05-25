#ifndef FIT_H
#define FIT_H

#include <vector>
#include <functional>

#include "BaseLikelihood.h"
#include "IProfilingStrategy.h"
#include "IContourExtractor.h"
#include "Math.h"
#include "ContourEngine.h"
#include "ContourObserver.h"

/**
 * @file Fit.h
 * @brief High-level maximum-likelihood fitting and confidence-contour API.
 *
 * This header exposes the user-facing fitting workflow built around
 * @ref MLFitter.  It combines a likelihood model, a Minuit-compatible fitting
 * backend, optional covariance fallbacks, and 2D contour extraction utilities.
 *
 * @see BaseLikelihood
 * @see FitResult
 * @see ContourEngine
 */

/**
 * @struct MLFitOptions
 * @brief Runtime options controlling the global maximum-likelihood fit.
 *
 * These options are forwarded to the underlying fit backend whenever the
 * corresponding value is enabled or strictly positive.  Values left at their
 * neutral defaults delegate the choice to the backend.
 */
struct MLFitOptions {
    /** Whether to request a HESSE covariance evaluation after minimization. */
    bool run_hesse = true;

    /** Whether to request MINOS uncertainties when supported by the backend. */
    bool request_minos = false;

    /** Enables verbose backend output. */
    bool verbose = false;

    /** Minuit strategy. A value of 0 keeps the backend default. */
    unsigned strategy = 0;

    /** Maximum number of function calls. A value of 0 keeps the backend default. */
    unsigned max_fcn = 0;

    /** Fit tolerance. A non-positive value keeps the backend default. */
    double tolerance = 0.0;

    /**
     * Enables a numerical profile-Hessian fallback when no usable covariance is
     * returned by the main fit.
     */
    bool allow_profile_hessian_fallback = true;

    /** Multiplicative scale applied to finite-difference steps in the fallback. */
    double profile_hessian_step_scale = 1.0;

    /** Relative eigenvalue floor used to regularize the fallback Hessian. */
    double profile_hessian_eig_floor_rel = 1e-8;

    /** Enables logging of the first likelihood evaluations for debugging. */
    bool trace_first_evals = false;

    /** Maximum number of initial likelihood evaluations printed when tracing. */
    std::size_t trace_max_evals = 25;
};

/**
 * @enum ProfileBackend
 * @brief Backend used to profile nuisance parameters during contour building.
 */
enum class ProfileBackend {
    /** Full numerical profiling with the Minuit backend. */
    MINUIT,

    /** Fast Laplace approximation for nuisance profiling, with Minuit fallback. */
    LAPLACE_NUISANCE
};

/**
 * @struct ContourOptions
 * @brief Runtime options controlling 2D contour computation.
 *
 * The contour is computed in the plane of two selected fit parameters.  The
 * remaining degrees of freedom are handled according to @ref profiling_method
 * and @ref profile_backend.
 */
struct ContourOptions {
    /** Strategy used to define fixed and profiled parameters on the contour. */
    ProfilingMethod profiling_method = ProfilingMethod::SLICE;

    /** Backend used to profile nuisance parameters. */
    ProfileBackend profile_backend = ProfileBackend::LAPLACE_NUISANCE;

    /** Primary algorithm used to extract the contour. */
    ContourAlgorithm primary_contour_method = ContourAlgorithm::MINUIT;

    /** Optional fallback contour extractor used if the primary extractor fails. */
    std::optional<ContourAlgorithm> fallback_contour_method;

    /** Resolution parameter forwarded to the selected contour extractor. */
    std::size_t resolution = 40;

    /** Optional callback receiving progress events during contour extraction. */
    ContourProgressCallback on_progress {};
};

/**
 * @class MLFitter
 * @brief High-level driver for likelihood minimization and 2D contour extraction.
 *
 * @ref MLFitter owns a profileable likelihood, performs the global maximum-
 * likelihood fit, stores the best-fit state, and subsequently builds confidence
 * contours for selected fit-parameter pairs.
 *
 * Typical usage:
 * - construct the fitter from a @ref LikelihoodContext and a model function,
 * - call @ref maximum_likelihood_fit,
 * - call @ref contour for parameter pairs of interest.
 *
 * @note Contour computation requires a successful prior call to
 *       @ref maximum_likelihood_fit.
 */
class MLFitter {
public:
    /**
     * @brief Constructs a fitter from a likelihood context and model function.
     *
     * A @ref BaseLikelihood is created internally.  The number of fit
     * parameters is inferred from @c ctx->fp_defs.
     *
     * @param ctx Likelihood context containing distributions, observations and
     *        parameter definitions.
     * @param model Model function mapping fit and nuisance parameters to
     *        predicted observables.
     * @param options Fit options controlling the minimization and covariance
     *        extraction.
     */
    MLFitter(std::shared_ptr<LikelihoodContext> ctx, const ModelFn& model, MLFitOptions options = {});

    /**
     * @brief Constructs a fitter from an already configured base likelihood.
     *
     * @param like Likelihood object used by the fitter.
     * @param options Fit options controlling the minimization and covariance
     *        extraction.
     *
     * @throws std::invalid_argument if @p like is null.
     */
    explicit MLFitter(std::shared_ptr<BaseLikelihood> like, MLFitOptions options = {});
    
    /**
     * @brief Runs the global maximum-likelihood fit.
     *
     * The fit parameters are initialized from @p p0, while nuisance parameters
     * keep the values stored in the likelihood parameter definitions.  The
     * returned @ref FitResult contains the best-fit values, profiled nuisance
     * estimates, minimum NLL, and covariance-derived uncertainties when
     * available.
     *
     * @param p0 Initial values for the fit-parameter block.
     * @return Global fit result.
     */
    FitResult maximum_likelihood_fit(const std::vector<double>& p0);

    /**
     * @brief Computes a 2D confidence contour after the global fit.
     *
     * @param x_id Index of the first fit parameter in the global fit-parameter
     *        block.
     * @param y_id Index of the second fit parameter in the global fit-parameter
     *        block.
     * @param z Gaussian-equivalent confidence level in standard deviations.
     * @param bounds Contour bounds ordered as @c {xmin, xmax, ymin, ymax}.
     * @param options Contour and profiling options.
     * @return Extracted contour at the requested level.
     *
     * @pre @ref maximum_likelihood_fit must have been called successfully.
     */
    Contour contour(std::size_t x_id, std::size_t y_id, double z, std::array<double, 4> bounds, ContourOptions options) const;

private:
    std::shared_ptr<BaseLikelihood> like_;  ///< Likelihood minimized by this fitter.
    MLFitOptions fit_options_;              ///< Stored options used for the global fit.
    FitResult master_fit_result;            ///< Cached result of the latest global fit.
    bool master_fit_success = false;        ///< Whether the latest global fit produced valid parameters.
};

#endif // FIT_H
