#ifndef WITHPROFILING_H
#define WITHPROFILING_H

#include "ILikelihood.h"
#include "Math.h"
#include "Indexing.h"
#include "GradientHelper.h"

/**
 * @file Profiler.h
 * @brief Generic likelihood profiling engine.
 *
 * This file declares the request/result structures and the @ref Profiler class
 * used to minimize likelihoods under fixed-parameter constraints.
 *
 * @see ILikelihood
 * @see IProfileableLikelihood
 * @see GradientHelper.h
 */

/**
 * @enum ProfilerMode
 * @brief Available algorithms for profiling free parameters.
 */
enum class ProfilerMode {
    MINUIT,             ///< Use the configured numerical minimizer backend.
    LAPLACE_NUISANCE    ///< Use the Laplace nuisance approximation when possible, with Minuit fallback.
};

/**
 * @struct ProfileRequest
 * @brief Input specification for a constrained profile minimization.
 *
 * The request partitions the full parameter vector into fixed and free
 * coordinates. Fixed parameters are assigned explicit values; free parameters
 * are minimized by the profiler using @ref start as the initial point.
 */
struct ProfileRequest {
    std::vector<std::size_t> free_params;           ///< Global indices of parameters released during profiling.
    std::map<std::size_t, double> fixed_params;     ///< Fixed parameter values keyed by global parameter index.
    std::vector<double> start;                      ///< Initial full parameter vector used as a warm start.
};

/**
 * @struct ProfileResult
 * @brief Output of a constrained profile minimization.
 */
struct ProfileResult {
    double nll_hat = 1e300;                     ///< Minimum or approximate profiled NLL value.
    std::map<std::size_t, double> theta_hat;    ///< Profiled parameter values keyed by global parameter index.
    bool converged = false;                     ///< True when the profiler accepted the result as reliable.
};

/**
 * @class Profiler
 * @brief Minimizes a likelihood over free parameters subject to fixed constraints.
 *
 * The profiler supports two modes:
 * - @ref ProfilerMode::MINUIT, which delegates to the configured fit backend,
 * - @ref ProfilerMode::LAPLACE_NUISANCE, which attempts a fast Laplace
 *   approximation for nuisance-only profiling and falls back to Minuit when
 *   the request is incompatible or the approximation fails.
 */
class Profiler {
public:
    /**
     * @brief Constructs a profiler.
     *
     * @param minimizer Fit backend used for numerical minimization and fallbacks.
     * @param mode Profiling algorithm selected for future calls.
     */
    Profiler(
        std::shared_ptr<fit_app::IFitBackend> minimizer,
        ProfilerMode mode = ProfilerMode::MINUIT
    );

    /**
     * @brief Profiles a likelihood according to a request.
     *
     * If no free parameters are provided, the method evaluates the likelihood
     * directly at the fixed point. Otherwise, it dispatches to the configured
     * profiling mode.
     *
     * @param base Likelihood to profile.
     * @param pr Fixed/free parameter specification and warm-start vector.
     *
     * @return Profile result containing the profiled NLL, argmin and convergence status.
     *
     * @throws std::runtime_error for inconsistent request dimensions or invalid
     *         fixed parameter indices.
     */
    ProfileResult profile(
        std::shared_ptr<ILikelihood> base,
        const ProfileRequest& pr
    ) const;

private:
    /**
     * @brief Profiles with the configured Minuit-like backend.
     *
     * @param base Likelihood to profile.
     * @param pr Fixed/free parameter specification and warm-start vector.
     *
     * @return Profile result accepted according to backend diagnostics.
     */
    ProfileResult profile_minuit(
        std::shared_ptr<ILikelihood> base,
        const ProfileRequest& pr
    ) const;

    /**
     * @brief Profiles nuisance parameters with the Laplace approximation.
     *
     * This fast path requires a likelihood implementing @ref IProfileableLikelihood
     * and supports only requests where fit parameters are fixed and nuisance
     * parameters are free.
     *
     * @param base Profileable likelihood to profile.
     * @param pr Fixed/free parameter specification.
     *
     * @return Approximate profiled result.
     *
     * @throws std::runtime_error if the request is incompatible with the Laplace
     *         nuisance profiler.
     */
    ProfileResult profile_laplace_nuisance(
        std::shared_ptr<ILikelihood> base,
        const ProfileRequest& pr
    ) const;

    std::shared_ptr<fit_app::IFitBackend> minimizer;    ///< Backend minimizer used for full numerical profiling.
    ProfilerMode mode = ProfilerMode::MINUIT;           ///< Selected profiling mode.
};

#endif // WITHPROFILING_H
