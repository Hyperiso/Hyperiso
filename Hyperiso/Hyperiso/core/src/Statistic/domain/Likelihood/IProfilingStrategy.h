#ifndef IPROFILINGSTRATEGY_H
#define IPROFILINGSTRATEGY_H

#include "Profiler.h"

/**
 * @file IProfilingStrategy.h
 * @brief Profiling strategies used to build two-dimensional likelihood scan requests.
 *
 * This file defines:
 * - @ref FitResult, the global best-fit point used as a profiling reference,
 * - @ref IProfilingStrategy, the strategy interface,
 * - @ref SliceProfilingStrategy, which fixes all parameters of interest except
 *   nuisance parameters,
 * - @ref ProjectionProfilingStrategy, which fixes only the two scanned
 *   coordinates and profiles all remaining dimensions.
 *
 * @see Profiler
 * @see ProfileRequest
 * @see ProfiledLikelihood2D
 */

/**
 * @struct FitResult
 * @brief Summary of a global likelihood fit.
 *
 * The result stores the best-fit parameters of interest, profiled nuisances and
 * uncertainty information used to seed subsequent profile scans.
 */
struct FitResult {
    std::vector<double> p_hat;      ///< Maximum-likelihood estimates for parameters of interest.
    std::vector<double> eta_hat;    ///< Nuisance estimates at the global maximum-likelihood point.
    std::vector<double> p_hat_std;  ///< Standard deviations of the parameters of interest.
    RealMatrix p_hat_correlations;  ///< Correlation matrix for the parameters of interest.
    double ell_hat {0.0};           ///< Minimum NLL value at the global best-fit point.
};

/**
 * @class IProfilingStrategy
 * @brief Interface for converting 2D scan coordinates into profiler requests.
 *
 * A profiling strategy owns the indices of the two coordinates to scan and a
 * reference fit result. Concrete strategies decide which remaining parameters
 * are fixed or released when building a @ref ProfileRequest.
 */
class IProfilingStrategy {
public:
    /**
     * @brief Constructs a profiling strategy.
     *
     * @param x_id Global parameter index used as the first scan coordinate.
     * @param y_id Global parameter index used as the second scan coordinate.
     * @param fr Reference global fit result used for central values and warm starts.
     */
    IProfilingStrategy(std::size_t x_id, std::size_t y_id, const FitResult& fr);

    /**
     * @brief Virtual destructor.
     */
    virtual ~IProfilingStrategy() = default;

    /**
     * @brief Builds a profiler request for one 2D scan point.
     *
     * @param px Value assigned to the first scanned coordinate.
     * @param py Value assigned to the second scanned coordinate.
     * @param current_argmin Warm-start map from free parameter indices to their
     *                       last profiled values.
     *
     * @return Fully specified profiling request with fixed parameters, free
     *         parameters and starting values.
     *
     * @throws std::runtime_error if required warm-start values are missing.
     */
    virtual ProfileRequest build_request(
        double px,
        double py,
        const std::map<std::size_t, double>& current_argmin
    ) const = 0;

    /**
     * @brief Creates the initial warm-start map for the strategy.
     *
     * @return Map from free parameter indices to their initial values.
     */
    virtual std::map<std::size_t, double> init_warm_start() const = 0;

    /**
     * @brief Returns the global index of the first scanned coordinate.
     *
     * @return First scan-coordinate index.
     */
    std::size_t get_x_id() { return x_id; };

    /**
     * @brief Returns the global index of the second scanned coordinate.
     *
     * @return Second scan-coordinate index.
     */
    std::size_t get_y_id() { return y_id; };

protected:
    std::size_t x_id, y_id; ///< Global indices of the two profiled scan coordinates.
    FitResult fr;           ///< Reference global fit used for central values and warm starts.
};

/**
 * @class SliceProfilingStrategy
 * @brief Builds slice-profile requests with fixed parameters of interest.
 *
 * In this strategy, all parameters of interest are fixed to their global
 * best-fit values, except the two scanned coordinates which are fixed to the
 * requested scan values. Only nuisance parameters are profiled.
 */
class SliceProfilingStrategy : public IProfilingStrategy {
public:
    /**
     * @brief Constructs a slice profiling strategy.
     *
     * @param x_id Global index of the first scan coordinate.
     * @param y_id Global index of the second scan coordinate.
     * @param fr Reference global fit result.
     */
    SliceProfilingStrategy(std::size_t x_id, std::size_t y_id, const FitResult& fr) : IProfilingStrategy(x_id, y_id, fr) {}

    /**
     * @copydoc IProfilingStrategy::build_request
     */
    ProfileRequest build_request(
        double px,
        double py,
        const std::map<std::size_t, double>& current_argmin
    ) const override;

    /**
     * @copydoc IProfilingStrategy::init_warm_start
     */
    std::map<std::size_t, double> init_warm_start() const override;
};

/**
 * @class ProjectionProfilingStrategy
 * @brief Builds projection-profile requests with all non-scanned dimensions free.
 *
 * In this strategy, only the two scanned coordinates are fixed. All other
 * parameters, including the remaining parameters of interest and all nuisance
 * parameters, are released and warm-started from the previous argmin.
 */
class ProjectionProfilingStrategy : public IProfilingStrategy {
public:
    /**
     * @brief Constructs a projection profiling strategy.
     *
     * @param x_id Global index of the first scan coordinate.
     * @param y_id Global index of the second scan coordinate.
     * @param fr Reference global fit result.
     */
    ProjectionProfilingStrategy(std::size_t x_id, std::size_t y_id, const FitResult& fr) : IProfilingStrategy(x_id, y_id, fr) {}

    /**
     * @copydoc IProfilingStrategy::build_request
     */
    ProfileRequest build_request(
        double px,
        double py,
        const std::map<std::size_t, double>& current_argmin
    ) const override;

    /**
     * @copydoc IProfilingStrategy::init_warm_start
     */
    std::map<std::size_t, double> init_warm_start() const override;
};

#endif // IPROFILINGSTRATEGY_H
