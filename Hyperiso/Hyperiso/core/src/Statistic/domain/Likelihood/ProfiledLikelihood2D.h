#ifndef PROFILEDLIKELIHOOD2D_H
#define PROFILEDLIKELIHOOD2D_H

#include "ILikelihood.h"
#include "Profiler.h"
#include "IProfilingStrategy.h"

/**
 * @file ProfiledLikelihood2D.h
 * @brief Two-dimensional profiled likelihood wrapper.
 *
 * This file declares @ref ProfiledLikelihood2D, a small adapter that evaluates
 * a base likelihood on a two-dimensional scan grid while delegating the actual
 * profiling of free parameters to @ref Profiler.
 *
 * @see ILikelihood
 * @see Profiler
 * @see IProfilingStrategy
 */

/**
 * @class ProfiledLikelihood2D
 * @brief Evaluates a profiled NLL as a function of two scan coordinates.
 *
 * The class combines:
 * - a base likelihood,
 * - a profiler backend,
 * - a profiling strategy defining which parameters are fixed or free.
 *
 * The last successful profiled point is retained as a warm start for the next
 * evaluation, which improves stability and performance in grid scans.
 */
class ProfiledLikelihood2D {
public:
    /**
     * @brief Default constructor.
     *
     * Creates an empty wrapper. The object must be assigned or reconstructed
     * with valid dependencies before @ref profiled_nll is called.
     */
    ProfiledLikelihood2D() = default;

    /**
     * @brief Constructs a two-dimensional profiled likelihood wrapper.
     *
     * @param base Base likelihood to evaluate.
     * @param profiler Profiler used to minimize over free dimensions.
     * @param profiling_strategy Strategy used to build profile requests for
     *                           each 2D scan point.
     */
    ProfiledLikelihood2D(
        std::shared_ptr<ILikelihood> base, 
        std::shared_ptr<Profiler> profiler, 
        std::shared_ptr<IProfilingStrategy> profiling_strategy
    );
    
    /**
     * @brief Evaluates the profiled NLL at one 2D scan point.
     *
     * @param px Value of the first scanned coordinate.
     * @param py Value of the second scanned coordinate.
     *
     * @return Profiled NLL value. A large penalty is returned if profiling
     *         fails without producing a finite NLL.
     */
    double profiled_nll(double px, double py);


    /**
     * @brief Returns the parameter definitions for the two scan coordinates.
     *
     * @return Array containing the definitions of the first and second scanned
     *         parameters, in that order.
     */
    std::array<fit_app::ParameterDefinition, 2> get_param_defs() const;

private:
    std::shared_ptr<ILikelihood> base;                          ///< Base likelihood being profiled.
    std::shared_ptr<Profiler> profiler;                         ///< Profiler used to minimize over free parameters.
    std::shared_ptr<IProfilingStrategy> profiling_strategy;     ///< Strategy defining fixed/free parameters for each scan point.
    std::map<std::size_t, double> last;                         ///< Warm-start map storing the last profiled argmin.
};

#endif // PROFILEDLIKELIHOOD2D_H
