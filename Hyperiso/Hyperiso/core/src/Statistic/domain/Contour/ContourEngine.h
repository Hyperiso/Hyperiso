#ifndef CONTOURENGINE_H
#define CONTOURENGINE_H

#include <chrono>

#include "ILikelihood.h"
#include "IProfilingStrategy.h"
#include "IContourExtractor.h"
#include "ProfiledLikelihood2D.h"
#include "Profiler.h"
#include "JointDistribution.h"
#include "WithGaussianConstraints.h"
#include "WithFallback.h"
#include "AMSContourExtractor.h"
#include "MnContourExtractor.h"
#include "GaussianMarginal.h"
#include "GaussianCopula.h"
#include "ContourObserver.h"
#include "Profiler.h"

/**
 * @file ContourEngine.h
 * @brief High-level engine for profiled two-dimensional likelihood contours.
 *
 * The contour engine combines three concerns:
 * - construction of a two-dimensional profiled likelihood,
 * - optional prior-like Gaussian constraints for projected parameters,
 * - extraction of the requested contour with a configurable algorithm.
 *
 * @see ProfiledLikelihood2D
 * @see Profiler
 * @see IContourExtractor
 */

/**
 * @enum ProfilingMethod
 * @brief Strategy used to reduce the full likelihood to a 2D profile.
 */
enum class ProfilingMethod {
    SLICE,                          ///< Fix all fit parameters except nuisances; scan only the selected axes.
    FREE_PROJECTION,                ///< Fix selected axes and profile all remaining parameters.
    PRIOR_CONSTRAINED_PROJECTION    ///< Profile remaining parameters with Gaussian constraints around the fit result.
};


/**
 * @enum ContourAlgorithm
 * @brief Algorithm family used to extract contour paths from the profiled field.
 */
enum class ContourAlgorithm {
    AMS,    ///< Adaptive/marching-squares contour extractor.
    MINUIT  ///< Minuit contour extractor.
};

/**
 * @struct ContourConfig
 * @brief Configuration object for @ref ContourEngine.
 */
struct ContourConfig {
    std::size_t x_id, y_id;                                     ///< Indices of the two fit parameters displayed on the contour axes.
    FitResult fr;                                               ///< Global fit result used for central values, uncertainties, and correlations.
    ProfilingMethod profiling_method;                           ///< Profiling strategy used for non-displayed parameters.
    ContourAlgorithm primary_contour_method;                    ///< Primary contour extraction algorithm.
    std::optional<ContourAlgorithm> fallback_contour_method;    ///< Optional fallback algorithm if the primary fails.
    ProfilerMode profile_backend;                               ///< Backend used to profile nuisance/free parameters.

    ContourProgressCallback on_progress {};                     ///< Optional progress callback invoked during contour computation.
};

/**
 * @class ContourEngine
 * @brief Orchestrates profiled likelihood contour computation.
 *
 * A @ref ContourEngine instance owns the profiled 2D likelihood view and the
 * contour extractor selected in @ref ContourConfig. Calling @ref compute_contour
 * evaluates the profiled likelihood field, converts the requested Gaussian-like
 * significance into the corresponding 2D likelihood-ratio level, and extracts
 * contour paths within the requested bounds.
 */
class ContourEngine {
public:
    /**
     * @brief Constructs a contour engine for a base likelihood.
     *
     * Depending on @ref ContourConfig::profiling_method, the base likelihood may
     * be wrapped in @ref WithGaussianConstraints before building the profiled 2D
     * likelihood. If a fallback algorithm is configured, the primary extractor is
     * wrapped in @ref WithFallback.
     *
     * @param base Full-dimensional likelihood to profile.
     * @param cfg Contour computation configuration.
     */
    ContourEngine(std::shared_ptr<ILikelihood> base, const ContourConfig& cfg);

    /**
     * @brief Computes a profiled 2D contour for a requested significance.
     *
     * @param z One-dimensional Gaussian-equivalent significance. The engine
     *          converts it to the appropriate two-dimensional chi-square level.
     * @param bounds Extraction domain as {xmin, xmax, ymin, ymax}.
     * @param resolution Requested extraction resolution.
     * @return Extracted contour at the converted likelihood-ratio level.
     *
     * @throws std::exception if profiling or contour extraction fails without a
     *         configured fallback being able to recover.
     */
    Contour compute_contour(double z, std::array<double, 4> bounds, std::size_t resolution);

private:
    /**
     * @brief Builds the Gaussian constraint distribution for hidden fit parameters.
     *
     * The distribution is derived from the fitted central values, standard
     * deviations, and correlation matrix after removing the displayed axes.
     *
     * @return Joint Gaussian distribution over the constrained projected axes.
     */
    std::shared_ptr<JointDistribution> build_constraints_distribution();

    /**
     * @brief Creates a concrete contour extractor from an algorithm selector.
     *
     * @param ca Requested contour algorithm.
     * @return Polymorphic contour extractor instance.
     *
     * @throws std::invalid_argument if \p ca is not recognized.
     */
    std::shared_ptr<IContourExtractor> build_contour_extractor(ContourAlgorithm ca);

    ContourConfig cfg;                              ///< Engine configuration.
    ProfiledLikelihood2D likelihood;                ///< Two-dimensional profiled likelihood view.
    std::shared_ptr<IContourExtractor> extractor;   ///< Selected contour extractor, possibly with fallback.
};

#endif // CONTOURENGINE_H
