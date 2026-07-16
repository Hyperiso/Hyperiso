#ifndef AMSCONTOUREXTRACTOR_H
#define AMSCONTOUREXTRACTOR_H

#include <utility>

#include "IContourExtractor.h"
#include "ContourObserver.h"
#include "Math.h"

/**
 * @file AMSContourExtractor.h
 * @brief Adaptive marching-squares contour extractor.
 *
 * This header declares @ref AMSContourExtractor, an implementation of
 * @ref IContourExtractor that samples a two-dimensional scalar field with a
 * marching-squares based algorithm and optionally reports progress events for
 * each extracted path and point.
 *
 * @see IContourExtractor
 * @see ContourObserver.h
 * @see WithFallback
 */

/**
 * @class AMSContourExtractor
 * @brief Extracts iso-contours using an adaptive marching-squares backend.
 *
 * The extractor converts a requested grid resolution into a marching-squares
 * recursion depth, evaluates the provided scalar field, and returns every path
 * found at the requested contour level.
 *
 * Progress reporting is optional. When a callback is provided, the extractor
 * emits @ref ContourProgressEventType::PathPoint events for every contour point
 * and @ref ContourProgressEventType::PathFinished events after each path.
 *
 * @note The implementation caps the internal recursion depth to keep extraction
 *       cost bounded and reports failure by returning a contour with
 *       `success == false` rather than throwing for recoverable extraction
 *       errors.
 */
class AMSContourExtractor: public IContourExtractor {
public:
    /**
     * @brief Constructs an AMS contour extractor.
     *
     * @param on_progress Optional callback invoked while contour paths are
     *        emitted. Pass an empty callback to disable progress reporting.
     */
    explicit AMSContourExtractor(ContourProgressCallback on_progress = {})
        : on_progress_(std::move(on_progress)) {}

        /**
     * @copydoc IContourExtractor::extract
     *
     * @return A contour whose @ref Contour::success flag is true only when at
     *         least one path was found.
     */
    Contour extract(const ScalarField2D& field, const ContourRequest& cr) override;
private:
    ContourProgressCallback on_progress_;   ///< Optional progress callback used during path emission.
};

#endif // AMSCONTOUREXTRACTOR_H
