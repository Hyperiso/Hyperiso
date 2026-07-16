#ifndef CONTOUR_OBSERVER_H
#define CONTOUR_OBSERVER_H

#include <cstddef>
#include <functional>
#include <string>

/**
 * @file ContourObserver.h
 * @brief Progress-event types and callback interface for contour extraction.
 *
 * The contour machinery can report coarse-grained lifecycle events and
 * fine-grained path-point events through @ref ContourProgressCallback. This is
 * intentionally UI-agnostic: consumers may log events, update a progress bar,
 * or stream contour points as they are produced.
 */

/**
 * @enum ContourProgressEventType
 * @brief Type of progress notification emitted during contour computation.
 */
enum class ContourProgressEventType {
    Started,        ///< Contour computation has started.
    PathPoint,      ///< A new point on a contour path is available.
    PathFinished,   ///< A contour path has been completed.
    Finished,       ///< Contour computation finished normally.
    Failed          ///< Contour computation failed and no valid result was produced.
};

/**
 * @struct ContourProgressEvent
 * @brief Payload passed to contour progress callbacks.
 *
 * Not every field is meaningful for every event type. For example, @ref x and
 * @ref y are mainly relevant for @ref ContourProgressEventType::PathPoint,
 * whereas @ref n_paths, @ref n_points, and @ref elapsed_seconds are typically
 * filled for final lifecycle events.
 */
struct ContourProgressEvent {
    ContourProgressEventType type = ContourProgressEventType::Started;  ///< Event category.

    double level = 0.0;                                                 ///< Contour level associated with the event.

    std::size_t path_id = 0;                                            ///< Index of the contour path being reported.
    std::size_t point_id = 0;                                           ///< Index of the point within the current path.

    double x = 0.0;                                                     ///< X coordinate of a reported contour point.
    double y = 0.0;                                                     ///< Y coordinate of a reported contour point.

    std::size_t n_paths = 0;                                            ///< Number of paths produced so far or in total.
    std::size_t n_points = 0;                                           ///< Number of points produced so far or in total.

    double elapsed_seconds = 0.0;                                       ///< Wall-clock time elapsed since computation start.

    std::string message;                                                ///< Optional human-readable status message.
};

/**
 * @brief Callback signature used to receive contour progress events.
 */
using ContourProgressCallback = std::function<void(const ContourProgressEvent&)>;

#endif  // CONTOUR_OBSERVER_H