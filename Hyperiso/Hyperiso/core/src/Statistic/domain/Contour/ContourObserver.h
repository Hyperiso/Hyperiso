#ifndef __CONTOUR_OBSERVER_H__
#define __CONTOUR_OBSERVER_H__

#include <cstddef>
#include <functional>
#include <string>

enum class ContourProgressEventType {
    Started,
    PathPoint,
    PathFinished,
    Finished,
    Failed
};

struct ContourProgressEvent {
    ContourProgressEventType type = ContourProgressEventType::Started;

    double level = 0.0;

    std::size_t path_id = 0;
    std::size_t point_id = 0;

    double x = 0.0;
    double y = 0.0;

    std::size_t n_paths = 0;
    std::size_t n_points = 0;

    double elapsed_seconds = 0.0;

    std::string message;
};

using ContourProgressCallback = std::function<void(const ContourProgressEvent&)>;

#endif