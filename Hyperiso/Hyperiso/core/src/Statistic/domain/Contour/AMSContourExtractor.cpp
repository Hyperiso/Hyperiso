#include "AMSContourExtractor.h"

#include <cmath>
#include <exception>
#include <iostream>

Contour AMSContourExtractor::extract(const ScalarField2D& field, const ContourRequest& cr) {
    Contour c;
    c.level = cr.level;
    c.success = false;

    try {
        const std::size_t safe_resolution = std::max<std::size_t>(2, cr.resolution);
        std::size_t depth = static_cast<std::size_t>(
            std::ceil(std::log2(static_cast<double>(safe_resolution)))
        );

        depth = std::min<std::size_t>(depth, 8);

        MarchingSquaresExtractor mse(field, cr.bounds, depth);
        c.paths = mse.find_iso_contour(cr.level);

        // Empty paths mean AMS did not find a contour. Let WithFallback try MINUIT.
        c.success = !c.paths.empty();
    } catch (const std::exception& e) {
        std::cout << "[AMS] extraction failed: " << e.what() << std::endl;
        c.paths.clear();
        c.success = false;
        return c;
    } catch (...) {
        std::cout << "[AMS] extraction failed with unknown exception" << std::endl;
        c.paths.clear();
        c.success = false;
        return c;
    }

    if (on_progress_) {
        std::size_t path_id = 0;
        for (const auto& path : c.paths) {
            std::size_t point_id = 0;

            for (const auto& pt : path) {
                ContourProgressEvent ev;
                ev.type = ContourProgressEventType::PathPoint;
                ev.level = c.level;
                ev.path_id = path_id;
                ev.point_id = point_id++;
                ev.x = pt.first;
                ev.y = pt.second;
                on_progress_(ev);
            }

            ContourProgressEvent ev;
            ev.type = ContourProgressEventType::PathFinished;
            ev.level = c.level;
            ev.path_id = path_id;
            ev.n_points = path.size();
            on_progress_(ev);

            ++path_id;
        }
    }

    return c;
}
