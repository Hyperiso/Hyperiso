// #include "AMSContourExtractor.h"

// Contour AMSContourExtractor::extract(const ScalarField2D& field, const ContourRequest& cr) {
//     Contour c;
//     c.level = cr.level;

//     std::size_t depth = std::ceil(std::log2(cr.resolution));
//     MarchingSquaresExtractor mse(field, cr.bounds, depth);
//     c.paths = mse.find_iso_contour(cr.level);
//     c.success = true;
    
//     return c; 
// }

#include "AMSContourExtractor.h"

Contour AMSContourExtractor::extract(const ScalarField2D& field, const ContourRequest& cr) {
    Contour c;
    c.level = cr.level;

    std::size_t depth = std::ceil(std::log2(cr.resolution));
    MarchingSquaresExtractor mse(field, cr.bounds, depth);
    c.paths = mse.find_iso_contour(cr.level);
    c.success = true;

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