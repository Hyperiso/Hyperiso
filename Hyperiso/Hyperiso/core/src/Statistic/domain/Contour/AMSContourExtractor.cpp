#include "AMSContourExtractor.h"

Contour AMSContourExtractor::extract(const ScalarField2D& field, const ContourRequest& cr) {
    Contour c;
    c.level = cr.level;

    std::size_t depth = std::ceil(std::log2(cr.resolution));
    MarchingSquaresExtractor mse(field, cr.bounds, depth);
    c.paths = mse.find_iso_contour(cr.level);
    c.success = true;
    
    return c; 
}