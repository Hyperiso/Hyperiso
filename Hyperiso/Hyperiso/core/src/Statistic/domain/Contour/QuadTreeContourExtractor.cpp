#include "QuadTreeContourExtractor.h"

void QuadTreeContourExtractor::extract(std::shared_ptr<Contour> contour) {
    auto f = [contour] () {
        return contour->likelihood->nll();
    };

    MarchingSquaresExtractor mse;
    mse.find_iso_contour(contour->level);
}