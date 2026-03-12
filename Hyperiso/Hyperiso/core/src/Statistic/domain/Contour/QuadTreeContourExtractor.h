#ifndef __QUADTREECONTOUREXTRACTOR_H__
#define __QUADTREECONTOUREXTRACTOR_H__

#include "IContourExtractor.h"
#include "ILikelihood.h"
#include "Math.h"

class QuadTreeContourExtractor: public IContourExtractor {
public:
    void extract(std::shared_ptr<Contour> contour) override;
};

#endif // __QUADTREECONTOUREXTRACTOR_H__
