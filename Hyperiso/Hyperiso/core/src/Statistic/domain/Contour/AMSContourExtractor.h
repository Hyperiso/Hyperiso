#ifndef __AMSCONTOUREXTRACTOR_H__
#define __AMSCONTOUREXTRACTOR_H__

#include "IContourExtractor.h"
#include "Math.h"

class AMSContourExtractor: public IContourExtractor {
public:
    Contour extract(const ScalarField2D& field, const ContourRequest& cr) override;
};

#endif // __AMSCONTOUREXTRACTOR_H__
