#ifndef __MNCONTOUREXTRACTOR_H__
#define __MNCONTOUREXTRACTOR_H__

#include "IContourExtractor.h"

class MnContourExtractor: public IContourExtractor {
public:
    Contour extract(const ScalarField2D& field, const ContourRequest& cr) override;
};

#endif // __MNCONTOUREXTRACTOR_H__
