#ifndef __AMSCONTOUREXTRACTOR_H__
#define __AMSCONTOUREXTRACTOR_H__

#include "IContourExtractor.h"
#include "ContourObserver.h"
#include "Math.h"

class AMSContourExtractor: public IContourExtractor {
public:
    explicit AMSContourExtractor(ContourProgressCallback on_progress = {})
        : on_progress_(std::move(on_progress)) {}
    Contour extract(const ScalarField2D& field, const ContourRequest& cr) override;
private:
    ContourProgressCallback on_progress_;
};

#endif // __AMSCONTOUREXTRACTOR_H__
