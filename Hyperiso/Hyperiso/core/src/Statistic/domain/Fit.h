#ifndef __FIT_H__
#define __FIT_H__

#include <vector>
#include <functional>
#include "BaseLikelihood.h"
#include "IProfilingStrategy.h"
#include "IContourExtractor.h"
#include "Math.h"

class MLFitter {
public:
    MLFitter(const LikelihoodContext& ctx, const ModelFn& model);

    FitResult maximum_likelihood_fit(const Vector& p0) const;
    Contour contour(double z, std::array<double, 4> bounds) const;

private:
    std::unique_ptr<BaseLikelihood> like_;
};

#endif // __FIT_H__
