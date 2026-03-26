#ifndef __FIT_H__
#define __FIT_H__

#include <vector>
#include <functional>
#include "BaseLikelihood.h"
#include "IProfilingStrategy.h"
#include "IContourExtractor.h"
#include "Math.h"
#include "ContourEngine.h"

struct ContourOptions {
    ProfilingMethod profiling_method = ProfilingMethod::SLICE;
    ContourAlgorithm primary_contour_method = ContourAlgorithm::MINUIT;
    std::optional<ContourAlgorithm> fallback_contour_method;
    std::size_t resolution = 100;
};

class MLFitter {
public:
    MLFitter(std::shared_ptr<LikelihoodContext> ctx, const ModelFn& model);

    FitResult maximum_likelihood_fit(const std::vector<double>& p0);
    Contour contour(std::size_t x_id, std::size_t y_id, double z, std::array<double, 4> bounds, ContourOptions options) const;

private:
    std::shared_ptr<BaseLikelihood> like_;
    FitResult master_fit_result;
    bool master_fit_success = false;
};

#endif // __FIT_H__
