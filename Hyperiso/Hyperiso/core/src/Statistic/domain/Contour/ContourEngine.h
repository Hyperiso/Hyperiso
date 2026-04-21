#ifndef __CONTOURENGINE_H__
#define __CONTOURENGINE_H__

#include <chrono>

#include "ILikelihood.h"
#include "IProfilingStrategy.h"
#include "IContourExtractor.h"
#include "ProfiledLikelihood2D.h"
#include "Profiler.h"
#include "JointDistribution.h"
#include "WithGaussianConstraints.h"
#include "WithFallback.h"
#include "AMSContourExtractor.h"
#include "MnContourExtractor.h"
#include "GaussianMarginal.h"
#include "GaussianCopula.h"
#include "ContourObserver.h"

enum class ProfilingMethod {
    SLICE,
    FREE_PROJECTION,
    PRIOR_CONSTRAINED_PROJECTION
};

enum class ContourAlgorithm {
    AMS,
    MINUIT
};

struct ContourConfig {
    std::size_t x_id, y_id;
    FitResult fr;
    ProfilingMethod profiling_method;
    ContourAlgorithm primary_contour_method;
    std::optional<ContourAlgorithm> fallback_contour_method;

    ContourProgressCallback on_progress {};
};

class ContourEngine {
public:
    ContourEngine(std::shared_ptr<ILikelihood> base, const ContourConfig& cfg);
    Contour compute_contour(double z, std::array<double, 4> bounds, std::size_t resolution);

private:
    std::shared_ptr<JointDistribution> build_constraints_distribution();
    std::shared_ptr<IContourExtractor> build_contour_extractor(ContourAlgorithm ca);

    ContourConfig cfg;
    ProfiledLikelihood2D likelihood;
    std::shared_ptr<IContourExtractor> extractor;
};

#endif // __CONTOURENGINE_H__
