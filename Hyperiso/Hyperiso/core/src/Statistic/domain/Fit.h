// #ifndef __FIT_H__
// #define __FIT_H__

// #include <vector>
// #include <functional>
// #include "BaseLikelihood.h"
// #include "IProfilingStrategy.h"
// #include "IContourExtractor.h"
// #include "Math.h"
// #include "ContourEngine.h"

// struct ContourOptions {
//     ProfilingMethod profiling_method = ProfilingMethod::SLICE;
//     ContourAlgorithm primary_contour_method = ContourAlgorithm::MINUIT;
//     std::optional<ContourAlgorithm> fallback_contour_method;
//     std::size_t resolution = 100;
// };

// class MLFitter {
// public:
//     MLFitter(std::shared_ptr<LikelihoodContext> ctx, const ModelFn& model);

//     FitResult maximum_likelihood_fit(const std::vector<double>& p0);
//     Contour contour(std::size_t x_id, std::size_t y_id, double z, std::array<double, 4> bounds, ContourOptions options) const;

// private:
//     std::shared_ptr<BaseLikelihood> like_;
//     FitResult master_fit_result;
//     bool master_fit_success = false;
// };

// #endif // __FIT_H__

#ifndef __FIT_H__
#define __FIT_H__

#include <vector>
#include <functional>
#include "BaseLikelihood.h"
#include "IProfilingStrategy.h"
#include "IContourExtractor.h"
#include "Math.h"
#include "ContourEngine.h"
#include "ContourObserver.h"


struct MLFitOptions {
    bool run_hesse = true;
    bool request_minos = false;
    bool verbose = false;
    unsigned strategy = 0;      // 0 => backend default
    unsigned max_fcn = 0;       // 0 => backend default
    double tolerance = 0.0;     // <=0 => backend default

    // Fallback si HESSE/Minuit ne donne pas de covariance
    bool allow_profile_hessian_fallback = true;
    double profile_hessian_step_scale = 1.0;
    double profile_hessian_eig_floor_rel = 1e-8;

    bool trace_first_evals = false;
    std::size_t trace_max_evals = 25;
};

// struct MLFitOptions {
//     bool run_hesse = true;
//     bool request_minos = false;
//     bool verbose = false;
//     unsigned strategy = 0;      // 0 => backend default
//     unsigned max_fcn = 0;       // 0 => backend default
//     double tolerance = 0.0;     // <=0 => backend default
// };

enum class ProfileBackend {
    MINUIT,
    LAPLACE_NUISANCE
};

struct ContourOptions {
    ProfilingMethod profiling_method = ProfilingMethod::SLICE;
    ProfileBackend profile_backend = ProfileBackend::LAPLACE_NUISANCE;
    ContourAlgorithm primary_contour_method = ContourAlgorithm::MINUIT;
    std::optional<ContourAlgorithm> fallback_contour_method;
    std::size_t resolution = 40;

    ContourProgressCallback on_progress {};
};

class MLFitter {
public:
    MLFitter(std::shared_ptr<LikelihoodContext> ctx, const ModelFn& model, MLFitOptions options = {});
    explicit MLFitter(std::shared_ptr<BaseLikelihood> like, MLFitOptions options = {});
    
    FitResult maximum_likelihood_fit(const std::vector<double>& p0);
    Contour contour(std::size_t x_id, std::size_t y_id, double z, std::array<double, 4> bounds, ContourOptions options) const;

private:
    std::shared_ptr<BaseLikelihood> like_;
    MLFitOptions fit_options_;
    FitResult master_fit_result;
    bool master_fit_success = false;
};

#endif // __FIT_H__
