#pragma once
#include <vector>
#include <functional>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_cdf.h>
#include "Likelihood.h"


struct FitResult {
    Vec p_hat; // MLE parameters
    Vec eta_hat; // profiled-at-MLE nuisances
    double ell_hat{0.0};
};


class MLEstimator {
public:
    using ModelFn = ProfiledLikelihood::ModelFn;


    MLEstimator(LikelihoodContext ctx, ModelFn model)
        : like_(std::move(ctx), std::move(model)) {}


    FitResult fit(const Vec& p0, const Vec& eta0, std::size_t max_iter=1000, double tol=1e-6) const;


    // T(p) = \tilde{ℓ}(p) − ℓ(\^p, \^η)
    double test_statistic(const Vec& p, const FitResult& fr, const Vec& eta_init) const {
        double lprof = like_.ell_profiled(p, eta_init);
        return lprof - fr.ell_hat;
    }


    // 1D scan: find interval {p | T(p) ≤ χ²_{1−α, dof}}, using coarse grid then refine
    std::pair<double,double> confidence_interval_1d(
        const FitResult& fr,
        double p_min, double p_max,
        int grid_points,
        double alpha, // e.g. 0.05 → 95%
        const Vec& eta_init
    ) const;


    const ProfiledLikelihood& like() const { return like_; }


private:
    ProfiledLikelihood like_;
};