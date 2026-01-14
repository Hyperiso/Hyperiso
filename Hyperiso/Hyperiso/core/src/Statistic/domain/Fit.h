#pragma once
#include <vector>
#include <functional>
#include "Likelihood.h"

struct FitResult {
    Vector p_hat; // MLE estimators
    Vector eta_hat; // profiled-at-MLE nuisances
    double ell_hat {0.0}; // minimum NLL
};

class MLEstimator {
public:
    using ModelFn = ProfiledLikelihood::ModelFn;

    MLEstimator(LikelihoodContext ctx, ModelFn model, std::size_t max_iter=500, double tol=1e-6)
        : like_(std::move(ctx), std::move(model)), max_iter(max_iter), tol(tol) {}

    FitResult fit(const Vector& p0) const;

    double wilks_T(const Vector& p, const FitResult& fr) const {
        double lprof = like_.nll_profiled(p);
        return lprof - fr.ell_hat;
    }

    // 1D scan: find interval {p | T(p) ≤ χ²_{1−α, dof}}, using coarse grid then refine
    // std::pair<double,double> confidence_interval_1d(
    //     const FitResult& fr,
    //     double p_min, double p_max,
    //     int grid_points,
    //     double alpha, // e.g. 0.05 → 95%
    //     const Vector& eta_init
    // ) const;

    const ProfiledLikelihood& like() const { return like_; }

private:
    ProfiledLikelihood like_;
    std::size_t max_iter;
    double tol;
};