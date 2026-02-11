#pragma once
#include <vector>
#include <functional>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_erf.h>
#include "Likelihood.h"
#include "contour.h"

struct FitResult {
    Vector p_hat; // MLE estimators
    Vector eta_hat; // profiled-at-MLE nuisances
    Vector p_hat_std;
    RealMatrix p_hat_correlations;
    double ell_hat {0.0}; // minimum NLL
};

class MLEstimator {
public:
    using ModelFn = ProfiledLikelihood::ModelFn;

    MLEstimator(LikelihoodContext ctx, ModelFn model, std::size_t max_iter=500, double tol=1e-6)
        : like_(std::move(ctx), std::move(model)), max_iter(max_iter), tol(tol) {}

    FitResult fit(const Vector& p0) const;
    std::set<std::vector<std::pair<double, double>>> contour(double z, std::array<double, 4> bounds, double ell_hat) const;

    double wilks_T(const Vector& p, double ell_hat) const {
        double lprof = like_.nll_profiled(p);
        return -2 * (lprof - ell_hat);
    }

    const ProfiledLikelihood& like() const { return like_; }

private:
    ProfiledLikelihood like_;
    std::size_t max_iter;
    double tol;
    int dim;
};