#pragma once
#include <vector>
#include <functional>
#include "LinearAlgebra.h"


struct LikelihoodContext {
    Vec Oexp; // experimental central values
    SPDMatrix SigmaO; // Cholesky of experimental covariance
    Vec eta_bar; // nuisance prior mean
    SPDMatrix SigmaEta; // Cholesky of nuisance covariance
};


// ℓ(p,η) = (O(p,η)-Oexp)^T Σ_O^{-1} (O(p,η)-Oexp) + (η-η̄)^T Σ_η^{-1} (η-η̄)
class ProfiledLikelihood {
public:
    using ModelFn = std::function<Vec(const Vec& p, const Vec& eta)>;


    explicit ProfiledLikelihood(LikelihoodContext ctx, ModelFn model)
    : ctx_(std::move(ctx)), model_(std::move(model)) {}


    double ell(const Vec& p, const Vec& eta) const {
        Vec r = residual_obs(p, eta);
        double q_obs = ctx_.SigmaO.quad_inv(r);
        Vec d = diff_eta(eta);
        double q_eta = ctx_.SigmaEta.quad_inv(d);
        return q_obs + q_eta;
    }


    // profile ℓ over η for fixed p
    double ell_profiled(const Vec& p, const Vec& eta0, std::size_t max_iter=500, double tol=1e-6) const;


private:
    Vec residual_obs(const Vec& p, const Vec& eta) const {
        Vec pred = model_(p, eta);
        if (pred.size()!=ctx_.Oexp.size()) throw std::invalid_argument("Model produced wrong dimension");
        Vec r(pred.size());
        for (std::size_t i=0;i<pred.size();++i) r[i] = pred[i]-ctx_.Oexp[i];
        return r;
    }
    Vec diff_eta(const Vec& eta) const {
        if (eta.size()!=ctx_.eta_bar.size()) throw std::invalid_argument("eta size mismatch");
        Vec d(eta.size());
        for (std::size_t i=0;i<eta.size();++i) d[i] = eta[i]-ctx_.eta_bar[i];
        return d;
    }


    LikelihoodContext ctx_;
    ModelFn model_;
};