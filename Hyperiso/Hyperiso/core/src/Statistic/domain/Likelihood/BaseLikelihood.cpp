#include "BaseLikelihood.h"

BaseLikelihood::BaseLikelihood(const ModelFn& model, LikelihoodContext ctx, size_t p_dim) : p_dim(p_dim) {
    this->ctx = std::move(ctx);
    this->model = model;
    this->dim_ = ctx.nuisance_dist->dim() + p_dim;
}

double BaseLikelihood::nll(const Vector& theta) const {
    Vector p(theta.begin(), theta.begin() + p_dim);
    Vector eta(theta.begin() + p_dim, theta.end());

    Vector res = model(p, eta);
    for (size_t i = 0; i < res.size(); i++)
        res[i] -= this->ctx.exp_obs_values[i];
    
    double ell_obs = this->ctx.exp_obs_dist->logpdf(res);
    double ell_nuis = this->ctx.nuisance_dist->logpdf(eta);
 
    return -(ell_obs + ell_nuis);
}
