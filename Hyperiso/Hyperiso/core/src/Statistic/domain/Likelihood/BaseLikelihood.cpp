#include "BaseLikelihood.h"

BaseLikelihood::BaseLikelihood(const ModelFn& model, std::shared_ptr<LikelihoodContext> ctx, size_t p_dim) : p_dim(p_dim) {
    LOG_INFO("1");
    this->ctx = std::move(ctx);
    this->model = model;
    LOG_INFO("2");
}

double BaseLikelihood::nll(const std::vector<double>& theta) const {
    std::vector<double>p(theta.begin(), theta.begin() + p_dim);
    std::vector<double>eta(theta.begin() + p_dim, theta.end());

    std::vector<double>res = model(p, eta);
    for (size_t i = 0; i < res.size(); i++)
        res[i] -= this->ctx->exp_obs_values[i];
    
    double ell_obs = this->ctx->exp_obs_dist->logpdf(res);
    double ell_nuis = this->ctx->nuisance_dist->logpdf(eta);
 
    return -(ell_obs + ell_nuis);
}

std::size_t BaseLikelihood::dim() const {
    return this->ctx->nuisance_dist->dim() + this->p_dim;
}
