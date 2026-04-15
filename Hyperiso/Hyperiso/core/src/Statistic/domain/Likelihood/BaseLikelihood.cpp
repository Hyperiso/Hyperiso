#include "BaseLikelihood.h"

BaseLikelihood::BaseLikelihood(const ModelFn& model, std::shared_ptr<LikelihoodContext> ctx, size_t p_dim) : p_dim(p_dim) {
    this->ctx = std::move(ctx);
    this->model = model;
}

// double BaseLikelihood::nll(const std::vector<double>& theta) const {
//     std::vector<double>p(theta.begin(), theta.begin() + p_dim);
//     std::vector<double>eta(theta.begin() + p_dim, theta.end());

//     std::vector<double>res = model(p, eta);
//     for (size_t i = 0; i < res.size(); i++)
//         res[i] -= this->ctx->exp_obs_values[i];
    
//     double ell_obs = this->ctx->exp_obs_dist->logpdf(res);
//     double ell_nuis = this->ctx->nuisance_dist->logpdf(eta);
 
//     return -(ell_obs + ell_nuis);
// }

// std::size_t BaseLikelihood::dim() const {
//     return this->ctx->nuisance_dist->dim() + this->p_dim;
// }


void BaseLikelihood::maybe_log_debug_eval(const std::vector<double>& theta,
                                          const std::vector<double>& p,
                                          const std::vector<double>& eta,
                                          const std::vector<double>& res,
                                          double ell_obs,
                                          double ell_nuis,
                                          double nll_value) const
{
    if (!debug_trace_enabled_) {
        return;
    }
    if (debug_eval_count_ >= debug_trace_max_evals_) {
        return;
    }

    if (!debug_have_ref_theta_) {
        debug_ref_theta_ = theta;
        debug_have_ref_theta_ = true;
    }
    if (!debug_have_ref_res_) {
        debug_ref_res_ = res;
        debug_have_ref_res_ = true;
    }

    const double theta_shift = max_abs_diff(theta, debug_ref_theta_);
    const double res_shift   = max_abs_diff(res, debug_ref_res_);

    std::cout << "[FITDBG] eval=" << debug_eval_count_
              << " nll=" << nll_value
              << " ell_obs=" << ell_obs
              << " ell_nuis=" << ell_nuis
              << " max|theta-theta0|=" << theta_shift
              << " max|res-res0|=" << res_shift;

    if (!p.empty()) {
        std::cout << " p=[";
        for (std::size_t i = 0; i < std::min<std::size_t>(p.size(), 4); ++i) {
            if (i) std::cout << ", ";
            std::cout << p[i];
        }
        if (p.size() > 4) std::cout << ", ...";
        std::cout << "]";
    }

    if (!eta.empty()) {
        std::cout << " eta=[";
        for (std::size_t i = 0; i < std::min<std::size_t>(eta.size(), 4); ++i) {
            if (i) std::cout << ", ";
            std::cout << eta[i];
        }
        if (eta.size() > 4) std::cout << ", ...";
        std::cout << "]";
    }

    std::cout << std::endl;
    ++debug_eval_count_;
}

double BaseLikelihood::nll(const std::vector<double>& theta) const {
    std::vector<double> p(theta.begin(), theta.begin() + p_dim);
    std::vector<double> eta(theta.begin() + p_dim, theta.end());

    std::vector<double> res = model(p, eta);
    for (size_t i = 0; i < res.size(); i++) {
        res[i] -= this->ctx->exp_obs_values[i];
    }

    double ell_obs  = this->ctx->exp_obs_dist->logpdf(res);
    double ell_nuis = this->ctx->nuisance_dist->logpdf(eta);
    double out = -(ell_obs + ell_nuis);

    maybe_log_debug_eval(theta, p, eta, res, ell_obs, ell_nuis, out);
    return out;
}

std::size_t BaseLikelihood::dim() const {
    return this->ctx->nuisance_dist->dim() + this->p_dim;
}