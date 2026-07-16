#include "BaseLikelihood.h"

BaseLikelihood::BaseLikelihood(const ModelFn& model, std::shared_ptr<LikelihoodContext> ctx, size_t p_dim) : p_dim(p_dim) {
    this->ctx = std::move(ctx);
    this->model = model;
}


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

    try {
        std::vector<double> res = model(p, eta);

        for (size_t i = 0; i < res.size(); i++) {
            if (!std::isfinite(res[i])) {
                return 1e100;
            }
            res[i] -= this->ctx->exp_obs_values[i];
            if (!std::isfinite(res[i])) {
                return 1e100;
            }
        }

        double ell_obs  = this->ctx->exp_obs_dist->logpdf(res);
        double ell_nuis = this->ctx->nuisance_dist->logpdf(eta);

        if (!std::isfinite(ell_obs) || !std::isfinite(ell_nuis)) {
            return 1e100;
        }

        double out = -(ell_obs + ell_nuis);
        if (!std::isfinite(out)) {
            return 1e100;
        }

        maybe_log_debug_eval(theta, p, eta, res, ell_obs, ell_nuis, out);
        return out;
    }
    catch (const std::exception& e) {
        if (debug_trace_enabled_) {
            std::cerr << "[FIT DEBUG] model/nll exception: " << e.what() << "\n";
        }
        return 1e100;
    }
    catch (...) {
        if (debug_trace_enabled_) {
            std::cerr << "[FIT DEBUG] model/nll unknown exception\n";
        }
        return 1e100;
    }
}

std::size_t BaseLikelihood::dim() const {
    return this->ctx->nuisance_dist->dim() + this->p_dim;
}

std::size_t BaseLikelihood::p_dimension() const {
    return p_dim;
}

std::size_t BaseLikelihood::eta_dimension() const {
    return ctx->nuis_defs.size();
}

std::vector<double> BaseLikelihood::central_p() const {
    std::vector<double> out;
    out.reserve(ctx->fp_defs.size());

    for (const auto& def : ctx->fp_defs) {
        out.push_back(def.value);
    }

    return out;
}

std::vector<double> BaseLikelihood::central_eta() const {
    std::vector<double> out;
    out.reserve(ctx->nuis_defs.size());

    for (const auto& def : ctx->nuis_defs) {
        out.push_back(def.value);
    }

    return out;
}

std::vector<double> BaseLikelihood::predict(
    const std::vector<double>& p,
    const std::vector<double>& eta
) const {
    return model(p, eta);
}

std::vector<double> BaseLikelihood::residuals(
    const std::vector<double>& p,
    const std::vector<double>& eta
) const {
    std::vector<double> pred = model(p, eta);

    if (pred.size() != ctx->exp_obs_values.size()) {
        throw std::runtime_error("BaseLikelihood::residuals: prediction/observation size mismatch");
    }

    for (std::size_t i = 0; i < pred.size(); ++i) {
        pred[i] -= ctx->exp_obs_values[i];
    }

    return pred;
}

double BaseLikelihood::nll_from_split(
    const std::vector<double>& p,
    const std::vector<double>& eta
) const {
    std::vector<double> theta;
    theta.reserve(p.size() + eta.size());
    theta.insert(theta.end(), p.begin(), p.end());
    theta.insert(theta.end(), eta.begin(), eta.end());

    return nll(theta);
}

RealMatrix BaseLikelihood::observable_curvature(
    const std::vector<double>& r
) const {
    RealMatrix W = ctx->exp_obs_dist->curvature(r);

    for (std::size_t i = 0; i < W.rows(); ++i) {
        for (std::size_t j = 0; j < W.cols(); ++j) {
            if (!std::isfinite(W.at(i, j))) {
                std::cout << "[WOBSDBG] non-finite W_obs("
                          << i << "," << j << ")"
                          << " r_i=" << r[i]
                          << " r_j=" << r[j]
                          << " exp_i=" << ctx->exp_obs_values[i]
                          << " exp_j=" << ctx->exp_obs_values[j]
                          << std::endl;
            }
        }
    }

    return W;
}

RealMatrix BaseLikelihood::nuisance_curvature(
    const std::vector<double>& eta
) const {
    return ctx->nuisance_dist->curvature(eta);
}