#include "Fit.h"

MLFitter::MLFitter(std::shared_ptr<LikelihoodContext> ctx, const ModelFn& model) {
    const std::size_t p_dim = ctx->fp_defs.size();
    this->like_ = std::make_shared<BaseLikelihood>(model, ctx, p_dim);
}

FitResult MLFitter::maximum_likelihood_fit(const std::vector<double>& p0) {
    const auto defs = like_->get_param_defs();
    const std::size_t p_dim = p0.size();
    const std::size_t dim = defs.size();

    std::vector<fit_app::ParameterDefinition> theta0 = defs;
    for (std::size_t i = 0; i < p_dim; ++i) {
        theta0[i].value = p0[i];
    }

    auto f = fit_app::LambdaObjectiveFunction(
        [this](const std::vector<double>& theta) {
            return like_->nll(theta);
        },
        0.5
    );

    fit_app::FitOptions opt;
    opt.run_hesse = true;
    opt.verbose = false;

    std::unique_ptr<fit_app::IFitBackend> minimizer = fit_app::make_minuit_backend();
    fit_app::BackendFitResult res = minimizer->minimize(f, theta0, opt);

    FitResult fr;
    fr.ell_hat = res.diagnostics.fmin;
    fr.p_hat.assign(res.values.begin(), res.values.begin() + p_dim);
    fr.eta_hat.assign(res.values.begin() + p_dim, res.values.end());

    // bloc p×p de la covariance complète
    RealMatrix H = res.covariance.inv();
    RealMatrix H_p_p(p_dim, p_dim);
    RealMatrix H_p_eta(p_dim, dim - p_dim);
    RealMatrix H_eta_eta(dim - p_dim, dim - p_dim);

    for (std::size_t i = 0; i < p_dim; ++i) {
        for (std::size_t j = 0; j < p_dim; ++j) {
            H_p_p.at(i, j) = H.at(i, j);
        }
    }

    for (std::size_t i = 0; i < p_dim; ++i) {
        for (std::size_t j = p_dim; j < dim; ++j) {
            H_p_eta.at(i, j - p_dim) = H.at(i, j);
        }
    }

    for (std::size_t i = p_dim; i < dim; ++i) {
        for (std::size_t j = p_dim; j < dim; ++j) {
            H_eta_eta.at(i - p_dim, j - p_dim) = H.at(i, j);
        }
    }

    RealMatrix H_prof = H_p_p - H_p_eta * H_eta_eta.inv() * H_p_eta.transpose();
    RealMatrix cov_prof = H_prof.inv();

    fr.p_hat_std.resize(p_dim);
    fr.p_hat_correlations = RealMatrix(p_dim, p_dim);

    for (std::size_t i = 0; i < p_dim; ++i) {
        fr.p_hat_std[i] = std::sqrt(std::max(0.0, cov_prof.at(i, i)));
    }
    for (std::size_t i = 0; i < p_dim; ++i) {
        for (std::size_t j = 0; j < p_dim; ++j) {
            const double si = fr.p_hat_std[i];
            const double sj = fr.p_hat_std[j];
            fr.p_hat_correlations.at(i, j) =
                (si > 0.0 && sj > 0.0) ? cov_prof.at(i, j) / (si * sj) : 0.0;
        }
    }

    this->master_fit_success = res.diagnostics.ok;
    this->master_fit_result = fr;
    return fr;
}

Contour MLFitter::contour(std::size_t x_id, std::size_t y_id, double z, std::array<double, 4> bounds, ContourOptions options) const {
    if (!this->master_fit_success)
        LOG_ERROR("InvalidState", "ML fit must have converged before contour computation is available.");

    ContourConfig cc;
    cc.fr = this->master_fit_result;
    cc.x_id = x_id;
    cc.y_id = y_id;
    cc.primary_contour_method = options.primary_contour_method;
    cc.fallback_contour_method = options.fallback_contour_method;
    cc.profiling_method = options.profiling_method;

    ContourEngine ce(this->like_, cc);
    Contour cl = ce.compute_contour(z, bounds, options.resolution);

    return cl;
}