#include "Likelihood.h"

ProfiledLikelihood::ProfiledLikelihood(LikelihoodContext ctx, ModelFn model) : ctx_(std::move(ctx)), model_(std::move(model)) {
    update_step_sizes();
    min_ctx_.max_iter = 500;
    min_ctx_.tol = 1e-3;
}

double ProfiledLikelihood::nll(const Vector &p, const Vector &eta) const {
    Vector r = residual_obs(p, eta);
    double ell_obs = ctx_.exp_obs_dist->logpdf(r);
    double ell_eta = ctx_.nuisance_dist->logpdf(eta);
    return ell_obs + ell_eta;
}

double ProfiledLikelihood::nll_profiled(const Vector &p) const {
    auto f = [this, p] (Vector eta) -> double {
        double v = nll(p, eta);
        return std::isfinite(v) ? v : 1e300;
    };

    return minimize(f, ctx_.nuisance_central_values, min_ctx_).min;
}

void ProfiledLikelihood::set_minimizer_max_iter(std::size_t max_iter) {
    if (max_iter < 1) 
        LOG_ERROR("Invalid Argument", "Maximum number of iterations should be at least 1.");

    this->min_ctx_.max_iter = max_iter;
}

void ProfiledLikelihood::set_minimizer_tolerance(double tol) {
    if (tol < 0 || tol > 1) 
        LOG_ERROR("Invalid Argument", "Tolerance should be between 0 and 1.");

    this->min_ctx_.tol = tol;
}

Vector ProfiledLikelihood::get_eta_steps() const {
    return min_ctx_.step_sizes;
}

Vector ProfiledLikelihood::get_eta_central_values() const {
    return ctx_.nuisance_central_values;
}

void ProfiledLikelihood::update_step_sizes() {
    std::size_t dim = ctx_.nuisance_central_values.size();
    Vector step_sizes = Vector(dim);

    auto step_from = [] (double x0) {
        const double rel = 0.05 * std::abs(x0);
        const double floor_abs = 1000.0 * std::numeric_limits<double>::epsilon() * (std::abs(x0) + 1.0);
        return std::max(rel, floor_abs);
    };

    for (size_t i = 0; i < dim; i++) {
        step_sizes[i] = step_from(ctx_.nuisance_central_values[i]);
    }

    min_ctx_.step_sizes = step_sizes;
}
