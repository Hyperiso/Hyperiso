#include "Fit.h"

static double step_from(double x0) {
    const double rel = 0.01 * std::abs(x0);
    const double floor_abs = 1000.0 * std::numeric_limits<double>::epsilon() * (std::abs(x0) + 1.0);
    return std::max(rel, floor_abs);
}

FitResult MLEstimator::fit(const Vector& p0) const {
    std::size_t p_dim = p0.size();

    auto f = [this, p_dim] (Vector p_and_eta) -> double {
        auto p_first = p_and_eta.begin();
        auto p_last = p_first + p_dim;
        std::vector<double> p (p_first, p_last);
        std::vector<double> eta (p_last, p_and_eta.end());
        double v = this->like_.nll(p, eta);
        return std::isfinite(v) ? v : 1e300;
    };

    MinimizationContext min_ctx;
    min_ctx.max_iter = max_iter;
    min_ctx.tol = tol;
    
    Vector step_sizes;
    for (double p : p0) {
        step_sizes.emplace_back(step_from(p));
    }

    Vector eta_steps = like_.get_eta_steps();
    step_sizes.insert(step_sizes.end(), eta_steps.begin(), eta_steps.end());
    min_ctx.step_sizes = step_sizes;

    Vector start = p0;
    Vector eta_start = like_.get_eta_central_values();
    start.insert(start.end(), eta_start.begin(), eta_start.end());

    MinimizationResult min_res = minimize(f, start, min_ctx);

    auto p_and_eta_hat = min_res.argmin;
    auto first = p_and_eta_hat.begin();
    auto last = first + p_dim;
    std::vector<double> p_hat (first, last);
    std::vector<double> eta_hat (last, p_and_eta_hat.end());

    return FitResult {p_hat, eta_hat, min_res.min};
}
