#include "Fit.h"

static double step_from(double x0) {
    const double rel = 0.1 * std::abs(x0);
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

    // {
    //     Vector eta0 = like_.get_eta_central_values();
    //     Vector p = p0;

    //     double f0 = like_.nll(p, eta0);

    //     double dp = step_from(p0[0]);
    //     p[0] = p0[0] + dp;
    //     double f_plus = like_.nll(p, eta0);

    //     p[0] = p0[0] - dp;
    //     double f_minus = like_.nll(p, eta0);

    //     p[0] = p0[0] + 100.0 * dp;
    //     double f_big = like_.nll(p, eta0);

    //     std::cout << "[DBG] nll(p0,eta0)=" << f0
    //             << " nll(p0+dp)=" << f_plus
    //             << " nll(p0-dp)=" << f_minus
    //             << " nll(p0+100dp)=" << f_big
    //             << " (dp=" << dp << ")"
    //             << std::endl;
    // }

    // --- DEBUG scan profiled ---
    // {
    //     Vector p = p0;
    //     double dp = step_from(p0[0]);

    //     for (int k = -5; k <= 5; ++k) {
    //         p[0] = p0[0] + k * dp;
    //         double v = like_.nll_profiled(p);
    //         std::cout << "[DBG] k=" << k << " p=" << p[0] << " nll_prof=" << v << "\n";
    //     }
    // }

    auto p = p0;
    auto eta = like_.get_eta_central_values();
    double a = like_.nll(p, eta);
    double b = like_.nll(p, eta);
    std::cout << std::setprecision(17)
            << "nll repeat diff = " << (a-b) << "\n";


    {
        std::vector<double> p(p0.begin(), p0.end());
        std::vector<double> eta = like_.get_eta_central_values();

        double n1 = like_.nll(p, eta);
        double n2 = like_.nll(p, eta);

        std::cout << std::setprecision(17)
                << "[DBG] fixed-point nll1=" << n1
                << " nll2=" << n2
                << " diff=" << (n1 - n2) << "\n";
    }

    // MinimizationResult min_res = minimize(f, start, min_ctx);
    MinimizationResult min_res = minimize_scaled(f, start, min_ctx);

    auto p_and_eta_hat = min_res.argmin;
    auto first = p_and_eta_hat.begin();
    auto last = first + p_dim;
    std::vector<double> p_hat (first, last);
    std::vector<double> eta_hat (last, p_and_eta_hat.end());

    double f_start = f(start);
    double f_end   = min_res.min;
    std::cout << "[DBG] f_start=" << f_start << " f_end=" << f_end << std::endl;

    return FitResult {p_hat, eta_hat, min_res.min}; 
}
