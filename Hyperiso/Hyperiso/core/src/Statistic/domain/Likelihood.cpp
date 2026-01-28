#include "Likelihood.h"

ProfiledLikelihood::ProfiledLikelihood(LikelihoodContext ctx, ModelFn model) : ctx_(std::move(ctx)), model_(std::move(model)) {
    update_step_sizes();
    min_ctx_.max_iter = 500;
    min_ctx_.tol = 1e-6;
}

// double ProfiledLikelihood::nll(const Vector &p, const Vector &eta) const {
//     Vector r = residual_obs(p, eta);
//     double ell_obs = ctx_.exp_obs_dist->logpdf(r);
//     double ell_eta = ctx_.nuisance_dist->logpdf(eta);

//     // printVector(r);
//     // printf("ll_obs = %.5e\n", ell_obs);
//     // printf("ll_eta = %.5e\n", ell_eta);

//     Vector pred1 = model_(p, eta);
//     Vector pred2 = model_(p, eta);

//     double md = 0.0;
//     for (size_t i = 0; i < pred1.size(); ++i) md = std::max(md, std::abs(pred1[i] - pred2[i]));
//     std::cout << std::setprecision(17) << "[DBG] max |pred1-pred2| = " << md << "\n";

//     static int k = 0;
//     if (k < 3) {
//         std::cout << std::setprecision(17)
//                   << "[DBG] call " << k
//                   << " r0=" << (r.empty()?0.0:r[0])
//                   << " ell_obs=" << ell_obs
//                   << " ell_eta=" << ell_eta
//                   << " nll=" << -(ell_obs + ell_eta)
//                   << "\n";
//         k++;
//     }
    
//     return -(ell_obs + ell_eta);
// }

double ProfiledLikelihood::nll(const Vector &p, const Vector &eta) const {
    Vector pred = model_(p, eta);

    static bool have_prev = false;
    static Vector pred_prev;
    if (have_prev) {
        double md = 0.0;
        size_t imax = 0;
        for (size_t i = 0; i < pred.size(); ++i) {
            double d = std::abs(pred[i] - pred_prev[i]);
            if (d > md) { md = d; imax = i; }
        }
        // std::cout << std::setprecision(17)
        //           << "[DBG] cross-call max|pred-pred_prev|=" << md
        //           << " at i=" << imax
        //           << " pred=" << pred[imax]
        //           << " prev=" << pred_prev[imax]
        //           << "\n";
    }
    pred_prev = pred;
    have_prev = true;

    // check exp_obs_values stability too
    static bool have_exp_prev = false;
    static Vector exp_prev;
    if (!have_exp_prev) { exp_prev = ctx_.exp_obs_values; have_exp_prev = true; }
    else {
        double md = 0.0;
        for (size_t i = 0; i < exp_prev.size(); ++i)
            md = std::max(md, std::abs(exp_prev[i] - ctx_.exp_obs_values[i]));
        if (md != 0.0) {
            std::cout << std::setprecision(17)
                      << "[DBG] exp_obs_values CHANGED! maxdiff=" << md << "\n";
        }
    }

    Vector r(pred.size());
    for (size_t i = 0; i < pred.size(); ++i) r[i] = pred[i] - ctx_.exp_obs_values[i];

    double ell_obs = ctx_.exp_obs_dist->logpdf(r);
    double ell_eta = ctx_.nuisance_dist->logpdf(eta);

    return -(ell_obs + ell_eta);
}

double ProfiledLikelihood::nll_profiled(const Vector &p) const {
    auto f = [this, p] (Vector eta) -> double {
        double v = nll(p, eta);
        return std::isfinite(v) ? v : 1e300;
    };
    return minimize_scaled(f, ctx_.nuisance_central_values, min_ctx_).min;
    // return minimize(f, ctx_.nuisance_central_values, min_ctx_).min;
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
        const double rel = 0.5 * std::abs(x0);
        const double floor_abs = 1000.0 * std::numeric_limits<double>::epsilon() * (std::abs(x0) + 1.0);
        return std::max(rel, floor_abs);
    };

    for (size_t i = 0; i < dim; i++) {
        step_sizes[i] = step_from(ctx_.nuisance_central_values[i]);
    }

    min_ctx_.step_sizes = step_sizes;
}
