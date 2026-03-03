#include "Likelihood.h"

ProfiledLikelihood::ProfiledLikelihood(LikelihoodContext ctx, ModelFn model) : ctx_(std::move(ctx)), model_(std::move(model)) {
    min_ctx_.final_tol = 1e-8;
    min_ctx_.switch_tol = 1e-3;
    min_ctx_.simplex_initial_step_size = 0.2;
    min_ctx_.simplex_max_iter = 10000; //TODO
}

double ProfiledLikelihood::nll(const Vector &p, const Vector &eta) const {
    Vector pred = model_(p, eta);

    // static bool have_prev = false;
    // static Vector pred_prev;
    // if (have_prev) {
    //     double md = 0.0;
    //     size_t imax = 0;
    //     for (size_t i = 0; i < pred.size(); ++i) {
    //         double d = std::abs(pred[i] - pred_prev[i]);
    //         if (d > md) { md = d; imax = i; }
    //     }
    // }
    // pred_prev = pred;
    // have_prev = true;

    // static bool have_exp_prev = false;
    // static Vector exp_prev;
    // if (!have_exp_prev) { exp_prev = ctx_.exp_obs_values; have_exp_prev = true; }
    // else {
    //     double md = 0.0;
    //     for (size_t i = 0; i < exp_prev.size(); ++i)
    //         md = std::max(md, std::abs(exp_prev[i] - ctx_.exp_obs_values[i]));
    //     if (md != 0.0) {
    //         std::cout << std::setprecision(17)
    //                   << "[DBG] exp_obs_values CHANGED! maxdiff=" << md << "\n";
    //     }
    // }

    Vector r(pred.size());
    for (size_t i = 0; i < pred.size(); ++i) r[i] = pred[i] - ctx_.exp_obs_values[i];

    Vector r_eta(eta.size());
    for (size_t i = 0; i < eta.size(); ++i) r_eta[i] = eta[i] - ctx_.nuisance_central_values[i];

    double ell_obs = ctx_.exp_obs_dist->logpdf(r);
    double ell_eta = ctx_.nuisance_dist->logpdf(eta);

    // DEBUG
    // Vector eta_ (eta);
    // for (size_t i = 0; i < eta.size(); i++) {
    //     eta_ = eta;
    //     eta_[i] += 0.2 * get_eta_standard_devs()[i];
    //     Vector pred = model_(p, eta_);
    //     Vector r(pred.size());
    //     for (size_t j = 0; j < pred.size(); ++j) r[j] = pred[j] - ctx_.exp_obs_values[j];
    //     double d_obs = std::abs((ctx_.exp_obs_dist->logpdf(r) - ell_obs));
    //     double d_nuis = std::abs((ctx_.nuisance_dist->logpdf(eta_) - ell_eta));
    //     std::cout << "d ell_obs = " << d_obs << std::endl;
    //     std::cout << "d ell_eta = " << d_nuis << std::endl;
    //     std::cout << "d ell_obs / ell_obs = " << d_obs / ell_obs << std::endl;
    //     std::cout << "d ell_eta / ell_eta = " << d_nuis / ell_eta << std::endl;
    // }

    // std::cout << '\n';

    // std::cout << "Residues on obs = [";
    // for (double res : r) std::cout << res << ", ";
    // std::cout << "]" << std::endl;
    // std::cout << "Residues on nuis = [";
    // for (double r_ : r_eta) std::cout << r_ << ", ";
    // std::cout << "]" << std::endl;
    // std::cout << "ell_obs = " << ell_obs << std::endl;
    // std::cout << "ell_nuis = " << ell_eta << std::endl;

    double ell_eta_wrong = ctx_.nuisance_dist->logpdf(eta);
double ell_eta_right = ctx_.nuisance_dist->logpdf(r_eta);

// static int c = 0;
// if (c++ < 5) {
//     std::cout << std::setprecision(17);
//     std::cout << "[TEST B] ell_eta_wrong=" << ell_eta_wrong
//               << " ell_eta_right=" << ell_eta_right << "\n";
// }

//     static bool once = true;
// if (once) {
//     once = false;
//     auto eta0 = ctx_.nuisance_central_values;
//     auto eta_zero = Vector(eta0.size(), 0.0);

//     auto r0 = Vector(eta0.size());
//     for (size_t i=0;i<eta0.size();++i) r0[i] = eta0[i] - eta0[i]; // =0

//     std::cout << std::setprecision(17);
//     std::cout << "[TEST A] logpdf(eta0)  = " << ctx_.nuisance_dist->logpdf(eta0) << "\n";
//     std::cout << "[TEST A] logpdf(0)     = " << ctx_.nuisance_dist->logpdf(eta_zero) << "\n";
//     std::cout << "[TEST A] logpdf(r_eta0)= " << ctx_.nuisance_dist->logpdf(r0) << " (r_eta0 should be 0)\n";
// }

    return -(ell_obs + ell_eta);
}

double ProfiledLikelihood::nll_profiled(const Vector &p) const {
    auto f = [this, p] (Vector eta) -> double {
        double v = nll(p, eta);
        return std::isfinite(v) ? v : 1e300;
    };

    MinimizationResult min_res = minimize_combined(f, ctx_.nuisance_central_values, get_eta_standard_devs(), min_ctx_);

    std::cout << std::setprecision(17);
    for (auto v : ctx_.nuisance_central_values) std::cout << v << " ";
    std::cout << std::endl;
    for (auto v : min_res.argmin) std::cout << v << " ";
    std::cout << std::endl;
    return min_res.min;
}

void ProfiledLikelihood::set_minimizer_max_iter(std::size_t max_iter) {
    if (max_iter < 1) 
        LOG_ERROR("Invalid Argument", "Maximum number of iterations should be at least 1.");

    // TODO : separate simplex and bfgs iterations
    this->min_ctx_.bfgs_max_iter = max_iter;
}

void ProfiledLikelihood::set_minimizer_tolerance(double tol) {
    if (tol < 0 || tol > 1) 
        LOG_ERROR("Invalid Argument", "Tolerance should be between 0 and 1.");

    this->min_ctx_.final_tol = tol;
}

Vector ProfiledLikelihood::get_eta_central_values() const {
    return ctx_.nuisance_central_values;
}

Vector ProfiledLikelihood::get_eta_standard_devs() const {
    return ctx_.nuisance_dist->get_stds();
}
