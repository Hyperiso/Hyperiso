#include "Profiler.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

static bool acceptable_for_profile(const fit_app::BackendFitResult& r) {
    if (!r.diagnostics.has_valid_parameters) return false;
    if (!std::isfinite(r.diagnostics.fmin)) return false;
    if (r.diagnostics.reached_call_limit) return false;

    if (r.diagnostics.ok) return true;

    if (!r.diagnostics.hesse_failed &&
        r.diagnostics.has_valid_covar &&
        r.diagnostics.has_posdef_covar &&
        r.diagnostics.edm < 1.0) {
        return true;
    }

    return false;
}

static ProfileResult full_theta_result_from_split(
    const std::vector<double>& p,
    const std::vector<double>& eta,
    double nll_hat,
    bool ok
) {
    ProfileResult out;
    out.nll_hat = nll_hat;
    out.converged = ok;

    for (std::size_t i = 0; i < p.size(); ++i) {
        out.theta_hat[i] = p[i];
    }

    for (std::size_t a = 0; a < eta.size(); ++a) {
        out.theta_hat[p.size() + a] = eta[a];
    }

    return out;
}

static ProfileResult direct_profile_no_free(
    const std::shared_ptr<ILikelihood>& base,
    const ProfileRequest& pr
) {
    ProfileResult out;
    out.converged = false;
    out.nll_hat = 1e300;

    std::vector<double> theta = pr.start;
    if (theta.empty()) {
        theta.assign(base->dim(), 0.0);
    }
    if (theta.size() != base->dim()) {
        throw std::runtime_error("direct_profile_no_free: bad start dimension");
    }

    for (const auto& [idx, val] : pr.fixed_params) {
        if (idx >= theta.size()) {
            throw std::runtime_error("direct_profile_no_free: fixed index out of range");
        }
        theta[idx] = val;
    }

    out.nll_hat = base->nll(theta);
    out.converged = std::isfinite(out.nll_hat);

    for (std::size_t i = 0; i < theta.size(); ++i) {
        out.theta_hat[i] = theta[i];
    }

    return out;
}

} // namespace

Profiler::Profiler(
    std::shared_ptr<fit_app::IFitBackend> minimizer,
    ProfilerMode mode
)
    : minimizer(std::move(minimizer)),
      mode(mode)
{}

ProfileResult Profiler::profile(
    std::shared_ptr<ILikelihood> base,
    const ProfileRequest& pr
) const {
    if (pr.free_params.empty()) {
        return direct_profile_no_free(base, pr);
    }

    if (mode == ProfilerMode::LAPLACE_NUISANCE) {
        try {
            return profile_laplace_nuisance(base, pr);
        } catch (const std::exception& e) {
            std::cout << "[LAPLACE PROFILE] failed, falling back to Minuit: "
                      << e.what() << std::endl;
            return profile_minuit(base, pr);
        }
    }

    return profile_minuit(base, pr);
}

ProfileResult Profiler::profile_laplace_nuisance(
    std::shared_ptr<ILikelihood> base,
    const ProfileRequest& pr
) const {
    auto like = std::dynamic_pointer_cast<IProfileableLikelihood>(base);

    if (!like) {
        throw std::runtime_error("Likelihood does not implement IProfileableLikelihood");
    }

    const std::size_t p_dim = like->p_dimension();
    const std::size_t eta_dim = like->eta_dimension();

    if (p_dim + eta_dim != like->dim()) {
        throw std::runtime_error("Inconsistent p/eta dimensions in likelihood");
    }

    // This fast profiler is meant for 2D slice contours: p is fixed and eta is profiled.
    // If the requested strategy also frees a fit parameter, fall back to full Minuit.
    for (std::size_t i : pr.free_params) {
        if (i < p_dim) {
            throw std::runtime_error(
                "Laplace nuisance profiler only supports fixed fit parameters; "
                "free fit parameters require Minuit fallback"
            );
        }
    }

    std::vector<double> p = like->central_p();

    for (const auto& [idx, val] : pr.fixed_params) {
        if (idx < p_dim) {
            p[idx] = val;
        }
    }

    LaplaceProfileOptions opts;
    // These defaults are intentionally conservative:
    // - auto-detect at most a few non-stationary nuisance directions;
    // - correct them with cheap damped Newton steps, not an inner Minuit.
    opts.stationarity_threshold = 5e-2;
    opts.max_refined_eta = 4;
    opts.max_refinement_iters = 2;
    opts.max_newton_step_in_sigma = 1.0;
    opts.use_direct_nll_for_final_value = false;

    const LaplaceProfileComputation comp =
        laplace_profile_eta_refined(*like, p, opts);

    return full_theta_result_from_split(
        p,
        comp.eta_hat,
        comp.nll_hat,
        comp.ok
    );
}


ProfileResult Profiler::profile_minuit(std::shared_ptr<ILikelihood> base, const ProfileRequest& pr) const {
    if (pr.fixed_params.size() + pr.free_params.size() != base->dim()) {
        LOG_ERROR("InvalidArgument", "Dimension mismatch in profiler. Fit dimension is", base->dim(),
                  ", found", pr.fixed_params.size(), "fixed params and", pr.free_params.size(), "free.");
    }

    auto unzipped = unzip(pr.fixed_params);
    std::vector<std::size_t> fixed_idx = unzipped.ids;
    std::vector<double> fixed_vals = unzipped.vals;

    auto f = fit_app::LambdaObjectiveFunction(
        [base](std::vector<double> theta) { return base->nll(theta); },
        0.5
    );

    auto make_defs_from = [&](const std::vector<double>& start) {
        auto defs = base->get_param_defs();
        if (!start.empty()) {
            if (start.size() != defs.size()) {
                LOG_ERROR("InvalidArgument", "ProfileRequest::start has wrong dimension.");
            }
            for (std::size_t i = 0; i < defs.size(); ++i) {
                defs[i].value = start[i];
            }
        }
        for (std::size_t k = 0; k < fixed_idx.size(); ++k) {
            defs[fixed_idx[k]].value = fixed_vals[k];
        }
        return defs;
    };

    auto central_seed = [&]() {
        auto defs = base->get_param_defs();
        std::vector<double> out(defs.size());
        for (std::size_t i = 0; i < defs.size(); ++i) out[i] = defs[i].value;
        for (std::size_t k = 0; k < fixed_idx.size(); ++k) out[fixed_idx[k]] = fixed_vals[k];
        return out;
    }();

    fit_app::FitOptions opt;
    opt.run_hesse = false;
    opt.verbose = false;
    opt.strategy = 2;
    opt.max_fcn = 30000;
    opt.tolerance = 0.2;

    fit_app::BackendFitResult best_raw;
    bool best_raw_set = false;

    auto try_one = [&](const std::vector<double>& seed,
                       unsigned strategy,
                       unsigned max_fcn,
                       double tolerance) -> fit_app::BackendFitResult {
        auto defs = make_defs_from(seed);
        fit_app::FitOptions local_opt = opt;
        local_opt.strategy = strategy;
        local_opt.max_fcn = max_fcn;
        local_opt.tolerance = tolerance;
        return minimizer->minimize_with_fixed(f, defs, local_opt, fixed_idx, fixed_vals);
    };

    // 1) warm start courant
    fit_app::BackendFitResult r1 = try_one(pr.start.empty() ? central_seed : pr.start, 2, 30000, 0.2);
    best_raw = r1;
    best_raw_set = true;

    if (!acceptable_for_profile(r1)) {
        // 2) restart depuis le seed global/central
        fit_app::BackendFitResult r2 = try_one(central_seed, 2, 60000, 0.2);
        if ((!best_raw_set) || (std::isfinite(r2.diagnostics.fmin) && r2.diagnostics.fmin < best_raw.diagnostics.fmin)) {
            best_raw = r2;
            best_raw_set = true;
        }

        // 3) dernier essai plus permissif
        if (!acceptable_for_profile(best_raw)) {
            fit_app::BackendFitResult r3 = try_one(central_seed, 1, 100000, 0.5);
            if (std::isfinite(r3.diagnostics.fmin) && r3.diagnostics.fmin < best_raw.diagnostics.fmin) {
                best_raw = r3;
            }
        }
    }

    const bool accepted = acceptable_for_profile(best_raw);

    if (!accepted) {
        std::cout << "\n=== [PROFILE FAIL] ===\n";
        std::cout << "fixed_idx = [ ";
        for (auto i : fixed_idx) std::cout << i << " ";
        std::cout << "]\n";

        std::cout << "fixed_vals = [ ";
        for (auto v : fixed_vals) std::cout << std::setprecision(17) << v << " ";
        std::cout << "]\n";

        std::cout << "fmin               = " << best_raw.diagnostics.fmin << "\n";
        std::cout << "edm                = " << best_raw.diagnostics.edm << "\n";
        std::cout << "nfcn               = " << best_raw.diagnostics.nfcn << "\n";
        std::cout << "ok                 = " << best_raw.diagnostics.ok << "\n";
        std::cout << "has_valid_params   = " << best_raw.diagnostics.has_valid_parameters << "\n";
        std::cout << "has_valid_covar    = " << best_raw.diagnostics.has_valid_covar << "\n";
        std::cout << "has_accurate_covar = " << best_raw.diagnostics.has_accurate_covar << "\n";
        std::cout << "has_posdef_covar   = " << best_raw.diagnostics.has_posdef_covar << "\n";
        std::cout << "made_posdef        = " << best_raw.diagnostics.made_posdef << "\n";
        std::cout << "hesse_failed       = " << best_raw.diagnostics.hesse_failed << "\n";
    }

    ProfileResult pres;
    pres.nll_hat = best_raw.diagnostics.fmin;
    for (std::size_t i : pr.free_params) {
        pres.theta_hat[i] = best_raw.values[i];
    }
    pres.converged = accepted;

    return pres;
}
