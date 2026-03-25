#include "Profiler.h"

Profiler::Profiler(std::shared_ptr<fit_app::IFitBackend> minimizer) 
    : minimizer(std::move(minimizer)) {}

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

//TODO : test
ProfileResult Profiler::profile(std::shared_ptr<ILikelihood> base, const ProfileRequest& pr) const {
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

// //TODO : use pr.start because if not "Profiling failed to converge, check results"
// ProfileResult Profiler::profile(std::shared_ptr<ILikelihood> base, const ProfileRequest& pr) const {
//     if (pr.fixed_params.size() + pr.free_params.size() != base->dim())
//         LOG_ERROR("InvalidArgument", "Dimension mismatch in profiler. Fit dimension is", base->dim(),
//                   ", found", pr.fixed_params.size(), "fixed params and", pr.free_params.size(), "free.");

//     auto unzipped = unzip(pr.fixed_params);
//     std::vector<std::size_t> fixed_idx = unzipped.ids;
//     std::vector<double> fixed_vals = unzipped.vals;

//     auto base_nll = [base](std::vector<double> theta) { return base->nll(theta); };

//     auto f = fit_app::LambdaObjectiveFunction(
//         base_nll,
//         0.5
//     );

//     fit_app::FitOptions opt;
//     opt.run_hesse = false;
//     opt.verbose = false;

//     auto defs = base->get_param_defs();
//     if (!pr.start.empty()) {
//         if (pr.start.size() != defs.size()) {
//             LOG_ERROR("InvalidArgument", "ProfileRequest::start has wrong dimension.");
//         }
//         for (std::size_t i = 0; i < defs.size(); ++i) {
//             defs[i].value = pr.start[i];
//         }
//     }

//     std::cout << "pr.start = (";
//     for (size_t i = 0; i < pr.start.size(); i++) {
//         std::cout << pr.start[i] << ", ";
//     }
//     std::cout << ")" << std::endl;

//     auto res = minimizer->minimize_with_fixed(f, defs, opt, fixed_idx, fixed_vals);

//     // if (!res.diagnostics.ok)
//     //     LOG_WARN("Profiling failed to converge, check results.");
//     const bool acceptable_for_profile =
//         res.diagnostics.has_valid_parameters &&
//         std::isfinite(res.diagnostics.fmin) &&
//         (
//             res.diagnostics.ok ||
//             (!res.diagnostics.hesse_failed && res.diagnostics.edm < 1e-4)
//         );

//     if (!acceptable_for_profile) {
//         std::cout << "\n=== [PROFILE FAIL] ===\n";
//         std::cout << "fixed_idx = [ ";
//         for (auto i : fixed_idx) std::cout << i << " ";
//         std::cout << "]\n";

//         std::cout << "fixed_vals = [ ";
//         for (auto v : fixed_vals) std::cout << std::setprecision(17) << v << " ";
//         std::cout << "]\n";

//         std::cout << "fmin               = " << res.diagnostics.fmin << "\n";
//         std::cout << "edm                = " << res.diagnostics.edm << "\n";
//         std::cout << "nfcn               = " << res.diagnostics.nfcn << "\n";
//         std::cout << "ok                 = " << res.diagnostics.ok << "\n";
//         std::cout << "has_valid_params   = " << res.diagnostics.has_valid_parameters << "\n";
//         std::cout << "has_valid_covar    = " << res.diagnostics.has_valid_covar << "\n";
//         std::cout << "has_accurate_covar = " << res.diagnostics.has_accurate_covar << "\n";
//         std::cout << "has_posdef_covar   = " << res.diagnostics.has_posdef_covar << "\n";
//         std::cout << "made_posdef        = " << res.diagnostics.made_posdef << "\n";
//         std::cout << "hesse_failed       = " << res.diagnostics.hesse_failed << "\n";
//     }

//     ProfileResult pres;
//     pres.nll_hat = res.diagnostics.fmin;
//     for (std::size_t i : pr.free_params)
//         pres.theta_hat[i] = res.values[i];
//     pres.converged = acceptable_for_profile;

//     return pres;
//     // ProfileResult pres;
//     // pres.nll_hat = res.diagnostics.fmin;
//     // for (std::size_t i : pr.free_params)
//     //     pres.theta_hat[i] = res.values[i];
//     // pres.converged = res.diagnostics.ok;

//     // return pres;
// }

// ProfileResult Profiler::profile(std::shared_ptr<ILikelihood> base, const ProfileRequest& pr) const {
//     if (pr.fixed_params.size() + pr.free_params.size() != base->dim()) 
//         LOG_ERROR("InvalidArgument", "Dimension mismatch in profiler. Fit dimension is", base->dim(), ", found", pr.fixed_params.size(), "fixed params and", pr.free_params.size(), "free.");

//     auto unzipped = unzip(pr.fixed_params);
//     std::vector<std::size_t> fixed_idx = unzipped.ids;
//     std::vector<double>fixed_vals = unzipped.vals;

//     auto f = fit_app::LambdaObjectiveFunction(
//         [base] (std::vector<double>theta) { return base->nll(theta); },
//         0.5
//     );
    
//     fit_app::FitOptions opt;
//     opt.run_hesse = false;
//     opt.verbose = false; // TODO : for now, remove later

//     auto res = minimizer->minimize_with_fixed(f, base->get_param_defs(), opt, fixed_idx, fixed_vals);

//     if (!res.diagnostics.ok)
//         LOG_WARN("Profiling failed to converge, check results.");

//     ProfileResult pres;
//     pres.nll_hat = res.diagnostics.fmin;
//     for (std::size_t i : pr.free_params) 
//         pres.theta_hat[i] = res.values[i];
//     pres.converged = res.diagnostics.ok;

//     return pres;
// }
