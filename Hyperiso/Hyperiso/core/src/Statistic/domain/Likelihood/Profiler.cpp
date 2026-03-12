#include "Profiler.h"

Profiler::Profiler(std::shared_ptr<fit_app::IFitBackend> minimizer) 
    : minimizer(std::move(minimizer)) {}

ProfileResult Profiler::profile(std::shared_ptr<ILikelihood> base, const ProfileRequest& pr) const {
    if (pr.fixed_params.size() + pr.free_params.size() != base->dim()) 
        LOG_ERROR("InvalidArgument", "Dimension mismatch in profiler.");

    auto unzipped = unzip(pr.fixed_params);
    std::vector<std::size_t> fixed_idx = unzipped.ids;
    Vector fixed_vals = unzipped.vals;

    auto f = fit_app::LambdaObjectiveFunction(
        [base] (Vector theta) { return base->nll(theta); },
        0.5
    );
    
    fit_app::FitOptions opt;
    opt.run_hesse = false;
    opt.verbose = true; // TODO : for now, remove later

    auto res = minimizer->minimize_with_fixed(f, base->get_param_defs(), opt, fixed_idx, fixed_vals);

    if (!res.diagnostics.ok)
        LOG_WARN("Profiling failed to converge, check results.");

    ProfileResult pres;
    pres.nll_hat = res.diagnostics.fmin;
    pres.theta_hat = zip(pr.free_params, res.values);
    pres.converged = res.diagnostics.ok;

    return pres;
}
