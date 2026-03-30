#include "ProfiledLikelihood2D.h"

ProfiledLikelihood2D::ProfiledLikelihood2D(
    std::shared_ptr<ILikelihood> base,
    std::shared_ptr<Profiler> profiler,
    std::shared_ptr<IProfilingStrategy> profiling_strategy
)
: base(std::move(base)),
  profiler(std::move(profiler)),
  profiling_strategy(std::move(profiling_strategy))
{
    this->last = this->profiling_strategy->init_warm_start();
}

double ProfiledLikelihood2D::profiled_nll(double px, double py) {
    ProfileRequest pr = profiling_strategy->build_request(px, py, this->last);
    ProfileResult res = profiler->profile(base, pr);

    if (res.converged) {
        this->last = res.theta_hat;
        return res.nll_hat;
    }

    if (std::isfinite(res.nll_hat)) {
        LOG_INFO("[PROFILED_NLL] weak profile at (", px, ',', py, "), using finite value", res.nll_hat);
        return res.nll_hat;
    }

    LOG_INFO("[PROFILED_NLL] rejected point (", px, ',', py, "), using fallback penalty");
    return 1e300;
}

std::array<fit_app::ParameterDefinition, 2> ProfiledLikelihood2D::get_param_defs() const {
    auto all_p_defs = this->base->get_param_defs();
    return {
        all_p_defs[this->profiling_strategy->get_x_id()],
        all_p_defs[this->profiling_strategy->get_y_id()]
    };
}