#include "ProfiledLikelihood2D.h"

ProfiledLikelihood2D::ProfiledLikelihood2D(std::shared_ptr<ILikelihood> base, std::shared_ptr<Profiler> profiler, std::shared_ptr<IProfilingStrategy> profiling_strategy)
: base(std::move(base)), profiler(std::move(profiler)), profiling_strategy(std::move(profiling_strategy))
{}

double ProfiledLikelihood2D::profiled_nll(double px, double py) {
    ProfileRequest pr = profiling_strategy->build_request(px, py, this->last);
    ProfileResult res = profiler->profile(base, pr);

    if (res.converged)
        this->last = res.theta_hat;
    else
        LOG_WARN("Profiling failed to converge. Check results.");

    return res.nll_hat;
}
