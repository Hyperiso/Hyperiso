#ifndef __WITHPROFILING_H__
#define __WITHPROFILING_H__

#include "ILikelihood.h"
#include "Math.h"
#include "Indexing.h"
#include "GradientHelper.h"

enum class ProfilerMode {
    MINUIT,
    LAPLACE_NUISANCE
};

struct ProfileRequest {
    std::vector<std::size_t> free_params;
    std::map<std::size_t, double> fixed_params;
    std::vector<double> start;
};

struct ProfileResult {
    double nll_hat = 1e300;
    std::map<std::size_t, double> theta_hat;
    bool converged = false;
};

class Profiler {
public:
    Profiler(
        std::shared_ptr<fit_app::IFitBackend> minimizer,
        ProfilerMode mode = ProfilerMode::MINUIT
    );

    ProfileResult profile(
        std::shared_ptr<ILikelihood> base,
        const ProfileRequest& pr
    ) const;

private:
    ProfileResult profile_minuit(
        std::shared_ptr<ILikelihood> base,
        const ProfileRequest& pr
    ) const;

    ProfileResult profile_laplace_nuisance(
        std::shared_ptr<ILikelihood> base,
        const ProfileRequest& pr
    ) const;

    std::shared_ptr<fit_app::IFitBackend> minimizer;
    ProfilerMode mode = ProfilerMode::MINUIT;
};

#endif // __WITHPROFILING_H__
