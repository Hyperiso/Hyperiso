#ifndef __WITHPROFILING_H__
#define __WITHPROFILING_H__

#include "ILikelihood.h"
#include "Math.h"
#include "Indexing.h"

struct ProfileRequest {
    std::vector<std::size_t> free_params;
    std::map<std::size_t, double> fixed_params;
    std::vector<double> start;
};

struct ProfileResult {
    double nll_hat;
    std::map<std::size_t, double> theta_hat;
    bool converged;
};

class Profiler {
public:
    Profiler(std::shared_ptr<fit_app::IFitBackend> minimizer);
    ProfileResult profile(std::shared_ptr<ILikelihood> base, const ProfileRequest& pr) const;

private:
    std::shared_ptr<fit_app::IFitBackend> minimizer;
};

#endif // __WITHPROFILING_H__
