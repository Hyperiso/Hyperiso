#ifndef __PROFILEDLIKELIHOOD2D_H__
#define __PROFILEDLIKELIHOOD2D_H__

#include "ILikelihood.h"
#include "Profiler.h"
#include "IProfilingStrategy.h"

class ProfiledLikelihood2D {
public:
    ProfiledLikelihood2D() = default;
    ProfiledLikelihood2D(
        std::shared_ptr<ILikelihood> base, 
        std::shared_ptr<Profiler> profiler, 
        std::shared_ptr<IProfilingStrategy> profiling_strategy
    );
    
    double profiled_nll(double px, double py);
    std::array<fit_app::ParameterDefinition, 2> get_param_defs() const;

private:
    std::shared_ptr<ILikelihood> base;
    std::shared_ptr<Profiler> profiler;
    std::shared_ptr<IProfilingStrategy> profiling_strategy;
};

#endif // __PROFILEDLIKELIHOOD2D_H__
