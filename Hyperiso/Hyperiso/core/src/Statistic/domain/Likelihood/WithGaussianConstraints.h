#ifndef __WITHGAUSSIANCONSTRAINTS_H__
#define __WITHGAUSSIANCONSTRAINTS_H__

#include "ILikelihood.h"
#include "JointDistribution.h"

class WithGaussianConstraints : public ILikelihood {
public:
    WithGaussianConstraints(std::shared_ptr<ILikelihood> base, std::shared_ptr<JointDistribution> constraints_dist, std::vector<std::size_t> constrained_params);
    double nll(const Vector& theta) const override;
    std::vector<fit_app::ParameterDefinition> get_param_defs() const override;
    std::size_t dim() const override;

private:
    std::shared_ptr<ILikelihood> base;
    std::shared_ptr<JointDistribution> constraints_dist;
    std::vector<std::size_t> constrained_params;
};

#endif // __WITHGAUSSIANCONSTRAINTS_H__
