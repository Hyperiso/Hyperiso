#ifndef __BASELIKELIHOOD_H__
#define __BASELIKELIHOOD_H__

#include "ILikelihood.h"
#include "JointDistribution.h"
#include "Math.h"
#include <memory>
#include <vector>

using ModelFn = std::function<Vector(const Vector& p, const Vector& eta)>;

struct LikelihoodContext {
    std::unique_ptr<JointDistribution> nuisance_dist;
    std::unique_ptr<JointDistribution> exp_obs_dist;
    Vector exp_obs_values;
    std::vector<fit_app::ParameterDefinition> nuis_defs;
    std::vector<fit_app::ParameterDefinition> fp_defs;
};

class BaseLikelihood : public ILikelihood {
public:
    BaseLikelihood(const ModelFn& model, LikelihoodContext ctx, size_t p_dim);
    double nll(const Vector& theta) const override;

    std::vector<fit_app::ParameterDefinition> get_param_defs() const override {
        auto p_defs = ctx.fp_defs;
        p_defs.insert(p_defs.end(), ctx.nuis_defs.begin(), ctx.nuis_defs.end());
        return p_defs;
    }

protected:
    LikelihoodContext ctx;
    ModelFn model;
    std::size_t p_dim;
};

#endif // __BASELIKELIHOOD_H__
