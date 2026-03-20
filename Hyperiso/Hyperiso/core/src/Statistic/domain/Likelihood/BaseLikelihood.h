#ifndef __BASELIKELIHOOD_H__
#define __BASELIKELIHOOD_H__

#include "ILikelihood.h"
#include "JointDistribution.h"
#include "Math.h"
#include <memory>
#include <vector>

using ModelFn = std::function<std::vector<double>(const std::vector<double>& p, const std::vector<double>& eta)>;

struct LikelihoodContext {
    std::unique_ptr<JointDistribution> nuisance_dist;
    std::unique_ptr<JointDistribution> exp_obs_dist;
    std::vector<double>exp_obs_values;
    std::vector<fit_app::ParameterDefinition> nuis_defs;
    std::vector<fit_app::ParameterDefinition> fp_defs;
};

class BaseLikelihood : public ILikelihood {
public:
    BaseLikelihood(const ModelFn& model, std::shared_ptr<LikelihoodContext> ctx, size_t p_dim);
    double nll(const std::vector<double>& theta) const override;

    std::vector<fit_app::ParameterDefinition> get_param_defs() const override {
        auto p_defs = ctx->fp_defs;
        p_defs.insert(p_defs.end(), ctx->nuis_defs.begin(), ctx->nuis_defs.end());
        return p_defs;
    }

    std::size_t dim() const override;

protected:
    std::shared_ptr<LikelihoodContext> ctx;
    ModelFn model;
    std::size_t p_dim;
};

#endif // __BASELIKELIHOOD_H__
