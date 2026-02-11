#ifndef __LIKELIHOOD_H__
#define __LIKELIHOOD_H__

#include <vector>
#include <functional>
#include "JointDistribution.h"
#include "Math.h"
#include "analysis.h"

struct LikelihoodContext {
    std::unique_ptr<JointDistribution> nuisance_dist;
    std::unique_ptr<JointDistribution> exp_obs_dist;
    Vector exp_obs_values;
    Vector nuisance_central_values;
};

class ProfiledLikelihood {
public:
    using ModelFn = std::function<Vector(const Vector& p, const Vector& eta)>;

    explicit ProfiledLikelihood(LikelihoodContext ctx, ModelFn model);

    double nll(const Vector& p, const Vector& eta) const;
    double nll_profiled(const Vector& p) const;

    void set_minimizer_max_iter(std::size_t max_iter);
    void set_minimizer_tolerance(double tol);

    Vector get_eta_central_values() const;

private:
    Vector residual_obs(const Vector& p, const Vector& eta) const {
        Vector pred = model_(p, eta);
        if (pred.size()!=ctx_.exp_obs_dist->dim()) throw std::invalid_argument("Model produced wrong dimension");
        Vector r(pred.size());
        for (std::size_t i=0;i<pred.size();++i) r[i] = pred[i] - ctx_.exp_obs_values[i];
        return r;
    }

    LikelihoodContext ctx_;
    ModelFn model_;
    MinimizationContext min_ctx_;
};

#endif // __LIKELIHOOD_H__