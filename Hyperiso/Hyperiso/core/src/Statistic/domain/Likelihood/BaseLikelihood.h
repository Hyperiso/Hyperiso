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

    void enable_debug_trace(std::size_t max_evals = 25) {
        debug_trace_enabled_ = true;
        debug_trace_max_evals_ = max_evals;
        debug_eval_count_ = 0;
        debug_have_ref_theta_ = false;
        debug_have_ref_res_ = false;
        debug_ref_theta_.clear();
        debug_ref_res_.clear();
    }

    void disable_debug_trace() {
        debug_trace_enabled_ = false;
    }
protected:
    std::shared_ptr<LikelihoodContext> ctx;
    ModelFn model;
    std::size_t p_dim;

    mutable bool debug_trace_enabled_ = false;
    mutable std::size_t debug_trace_max_evals_ = 0;
    mutable std::size_t debug_eval_count_ = 0;

    mutable bool debug_have_ref_theta_ = false;
    mutable bool debug_have_ref_res_ = false;
    mutable std::vector<double> debug_ref_theta_;
    mutable std::vector<double> debug_ref_res_;
    
private:
    static double max_abs_diff(const std::vector<double>& a,
                               const std::vector<double>& b) {
        if (a.size() != b.size()) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        double out = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i) {
            out = std::max(out, std::abs(a[i] - b[i]));
        }
        return out;
    }

    void maybe_log_debug_eval(const std::vector<double>& theta,
                              const std::vector<double>& p,
                              const std::vector<double>& eta,
                              const std::vector<double>& res,
                              double ell_obs,
                              double ell_nuis,
                              double nll_value) const;
};

#endif // __BASELIKELIHOOD_H__
