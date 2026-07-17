#include "Include.h"
#include "ChiSquaredLikelihood.h"
#include "Fit.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace {

fit_app::ParameterDefinition make_def(const std::string& name, double value) {
    fit_app::ParameterDefinition def;
    def.name = name;
    def.value = value;
    def.step_hint = 0.1;
    def.limits = std::make_pair(-5.0, 5.0);
    return def;
}

} // namespace

int main() {
    std::cout << "== Fit profile-Hessian UNIT ==\n";

    auto ctx = std::make_shared<LikelihoodContext>();
    ctx->exp_obs_values = {0.0, 0.0};
    for (std::size_t i = 0; i < 6; ++i) {
        ctx->fp_defs.push_back(make_def("p" + std::to_string(i), 0.0));
    }

    RealMatrix covariance_inv(2, 2);
    covariance_inv.at(0, 0) = 1.0;
    covariance_inv.at(0, 1) = 0.0;
    covariance_inv.at(1, 0) = 0.0;
    covariance_inv.at(1, 1) = 1.0;

    ModelFn model = [](const std::vector<double>& p,
                       const std::vector<double>& eta) {
        assert(eta.empty());
        assert(p.size() == 6);
        return std::vector<double>{p[0], p[1]};
    };

    auto likelihood = std::make_shared<ChiSquaredLikelihood>(
        model,
        ctx,
        6,
        covariance_inv
    );

    MLFitOptions options;
    options.run_hesse = true;
    options.allow_profile_hessian_fallback = true;
    options.max_fcn = 20000;
    options.tolerance = 0.1;

    MLFitter fitter(likelihood, options);
    FitResult result = fitter.maximum_likelihood_fit(
        {0.4, -0.3, 0.0, 0.0, 0.0, 0.0}
    );

    assert(result.p_hat.size() == 6);
    assert(result.p_hat_std.size() == 6);
    for (double sigma : result.p_hat_std) {
        assert(std::isfinite(sigma));
    }

    std::cout << "Fit profile-Hessian UNIT tests passed!\n";
    return 0;
}
