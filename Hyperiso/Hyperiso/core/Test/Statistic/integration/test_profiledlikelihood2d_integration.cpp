#include "Include.h"
#include "ChiSquaredLikelihood.h"
#include "ProfiledLikelihood2D.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

static bool approx(double a, double b, double eps = 1e-10) {
    return std::fabs(a - b) <= eps;
}

static fit_app::ParameterDefinition make_def(const std::string& name, double value) {
    fit_app::ParameterDefinition d;
    d.name = name;
    d.value = value;
    d.step_hint = 1.0;
    return d;
}

int main() {
    std::cout << "== ProfiledLikelihood2D INTEGRATION ==\n";

    auto ctx = std::make_shared<LikelihoodContext>();
    ctx->exp_obs_values = {0.0, 0.0};
    ctx->fp_defs = {make_def("a", 0.0), make_def("b", 0.0)};
    ctx->nuis_defs = {};

    RealMatrix Cinv(2, 2);
    Cinv.at(0, 0) = 1.0;
    Cinv.at(0, 1) = 0.0;
    Cinv.at(1, 0) = 0.0;
    Cinv.at(1, 1) = 4.0;

    ModelFn model = [](const std::vector<double>& p,
                       const std::vector<double>& eta) {
        assert(eta.empty());
        return std::vector<double>{p[0] - 1.0, p[1] + 2.0};
    };

    auto base = std::make_shared<ChiSquaredLikelihood>(model, ctx, 2, Cinv);

    FitResult fr;
    fr.p_hat = {1.0, -2.0};
    fr.eta_hat = {};
    fr.p_hat_std = {0.1, 0.1};
    fr.p_hat_correlations = RealMatrix(2, 2);
    fr.ell_hat = 0.0;

    auto profiler = std::make_shared<Profiler>(nullptr, ProfilerMode::MINUIT);
    auto strategy = std::make_shared<SliceProfilingStrategy>(0, 1, fr);
    ProfiledLikelihood2D profiled(base, profiler, strategy);

    {
        const double at_min = profiled.profiled_nll(1.0, -2.0);
        assert(approx(at_min, 0.0));
    }

    {
        const double off = profiled.profiled_nll(2.0, -1.5);
        const double expected = base->nll(std::vector<double>{2.0, -1.5});
        assert(approx(off, expected));
        assert(off > 0.0);
    }

    {
        auto defs = profiled.get_param_defs();
        assert(defs[0].name == "a");
        assert(defs[1].name == "b");
    }

    std::cout << "\nProfiledLikelihood2D INTEGRATION tests passed!\n";
    return 0;
}
