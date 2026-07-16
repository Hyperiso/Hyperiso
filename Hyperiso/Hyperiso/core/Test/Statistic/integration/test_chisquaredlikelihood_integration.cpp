#include "Include.h"
#include "ChiSquaredLikelihood.h"
#include "Profiler.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

static bool approx(double a, double b, double eps = 1e-10) {
    return std::fabs(a - b) <= eps;
}

static fit_app::ParameterDefinition make_def(
    const std::string& name,
    double value,
    double step_hint = 1.0
) {
    fit_app::ParameterDefinition d;
    d.name = name;
    d.value = value;
    d.step_hint = step_hint;
    return d;
}

static std::shared_ptr<LikelihoodContext> make_context() {
    auto ctx = std::make_shared<LikelihoodContext>();
    ctx->exp_obs_values = {5.0, -1.0};
    ctx->fp_defs = {make_def("a", 1.0), make_def("b", 2.0)};
    ctx->nuis_defs = {};
    return ctx;
}

int main() {
    std::cout << "== ChiSquaredLikelihood INTEGRATION ==\n";

    auto ctx = make_context();

    RealMatrix Cinv(2, 2);
    Cinv.at(0, 0) = 4.0;
    Cinv.at(0, 1) = 1.0;
    Cinv.at(1, 0) = 1.0;
    Cinv.at(1, 1) = 2.0;

    ModelFn model = [](const std::vector<double>& p,
                       const std::vector<double>& eta) {
        assert(eta.empty());
        return std::vector<double>{2.0 * p[0] + p[1], p[0] - p[1]};
    };

    auto like = std::make_shared<ChiSquaredLikelihood>(model, ctx, 2, Cinv);

    {
        const std::vector<double> theta{2.0, 1.0};
        const std::vector<double> pred{5.0, 1.0};
        const std::vector<double> r{0.0, 2.0};
        const double q = r[0] * (4.0 * r[0] + 1.0 * r[1])
                       + r[1] * (1.0 * r[0] + 2.0 * r[1]);

        assert(approx(like->nll(theta), 0.5 * q));
    }

    {
        Profiler profiler(nullptr, ProfilerMode::MINUIT);

        ProfileRequest pr;
        pr.free_params = {};
        pr.fixed_params = {{0, 2.0}, {1, 1.0}};
        pr.start = {}; // direct path builds zero vector, then applies fixed parameters.

        ProfileResult r = profiler.profile(like, pr);
        assert(r.converged);
        assert(approx(r.nll_hat, like->nll(std::vector<double>{2.0, 1.0})));
        assert(r.theta_hat.size() == 2);
        assert(approx(r.theta_hat.at(0), 2.0));
        assert(approx(r.theta_hat.at(1), 1.0));
    }

    {
        RealMatrix W = like->observable_curvature(std::vector<double>{1.0, 2.0});
        assert(W.rows() == 2 && W.cols() == 2);
        assert(approx(W.at(0, 0), 4.0));
        assert(approx(W.at(0, 1), 1.0));
        assert(approx(W.at(1, 0), 1.0));
        assert(approx(W.at(1, 1), 2.0));
    }

    std::cout << "\nChiSquaredLikelihood INTEGRATION tests passed!\n";
    return 0;
}
