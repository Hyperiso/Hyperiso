#include "Include.h"
#include "ProfiledLikelihood2D.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

static bool approx(double a, double b, double eps = 1e-12) {
    return std::fabs(a - b) <= eps;
}

static fit_app::ParameterDefinition make_def(const std::string& name, double value) {
    fit_app::ParameterDefinition d;
    d.name = name;
    d.value = value;
    d.step_hint = 1.0;
    return d;
}

class TwoDimLikelihood : public ILikelihood {
public:
    mutable std::vector<double> last_theta;
    mutable int nll_calls = 0;

    double nll(const std::vector<double>& theta) const override {
        nll_calls++;
        last_theta = theta;
        return 0.5 * ((theta[0] - 1.0) * (theta[0] - 1.0)
                    + 4.0 * (theta[1] + 2.0) * (theta[1] + 2.0));
    }

    std::vector<fit_app::ParameterDefinition> get_param_defs() const override {
        return {make_def("x", 1.0), make_def("y", -2.0)};
    }

    std::size_t dim() const override { return 2; }
};

static FitResult make_fit_result() {
    FitResult fr;
    fr.p_hat = {1.0, -2.0};
    fr.eta_hat = {};
    fr.p_hat_std = {0.1, 0.2};
    fr.p_hat_correlations = RealMatrix(2, 2);
    fr.ell_hat = 0.0;
    return fr;
}

int main() {
    std::cout << "== ProfiledLikelihood2D UNIT ==\n";

    {
        auto base = std::make_shared<TwoDimLikelihood>();
        auto profiler = std::make_shared<Profiler>(nullptr, ProfilerMode::MINUIT);
        auto strategy = std::make_shared<SliceProfilingStrategy>(0, 1, make_fit_result());

        ProfiledLikelihood2D profiled(base, profiler, strategy);

        const double out = profiled.profiled_nll(2.0, -1.0);
        const double expected = base->nll(std::vector<double>{2.0, -1.0});

        assert(approx(out, expected));
        assert(base->nll_calls >= 2); // one direct expected call plus the profiled call.
        assert((base->last_theta == std::vector<double>{2.0, -1.0}));
    }

    {
        auto base = std::make_shared<TwoDimLikelihood>();
        auto profiler = std::make_shared<Profiler>(nullptr, ProfilerMode::MINUIT);
        auto strategy = std::make_shared<SliceProfilingStrategy>(0, 1, make_fit_result());

        ProfiledLikelihood2D profiled(base, profiler, strategy);
        auto defs = profiled.get_param_defs();

        assert(defs[0].name == "x");
        assert(defs[1].name == "y");
        assert(approx(defs[0].value, 1.0));
        assert(approx(defs[1].value, -2.0));
    }

    {
        auto base = std::make_shared<TwoDimLikelihood>();
        auto profiler = std::make_shared<Profiler>(nullptr, ProfilerMode::MINUIT);
        auto strategy = std::make_shared<ProjectionProfilingStrategy>(0, 1, make_fit_result());

        ProfiledLikelihood2D profiled(base, profiler, strategy);

        const double v1 = profiled.profiled_nll(0.0, 0.0);
        const double v2 = profiled.profiled_nll(1.0, -2.0);

        assert(std::isfinite(v1));
        assert(std::isfinite(v2));
        assert(v2 < v1);
    }

    std::cout << "\nProfiledLikelihood2D UNIT tests passed!\n";
    return 0;
}
