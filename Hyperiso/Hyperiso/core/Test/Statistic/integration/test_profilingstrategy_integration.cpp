#include "Include.h"
#include "IProfilingStrategy.h"
#include "Profiler.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
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

class QuadraticLikelihood : public ILikelihood {
public:
    std::vector<fit_app::ParameterDefinition> defs = {
        make_def("p0", 0.0),
        make_def("p1", 0.0),
        make_def("eta0", 0.0)
    };

    mutable std::vector<double> last_theta;

    double nll(const std::vector<double>& theta) const override {
        last_theta = theta;
        double out = 0.0;
        for (std::size_t i = 0; i < theta.size(); ++i) {
            out += static_cast<double>(i + 1) * theta[i] * theta[i];
        }
        return out;
    }

    std::vector<fit_app::ParameterDefinition> get_param_defs() const override {
        return defs;
    }

    std::size_t dim() const override { return defs.size(); }
};

int main() {
    std::cout << "== IProfilingStrategy INTEGRATION ==\n";

    FitResult fr;
    fr.p_hat = {1.0, 2.0};
    fr.eta_hat = {0.5};
    fr.p_hat_std = {0.1, 0.2};
    fr.p_hat_correlations = RealMatrix(2, 2);
    fr.ell_hat = 0.0;

    auto like = std::make_shared<QuadraticLikelihood>();
    Profiler profiler(nullptr, ProfilerMode::MINUIT);

    {
        SliceProfilingStrategy s(0, 1, fr);
        auto warm = s.init_warm_start();
        assert(warm.size() == 1);
        assert(approx(warm.at(2), 0.5));

        ProfileRequest pr = s.build_request(3.0, -4.0, warm);
        assert((pr.free_params == std::vector<std::size_t>{2}));
        assert(approx(pr.fixed_params.at(0), 3.0));
        assert(approx(pr.fixed_params.at(1), -4.0));
        assert(approx(pr.start[2], 0.5));
    }

    {
        // Projection with a 2D fit-result and no nuisance dimension fixes both available coordinates.
        FitResult fr2;
        fr2.p_hat = {1.0, 2.0};
        fr2.eta_hat = {};
        fr2.p_hat_std = {0.1, 0.2};
        fr2.p_hat_correlations = RealMatrix(2, 2);

        ProjectionProfilingStrategy s(0, 1, fr2);
        ProfileRequest pr = s.build_request(3.0, -4.0, s.init_warm_start());
        assert(pr.free_params.empty());

        auto like2 = std::make_shared<QuadraticLikelihood>();
        like2->defs.resize(2);
        ProfileResult r = profiler.profile(like2, pr);
        assert(r.converged);
        assert(approx(r.nll_hat, 3.0 * 3.0 + 2.0 * (-4.0) * (-4.0)));
        assert((like2->last_theta == std::vector<double>{3.0, -4.0}));
    }

    std::cout << "\nIProfilingStrategy INTEGRATION tests passed!\n";
    return 0;
}
