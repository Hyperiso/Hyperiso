#include "Include.h"
#include "Profiler.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

static bool approx(double a, double b, double eps = 1e-4) {
    return std::fabs(a - b) <= eps;
}

static fit_app::ParameterDefinition make_def(const std::string& name,
                                             double value,
                                             double step_hint = 1.0) {
    fit_app::ParameterDefinition d;
    d.name = name;
    d.value = value;
    d.step_hint = step_hint;
    return d;
}

class QuadraticProfileableLikelihood : public IProfileableLikelihood {
public:
    double y = 1.5;
    std::vector<double> a = {2.0, -3.0};
    std::vector<fit_app::ParameterDefinition> defs;

    QuadraticProfileableLikelihood() {
        defs = {
            make_def("p0", 1.0, 0.1),
            make_def("eta0", 0.0, 1.0),
            make_def("eta1", 0.0, 1.0)
        };
    }

    double nll(const std::vector<double>& theta) const override {
        return nll_from_split({theta[0]}, {theta[1], theta[2]});
    }

    std::vector<fit_app::ParameterDefinition> get_param_defs() const override { return defs; }
    std::size_t dim() const override { return 3; }
    std::size_t p_dimension() const override { return 1; }
    std::size_t eta_dimension() const override { return 2; }
    std::vector<double> central_p() const override { return {defs[0].value}; }
    std::vector<double> central_eta() const override { return {defs[1].value, defs[2].value}; }

    std::vector<double> predict(const std::vector<double>& p,
                                const std::vector<double>& eta) const override {
        return {p[0] + a[0] * eta[0] + a[1] * eta[1]};
    }

    std::vector<double> residuals(const std::vector<double>& p,
                                  const std::vector<double>& eta) const override {
        return {predict(p, eta)[0] - y};
    }

    double nll_from_split(const std::vector<double>& p,
                          const std::vector<double>& eta) const override {
        const double r = residuals(p, eta)[0];
        return 0.5 * r * r + 0.5 * (eta[0] * eta[0] + eta[1] * eta[1]);
    }

    RealMatrix observable_curvature(const std::vector<double>&) const override {
        RealMatrix W(1, 1);
        W.at(0, 0) = 1.0;
        return W;
    }

    RealMatrix nuisance_curvature(const std::vector<double>&) const override {
        RealMatrix W(2, 2);
        W.at(0, 0) = 1.0;
        W.at(0, 1) = 0.0;
        W.at(1, 0) = 0.0;
        W.at(1, 1) = 1.0;
        return W;
    }
};

int main() {
    std::cout << "== Profiler INTEGRATION ==\n";

    auto like = std::make_shared<QuadraticProfileableLikelihood>();
    Profiler profiler(nullptr, ProfilerMode::LAPLACE_NUISANCE);

    ProfileRequest pr;
    pr.fixed_params = {{0, 1.0}};     // fixed fitted parameter p0
    pr.free_params = {1, 2};          // profile nuisance parameters eta0, eta1
    pr.start = {1.0, 0.0, 0.0};

    ProfileResult result = profiler.profile(like, pr);

    const double r0 = 1.0 - like->y;
    const double norm2 = 2.0 * 2.0 + (-3.0) * (-3.0);
    const double denom = 1.0 + norm2;
    const double eta0_hat = -r0 * 2.0 / denom;
    const double eta1_hat = -r0 * (-3.0) / denom;
    const double nll_hat = 0.5 * r0 * r0 / denom;

    assert(result.converged);
    assert(result.theta_hat.size() == 3);
    assert(approx(result.theta_hat.at(0), 1.0, 1e-12));
    assert(approx(result.theta_hat.at(1), eta0_hat, 3e-4));
    assert(approx(result.theta_hat.at(2), eta1_hat, 3e-4));
    assert(approx(result.nll_hat, nll_hat, 3e-4));

    const double direct = like->nll({
        result.theta_hat.at(0),
        result.theta_hat.at(1),
        result.theta_hat.at(2)
    });
    assert(std::isfinite(direct));
    assert(approx(direct, nll_hat, 3e-4));

    std::cout << "\nProfiler INTEGRATION tests passed!\n";
    return 0;
}
