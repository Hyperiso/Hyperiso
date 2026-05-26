#include "Include.h"
#include "GradientHelper.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

static bool approx(double a, double b, double eps = 1e-5) {
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
    std::cout << "== GradientHelper INTEGRATION ==\n";

    QuadraticProfileableLikelihood like;
    const std::vector<double> p{1.0};

    const double r0 = p[0] - like.y;
    const double norm2 = 2.0 * 2.0 + (-3.0) * (-3.0);
    const double denom = 1.0 + norm2;

    const std::vector<double> expected_eta{
        -r0 * 2.0 / denom,
        -r0 * (-3.0) / denom
    };
    const double expected_nll = 0.5 * r0 * r0 / denom;

    {
        LaplaceProfileComputation comp = laplace_profile_eta(like, p);
        assert(comp.ok);
        assert(comp.eta_hat.size() == 2);
        assert(approx(comp.eta_hat[0], expected_eta[0], 2e-4));
        assert(approx(comp.eta_hat[1], expected_eta[1], 2e-4));
        assert(approx(comp.nll_hat, expected_nll, 2e-4));
    }

    {
        LaplaceProfileOptions opts;
        opts.stationarity_threshold = 1e9; // pure Laplace, no Newton refinement requested.
        opts.max_refinement_iters = 3;
        opts.use_direct_nll_for_final_value = true;

        LaplaceProfileComputation comp = laplace_profile_eta_refined(like, p, opts);
        assert(comp.ok);
        assert(approx(comp.eta_hat[0], expected_eta[0], 2e-4));
        assert(approx(comp.eta_hat[1], expected_eta[1], 2e-4));
        assert(approx(comp.nll_hat, like.nll_from_split(p, comp.eta_hat), 2e-4));
    }

    {
        auto ranked = rank_eta_stationarity(like, p, std::vector<double>{0.0, 0.0});
        assert(!ranked.empty());
        assert(ranked.front().first > 0.0);

        LaplaceProfileOptions opts;
        opts.stationarity_threshold = 0.0;
        opts.max_refined_eta = 1;
        auto selected = select_nonstationary_eta_indices(like, p, std::vector<double>{0.0, 0.0}, opts);
        assert(selected.size() == 1);
        assert(selected[0] == 0 || selected[0] == 1);
    }

    std::cout << "\nGradientHelper INTEGRATION tests passed!\n";
    return 0;
}
