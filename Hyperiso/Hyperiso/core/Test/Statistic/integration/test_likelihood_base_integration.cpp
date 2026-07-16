#include "Include.h"
#include "BaseLikelihood.h"
#include "ICopula.h"
#include "IMarginalDistribution.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <gsl/gsl_cdf.h>

static bool approx(double a, double b, double eps = 1e-8) {
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

class StandardQuadraticMarginal : public IMarginalDistribution {
public:
    std::vector<double> rvs(std::size_t n) override { return std::vector<double>(n, 0.0); }
    double logpdf(double x) override { return -0.5 * x * x; }
    PDFDiff f_df_ddf(double x) override {
        const double f = std::exp(-0.5 * x * x);
        return {f, -x * f, (x * x - 1.0) * f};
    }
    double cdf(double x) override { return gsl_cdf_ugaussian_P(x); }
    double ppf(double p) override { return gsl_cdf_ugaussian_Pinv(std::clamp(p, 1e-13, 1.0 - 1e-13)); }
    double mean() override { return 0.0; }
    double std() override { return 1.0; }
};

class IndependentCopula : public ICopula {
public:
    explicit IndependentCopula(std::size_t d) : d_(d) {}
    std::vector<std::vector<double>> sample_u(std::size_t n) override {
        std::vector<std::vector<double>> out(n, std::vector<double>(d_, 0.5));
        return out;
    }
    std::vector<double> sample_u() override { return std::vector<double>(d_, 0.5); }
    double log_density(std::vector<double>) override { return 0.0; }
    RealMatrix dlog_density(std::vector<double>) override { return RealMatrix(d_, 1); }
    RealMatrix ddlog_density(std::vector<double>) override { return RealMatrix(d_, d_); }
    LogDensityDiff log_c_dc_ddc(std::vector<double>) override {
        return {0.0, RealMatrix(d_, 1), RealMatrix(d_, d_)};
    }
private:
    std::size_t d_;
};

static std::unique_ptr<JointDistribution> make_standard_joint(std::size_t d) {
    std::vector<std::unique_ptr<IMarginalDistribution>> marginals;
    for (std::size_t i = 0; i < d; ++i) {
        marginals.push_back(std::make_unique<StandardQuadraticMarginal>());
    }
    return std::make_unique<JointDistribution>(std::move(marginals), std::make_unique<IndependentCopula>(d));
}

static std::shared_ptr<LikelihoodContext> make_context() {
    auto ctx = std::make_shared<LikelihoodContext>();
    ctx->exp_obs_dist = make_standard_joint(2);
    ctx->nuisance_dist = make_standard_joint(2);
    ctx->exp_obs_values = {1.0, -2.0};
    ctx->fp_defs = {make_def("alpha", 0.5), make_def("beta", -1.0)};
    ctx->nuis_defs = {make_def("n0", 0.0), make_def("n1", 0.0)};
    return ctx;
}

int main() {
    std::cout << "== BaseLikelihood INTEGRATION ==\n";

    auto ctx = make_context();

    // Linear observable model with two fit parameters and two nuisance parameters.
    ModelFn model = [](const std::vector<double>& p,
                       const std::vector<double>& eta) {
        return std::vector<double>{
            p[0] + 0.5 * p[1] + eta[0],
            2.0 * p[0] - p[1] + eta[1]
        };
    };

    BaseLikelihood like(model, ctx, 2);

    const std::vector<double> p_best{0.0, -2.0};
    const std::vector<double> eta_best{0.0, 0.0};
    const std::vector<double> p_shift{1.0, -2.0};
    const std::vector<double> eta_shift{0.5, -0.25};

    const double nll_best = like.nll_from_split(p_best, eta_best);
    const double nll_shift = like.nll_from_split(p_shift, eta_shift);

    assert(std::isfinite(nll_best));
    assert(std::isfinite(nll_shift));
    assert(nll_best < nll_shift);

    auto r_best = like.residuals(p_best, eta_best);
    assert(r_best.size() == 2);
    assert(approx(r_best[0], -2.0, 1e-12));
    assert(approx(r_best[1], 4.0, 1e-12));

    RealMatrix Wobs = like.observable_curvature(r_best);
    RealMatrix Weta = like.nuisance_curvature(eta_best);
    assert(Wobs.rows() == 2 && Wobs.cols() == 2);
    assert(Weta.rows() == 2 && Weta.cols() == 2);
    assert(Wobs.at(0, 0) > 0.0);
    assert(Wobs.at(1, 1) > 0.0);
    assert(Weta.at(0, 0) > 0.0);
    assert(Weta.at(1, 1) > 0.0);

    auto defs = like.get_param_defs();
    assert(defs.size() == like.dim());
    assert(defs[0].name == "alpha");
    assert(defs[1].name == "beta");
    assert(defs[2].name == "n0");
    assert(defs[3].name == "n1");

    std::cout << "\nBaseLikelihood INTEGRATION tests passed!\n";
    return 0;
}
