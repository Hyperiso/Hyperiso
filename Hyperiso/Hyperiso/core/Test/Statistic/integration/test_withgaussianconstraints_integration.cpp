#include "Include.h"
#include "ChiSquaredLikelihood.h"
#include "WithGaussianConstraints.h"
#include "ICopula.h"
#include "IMarginalDistribution.h"

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

class QuadraticMarginal : public IMarginalDistribution {
public:
    explicit QuadraticMarginal(double mu = 0.0, double sigma = 1.0)
        : mu_(mu), sigma_(sigma)
    {}

    std::vector<double> rvs(std::size_t n) override { return std::vector<double>(n, mu_); }

    double logpdf(double x) override {
        const double z = (x - mu_) / sigma_;
        return -0.5 * z * z;
    }

    PDFDiff f_df_ddf(double x) override {
        const double z = (x - mu_) / sigma_;
        const double s2 = sigma_ * sigma_;
        const double f = std::exp(-0.5 * z * z);
        return {f, -(x - mu_) / s2 * f, (z * z - 1.0) / s2 * f};
    }

    double cdf(double) override { return 0.5; }
    double ppf(double) override { return mu_; }
    double mean() override { return mu_; }
    double std() override { return sigma_; }

private:
    double mu_, sigma_;
};

class IndependentCopula : public ICopula {
public:
    explicit IndependentCopula(std::size_t d) : d_(d) {}
    std::vector<std::vector<double>> sample_u(std::size_t n) override {
        return std::vector<std::vector<double>>(n, std::vector<double>(d_, 0.5));
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

static std::shared_ptr<JointDistribution> make_constraint_dist(double mu, double sigma) {
    std::vector<std::unique_ptr<IMarginalDistribution>> marginals;
    marginals.push_back(std::make_unique<QuadraticMarginal>(mu, sigma));
    return std::make_shared<JointDistribution>(
        std::move(marginals),
        std::make_unique<IndependentCopula>(1)
    );
}

int main() {
    std::cout << "== WithGaussianConstraints INTEGRATION ==\n";

    auto ctx = std::make_shared<LikelihoodContext>();
    ctx->exp_obs_values = {0.0};
    ctx->fp_defs = {make_def("x", 0.0), make_def("y", 0.0)};
    ctx->nuis_defs = {};

    RealMatrix Cinv(1, 1);
    Cinv.at(0, 0) = 2.0;

    ModelFn model = [](const std::vector<double>& p,
                       const std::vector<double>& eta) {
        assert(eta.empty());
        return std::vector<double>{p[0] + p[1]};
    };

    auto base = std::make_shared<ChiSquaredLikelihood>(model, ctx, 2, Cinv);
    auto constraint = make_constraint_dist(1.0, 2.0);
    WithGaussianConstraints wrapped(base, constraint, {1});

    {
        const std::vector<double> theta{3.0, 1.0};
        const double base_expected = 0.5 * 2.0 * std::pow(4.0, 2.0);
        const double constraint_expected = 0.0;
        assert(approx(wrapped.nll(theta), base_expected + constraint_expected));
    }

    {
        const std::vector<double> theta{3.0, 5.0};
        const double base_expected = 0.5 * 2.0 * std::pow(8.0, 2.0);
        const double constraint_expected = 0.5 * std::pow((5.0 - 1.0) / 2.0, 2.0);
        assert(approx(wrapped.nll(theta), base_expected + constraint_expected));
    }

    {
        assert(wrapped.dim() == 2);
        auto defs = wrapped.get_param_defs();
        assert(defs.size() == 2);
        assert(defs[0].name == "x");
        assert(defs[1].name == "y");
    }

    std::cout << "\nWithGaussianConstraints INTEGRATION tests passed!\n";
    return 0;
}
