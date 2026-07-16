#include "Include.h"
#include "WithGaussianConstraints.h"
#include "ICopula.h"
#include "IMarginalDistribution.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
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

class QuadraticMarginal : public IMarginalDistribution {
public:
    explicit QuadraticMarginal(double mu = 0.0, double sigma = 1.0)
        : mu_(mu), sigma_(sigma)
    {}

    std::vector<double> rvs(std::size_t n) override {
        return std::vector<double>(n, mu_);
    }

    double logpdf(double x) override {
        const double z = (x - mu_) / sigma_;
        return -0.5 * z * z;
    }

    PDFDiff f_df_ddf(double x) override {
        const double z = (x - mu_) / sigma_;
        const double s2 = sigma_ * sigma_;
        const double f = std::exp(-0.5 * z * z);
        const double df = -(x - mu_) / s2 * f;
        const double ddf = (z * z - 1.0) / s2 * f;
        return {f, df, ddf};
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

    std::vector<double> sample_u() override {
        return std::vector<double>(d_, 0.5);
    }

    double log_density(std::vector<double>) override { return 0.0; }
    RealMatrix dlog_density(std::vector<double>) override { return RealMatrix(d_, 1); }
    RealMatrix ddlog_density(std::vector<double>) override { return RealMatrix(d_, d_); }
    LogDensityDiff log_c_dc_ddc(std::vector<double>) override {
        return {0.0, RealMatrix(d_, 1), RealMatrix(d_, d_)};
    }

private:
    std::size_t d_;
};

static std::shared_ptr<JointDistribution> make_constraint_dist(
    const std::vector<double>& mu,
    const std::vector<double>& sigma
) {
    std::vector<std::unique_ptr<IMarginalDistribution>> marginals;
    for (std::size_t i = 0; i < mu.size(); ++i) {
        marginals.push_back(std::make_unique<QuadraticMarginal>(mu[i], sigma[i]));
    }
    return std::make_shared<JointDistribution>(
        std::move(marginals),
        std::make_unique<IndependentCopula>(mu.size())
    );
}

class BaseSpyLikelihood : public ILikelihood {
public:
    mutable int nll_calls = 0;
    mutable std::vector<double> last_theta;

    double nll(const std::vector<double>& theta) const override {
        nll_calls++;
        last_theta = theta;
        return 7.0;
    }

    std::vector<fit_app::ParameterDefinition> get_param_defs() const override {
        return {make_def("a", 1.0), make_def("b", 2.0), make_def("c", 3.0)};
    }

    std::size_t dim() const override { return 3; }
};

int main() {
    std::cout << "== WithGaussianConstraints UNIT ==\n";

    {
        auto base = std::make_shared<BaseSpyLikelihood>();
        auto dist = make_constraint_dist({0.0, 10.0}, {1.0, 2.0});
        WithGaussianConstraints wrapped(base, dist, {0, 2});

        const std::vector<double> theta{2.0, 99.0, 12.0};
        const double expected_constraint = 0.5 * 2.0 * 2.0
                                         + 0.5 * std::pow((12.0 - 10.0) / 2.0, 2.0);
        const double out = wrapped.nll(theta);

        assert(approx(out, 7.0 + expected_constraint));
        assert(base->nll_calls == 1);
        assert((base->last_theta == theta));
    }

    {
        auto base = std::make_shared<BaseSpyLikelihood>();
        auto dist = make_constraint_dist({0.0}, {1.0});
        WithGaussianConstraints wrapped(base, dist, {1});

        assert(wrapped.dim() == base->dim());
        auto defs = wrapped.get_param_defs();
        assert(defs.size() == 3);
        assert(defs[0].name == "a");
        assert(defs[1].name == "b");
        assert(defs[2].name == "c");
    }

    {
        // The current implementation indexes theta with operator[] rather than at(),
        // so an out-of-range constrained index is undefined behavior and should not
        // be asserted as throwing. The well-defined edge case is an empty constraint
        // list, which should simply forward the base NLL unchanged.
        auto base = std::make_shared<BaseSpyLikelihood>();
        auto dist = make_constraint_dist({}, {});
        WithGaussianConstraints wrapped(base, dist, {});

        const std::vector<double> theta{1.0, 2.0, 3.0};
        const double out = wrapped.nll(theta);

        assert(approx(out, 7.0));
        assert(base->nll_calls == 1);
        assert((base->last_theta == theta));
    }

    std::cout << "\nWithGaussianConstraints UNIT tests passed!\n";
    return 0;
}
