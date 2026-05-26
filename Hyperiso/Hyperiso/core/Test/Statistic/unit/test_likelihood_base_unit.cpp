#include "Include.h"
#include "BaseLikelihood.h"
#include "ICopula.h"
#include "IMarginalDistribution.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <gsl/gsl_cdf.h>

static bool approx(double a, double b, double eps = 1e-9) {
    return std::fabs(a - b) <= eps;
}

static bool equal_vec(const std::vector<double>& a,
                      const std::vector<double>& b,
                      double eps = 1e-12) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (!approx(a[i], b[i], eps)) return false;
    }
    return true;
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
    std::vector<double> rvs(std::size_t n) override {
        return std::vector<double>(n, 0.0);
    }

    double logpdf(double x) override {
        return -0.5 * x * x;
    }

    PDFDiff f_df_ddf(double x) override {
        const double f = std::exp(-0.5 * x * x);
        const double df = -x * f;
        const double ddf = (x * x - 1.0) * f;
        return {f, df, ddf};
    }

    double cdf(double x) override {
        return gsl_cdf_ugaussian_P(x);
    }

    double ppf(double p) override {
        return gsl_cdf_ugaussian_Pinv(std::clamp(p, 1e-13, 1.0 - 1e-13));
    }

    double mean() override { return 0.0; }
    double std() override { return 1.0; }
};

class IndependentCopula : public ICopula {
public:
    explicit IndependentCopula(std::size_t d) : d_(d) {}

    std::vector<std::vector<double>> sample_u(std::size_t n) override {
        std::vector<std::vector<double>> out;
        out.reserve(n);
        for (std::size_t i = 0; i < n; ++i) out.push_back(sample_u());
        return out;
    }

    std::vector<double> sample_u() override {
        return std::vector<double>(d_, 0.5);
    }

    double log_density(std::vector<double>) override { return 0.0; }

    RealMatrix dlog_density(std::vector<double>) override {
        return RealMatrix(d_, 1);
    }

    RealMatrix ddlog_density(std::vector<double>) override {
        return RealMatrix(d_, d_);
    }

    LogDensityDiff log_c_dc_ddc(std::vector<double>) override {
        return {0.0, RealMatrix(d_, 1), RealMatrix(d_, d_)};
    }

private:
    std::size_t d_;
};

static std::unique_ptr<JointDistribution> make_standard_joint(std::size_t d) {
    std::vector<std::unique_ptr<IMarginalDistribution>> marginals;
    marginals.reserve(d);
    for (std::size_t i = 0; i < d; ++i) {
        marginals.push_back(std::make_unique<StandardQuadraticMarginal>());
    }
    return std::make_unique<JointDistribution>(
        std::move(marginals),
        std::make_unique<IndependentCopula>(d)
    );
}

static std::shared_ptr<LikelihoodContext> make_context() {
    auto ctx = std::make_shared<LikelihoodContext>();
    ctx->exp_obs_dist = make_standard_joint(2);
    ctx->nuisance_dist = make_standard_joint(2);
    ctx->exp_obs_values = {10.0, 20.0};
    ctx->fp_defs = {make_def("p0", 1.0), make_def("p1", 2.0)};
    ctx->nuis_defs = {make_def("eta0", 0.1), make_def("eta1", -0.2)};
    return ctx;
}

int main() {
    std::cout << "== BaseLikelihood UNIT ==\n";

    {
        auto ctx = make_context();
        ModelFn model = [](const std::vector<double>& p,
                           const std::vector<double>& eta) {
            return std::vector<double>{p[0] + eta[0], p[1] + 2.0 * eta[1]};
        };

        BaseLikelihood like(model, ctx, 2);

        assert(like.p_dimension() == 2);
        assert(like.eta_dimension() == 2);
        assert(like.dim() == 4);
        assert(equal_vec(like.central_p(), std::vector<double>{1.0, 2.0}));
        assert(equal_vec(like.central_eta(), std::vector<double>{0.1, -0.2}));

        auto defs = like.get_param_defs();
        assert(defs.size() == 4);
        assert(defs[0].name == "p0");
        assert(defs[1].name == "p1");
        assert(defs[2].name == "eta0");
        assert(defs[3].name == "eta1");

        const std::vector<double> p{3.0, 4.0};
        const std::vector<double> eta{0.5, -1.0};
        const std::vector<double> theta{3.0, 4.0, 0.5, -1.0};

        assert(equal_vec(like.predict(p, eta), std::vector<double>{3.5, 2.0}));
        assert(equal_vec(like.residuals(p, eta), std::vector<double>{-6.5, -18.0}));

        const double ell_obs = -0.5 * (6.5 * 6.5 + 18.0 * 18.0);
        const double ell_nuis = -0.5 * (0.5 * 0.5 + 1.0 * 1.0);
        const double expected_nll = -(ell_obs + ell_nuis);

        assert(approx(like.nll(theta), expected_nll, 1e-12));
        assert(approx(like.nll_from_split(p, eta), expected_nll, 1e-12));
    }

    {
        auto ctx = make_context();
        ModelFn model = [](const std::vector<double>& p,
                           const std::vector<double>& eta) {
            return std::vector<double>{p[0] + eta[0], p[1] + eta[1]};
        };
        BaseLikelihood like(model, ctx, 2);

        RealMatrix Wobs = like.observable_curvature(std::vector<double>{0.2, -0.7});
        assert(Wobs.rows() == 2 && Wobs.cols() == 2);
        assert(approx(Wobs.at(0, 0), 1.0, 1e-10));
        assert(approx(Wobs.at(1, 1), 1.0, 1e-10));
        assert(approx(Wobs.at(0, 1), 0.0, 1e-10));
        assert(approx(Wobs.at(1, 0), 0.0, 1e-10));

        RealMatrix Weta = like.nuisance_curvature(std::vector<double>{0.2, -0.7});
        assert(Weta.rows() == 2 && Weta.cols() == 2);
        assert(approx(Weta.at(0, 0), 1.0, 1e-10));
        assert(approx(Weta.at(1, 1), 1.0, 1e-10));
    }

    {
        auto ctx = make_context();
        ModelFn bad_model = [](const std::vector<double>&,
                               const std::vector<double>&) {
            return std::vector<double>{std::numeric_limits<double>::infinity(), 0.0};
        };
        BaseLikelihood like(bad_model, ctx, 2);
        assert(like.nll(std::vector<double>{1.0, 2.0, 0.0, 0.0}) == 1e100);
    }

    {
        auto ctx = make_context();
        ModelFn throwing_model = [](const std::vector<double>&,
                                    const std::vector<double>&) -> std::vector<double> {
            throw std::runtime_error("intentional model failure");
        };
        BaseLikelihood like(throwing_model, ctx, 2);
        assert(like.nll(std::vector<double>{1.0, 2.0, 0.0, 0.0}) == 1e100);
    }

    {
        auto ctx = make_context();
        ModelFn wrong_size_model = [](const std::vector<double>&,
                                      const std::vector<double>&) {
            return std::vector<double>{1.0};
        };
        BaseLikelihood like(wrong_size_model, ctx, 2);

        bool threw = false;
        try {
            (void)like.residuals(std::vector<double>{1.0, 2.0}, std::vector<double>{0.0, 0.0});
        } catch (const std::runtime_error&) {
            threw = true;
        }
        assert(threw);

        assert(like.nll(std::vector<double>{1.0, 2.0, 0.0, 0.0}) == 1e100);
    }

    {
        auto ctx = make_context();
        ModelFn model = [](const std::vector<double>& p,
                           const std::vector<double>& eta) {
            return std::vector<double>{p[0] + eta[0], p[1] + eta[1]};
        };
        BaseLikelihood like(model, ctx, 2);
        like.enable_debug_trace(2);
        assert(std::isfinite(like.nll(std::vector<double>{1.0, 2.0, 0.0, 0.0})));
        assert(std::isfinite(like.nll(std::vector<double>{1.1, 2.1, 0.1, 0.1})));
        like.disable_debug_trace();
        assert(std::isfinite(like.nll(std::vector<double>{1.2, 2.2, 0.2, 0.2})));
    }

    std::cout << "\nBaseLikelihood UNIT tests passed!\n";
    return 0;
}
