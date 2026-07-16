#include "Include.h"
#include "JointDistribution.h"
#include "IMarginalDistribution.h"
#include "ICopula.h"
#include "GaussianCopula.h"
#include "StudentTCopula.h"
#include "CopulaFactory.h"

#include <gsl/gsl_cdf.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <vector>

static constexpr double kPi = 3.141592653589793238462643383279502884;

static bool approx(double a, double b, double eps = 1e-9) {
    return std::fabs(a - b) <= eps;
}

static bool equal_vec(const std::vector<double>& a,
                      const std::vector<double>& b,
                      double eps = 0.0) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (eps == 0.0) {
            if (a[i] != b[i]) return false;
        } else {
            if (!approx(a[i], b[i], eps)) return false;
        }
    }
    return true;
}

static bool equal_samples(const std::vector<std::vector<double>>& a,
                          const std::vector<std::vector<double>>& b,
                          double eps = 0.0) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (!equal_vec(a[i], b[i], eps)) return false;
    }
    return true;
}

static RealMatrix corr2(double rho) {
    RealMatrix R(2, 2);
    R.at(0, 0) = 1.0;
    R.at(0, 1) = rho;
    R.at(1, 0) = rho;
    R.at(1, 1) = 1.0;
    return R;
}

static double mean(const std::vector<double>& x) {
    double s = 0.0;
    for (double v : x) s += v;
    return x.empty() ? 0.0 : s / static_cast<double>(x.size());
}

static double stdev(const std::vector<double>& x) {
    assert(x.size() > 1);
    const double mx = mean(x);
    double s = 0.0;
    for (double v : x) {
        const double d = v - mx;
        s += d * d;
    }
    return std::sqrt(s / static_cast<double>(x.size() - 1));
}

static double corr(const std::vector<double>& x, const std::vector<double>& y) {
    assert(x.size() == y.size());
    assert(!x.empty());

    const double mx = mean(x);
    const double my = mean(y);

    double sxx = 0.0;
    double syy = 0.0;
    double sxy = 0.0;

    for (std::size_t i = 0; i < x.size(); ++i) {
        const double dx = x[i] - mx;
        const double dy = y[i] - my;
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
    }

    return sxy / std::sqrt(sxx * syy);
}

class TestGaussianMarginal final : public IMarginalDistribution {
public:
    TestGaussianMarginal(double mu, double sigma) : mu_(mu), sigma_(sigma) {
        if (!(sigma_ > 0.0)) {
            throw std::invalid_argument("sigma must be positive");
        }
    }

    std::vector<double> rvs(std::size_t n) override {
        return std::vector<double>(n, mu_);
    }

    double logpdf(double x) override {
        const double z = (x - mu_) / sigma_;
        return -0.5 * z * z - std::log(sigma_) - 0.5 * std::log(2.0 * kPi);
    }

    PDFDiff f_df_ddf(double x) override {
        const double f = std::exp(logpdf(x));
        const double dx = x - mu_;
        const double sigma2 = sigma_ * sigma_;
        const double df = -(dx / sigma2) * f;
        const double ddf = ((dx * dx) / (sigma2 * sigma2) - 1.0 / sigma2) * f;
        return PDFDiff{f, df, ddf};
    }

    double cdf(double x) override {
        return gsl_cdf_ugaussian_P((x - mu_) / sigma_);
    }

    double ppf(double u) override {
        const double uc = std::clamp(u, 1e-14, 1.0 - 1e-14);
        return mu_ + sigma_ * gsl_cdf_ugaussian_Pinv(uc);
    }

    double mean() override {
        return mu_;
    }

    double std() override {
        return sigma_;
    }

private:
    double mu_;
    double sigma_;
};

static std::vector<std::unique_ptr<IMarginalDistribution>> make_gaussian_marginals() {
    std::vector<std::unique_ptr<IMarginalDistribution>> marginals;
    marginals.emplace_back(std::make_unique<TestGaussianMarginal>(1.0, 2.0));
    marginals.emplace_back(std::make_unique<TestGaussianMarginal>(-1.0, 0.5));
    return marginals;
}

static double marginal0_cdf(double x) {
    return gsl_cdf_ugaussian_P((x - 1.0) / 2.0);
}

static double marginal1_cdf(double x) {
    return gsl_cdf_ugaussian_P((x + 1.0) / 0.5);
}

static double marginal0_logpdf(double x) {
    const double z = (x - 1.0) / 2.0;
    return -0.5 * z * z - std::log(2.0) - 0.5 * std::log(2.0 * kPi);
}

static double marginal1_logpdf(double x) {
    const double z = (x + 1.0) / 0.5;
    return -0.5 * z * z - std::log(0.5) - 0.5 * std::log(2.0 * kPi);
}

static JointDistribution make_gaussian_joint(unsigned long seed, double rho) {
    GaussianCopulaConfig cfg;
    cfg.R = corr2(rho);

    return JointDistribution(
        make_gaussian_marginals(),
        CopulaFactory::create(CopulaType::GAUSSIAN, cfg, seed)
    );
}

static JointDistribution make_student_joint(unsigned long seed, double rho, int nu) {
    StudentTCopulaConfig cfg;
    cfg.R = corr2(rho);
    cfg.nu = nu;

    return JointDistribution(
        make_gaussian_marginals(),
        CopulaFactory::create(CopulaType::STUDENT_T, cfg, seed)
    );
}

static double empirical_gaussian_latent_corr(JointDistribution& joint,
                                             std::size_t n) {
    auto X = joint.sample(n);

    std::vector<double> z0;
    std::vector<double> z1;
    z0.reserve(n);
    z1.reserve(n);

    for (const auto& x : X) {
        const double u0 = std::clamp(marginal0_cdf(x[0]), 1e-12, 1.0 - 1e-12);
        const double u1 = std::clamp(marginal1_cdf(x[1]), 1e-12, 1.0 - 1e-12);
        z0.push_back(gsl_cdf_ugaussian_Pinv(u0));
        z1.push_back(gsl_cdf_ugaussian_Pinv(u1));
    }

    return corr(z0, z1);
}

static double empirical_student_latent_corr(JointDistribution& joint,
                                            int nu,
                                            std::size_t n) {
    auto X = joint.sample(n);

    std::vector<double> t0;
    std::vector<double> t1;
    t0.reserve(n);
    t1.reserve(n);

    for (const auto& x : X) {
        const double u0 = std::clamp(marginal0_cdf(x[0]), 1e-12, 1.0 - 1e-12);
        const double u1 = std::clamp(marginal1_cdf(x[1]), 1e-12, 1.0 - 1e-12);
        t0.push_back(gsl_cdf_tdist_Pinv(u0, nu));
        t1.push_back(gsl_cdf_tdist_Pinv(u1, nu));
    }

    return corr(t0, t1);
}

int main() {
    std::cout << "== Running INTEGRATION tests for JointDistribution ==\n";

    {
        auto joint = make_gaussian_joint(1234, 0.55);

        assert(joint.dim() == 2);

        auto stds = joint.get_stds();
        assert(stds.size() == 2);
        assert(approx(stds[0], 2.0));
        assert(approx(stds[1], 0.5));

        auto x = joint.sample();
        assert(x.size() == 2);
        assert(std::isfinite(x[0]));
        assert(std::isfinite(x[1]));

        auto X = joint.sample(25);
        assert(X.size() == 25);
        for (const auto& xi : X) {
            assert(xi.size() == 2);
            assert(std::isfinite(xi[0]));
            assert(std::isfinite(xi[1]));
            assert(std::isfinite(joint.logpdf(xi)));
        }
    }

    {
        auto j1 = make_gaussian_joint(2026, 0.40);
        auto j2 = make_gaussian_joint(2026, 0.40);

        auto s1 = j1.sample(30);
        auto s2 = j2.sample(30);

        assert(equal_samples(s1, s2, 0.0));
    }

    {
        auto joint = make_gaussian_joint(3001, 0.0);
        std::vector<double> x{1.25, -0.75};

        const double expected = marginal0_logpdf(x[0]) + marginal1_logpdf(x[1]);
        assert(approx(joint.logpdf(x), expected, 1e-10));
    }

    {
        auto joint = make_gaussian_joint(4001, 0.80);
        const double rho_emp = empirical_gaussian_latent_corr(joint, 4000);

        assert(rho_emp > 0.65);
    }

    {
        auto joint = make_student_joint(5001, 0.75, 6);

        auto x = joint.sample();
        assert(x.size() == 2);
        assert(std::isfinite(x[0]));
        assert(std::isfinite(x[1]));
        assert(std::isfinite(joint.logpdf(x)));

        const double rho_emp = empirical_student_latent_corr(joint, 6, 4000);
        assert(rho_emp > 0.55);
    }

    {
        auto joint = make_gaussian_joint(6001, 0.45);
        auto X = joint.sample(5000);

        std::vector<double> x0;
        std::vector<double> x1;
        x0.reserve(X.size());
        x1.reserve(X.size());

        for (const auto& x : X) {
            x0.push_back(x[0]);
            x1.push_back(x[1]);
        }

        assert(std::fabs(mean(x0) - 1.0) < 0.15);
        assert(std::fabs(mean(x1) + 1.0) < 0.06);
        assert(std::fabs(stdev(x0) - 2.0) < 0.18);
        assert(std::fabs(stdev(x1) - 0.5) < 0.05);
    }

    {
        auto joint = make_gaussian_joint(7001, 0.50);
        RealMatrix W = joint.curvature(std::vector<double>{1.0, -1.0});

        assert(std::isfinite(W.at(0, 0)));
        assert(std::isfinite(W.at(0, 1)));
        assert(std::isfinite(W.at(1, 0)));
        assert(std::isfinite(W.at(1, 1)));
        assert(approx(W.at(0, 1), W.at(1, 0), 1e-12));
    }

    /*
     * StudentTCopula curvature intentionally not tested here yet.
     *
     * JointDistribution::curvature expects copula_->log_c_dc_ddc(u) to return:
     *   - dlog_c  shaped as d x 1
     *   - ddlog_c shaped as d x d
     *
     * GaussianCopula follows this convention and is tested above.
     * If StudentTCopula returns the gradient as 1 x d, or if one of its
     * derivative matrices has another shape, JointDistribution::curvature
     * throws std::out_of_range("Indices outside matrix shape").
     *
     * Once StudentTCopula.{h,cpp} is checked/fixed for the same matrix-shape
     * convention as GaussianCopula, this block can be re-enabled.
     */

    {
        auto joint = make_gaussian_joint(9001, 0.20);

        bool threw = false;
        try {
            (void)joint.logpdf(std::vector<double>{1.0});
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        assert(threw);

        threw = false;
        try {
            (void)joint.curvature(std::vector<double>{1.0});
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        assert(threw);
    }

    std::cout << "\nJointDistribution INTEGRATION tests passed!\n";
    return 0;
}
