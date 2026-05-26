#include "Include.h"
#include "JointDistribution.h"
#include "IMarginalDistribution.h"
#include "ICopula.h"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

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

class LinearTestMarginal final : public IMarginalDistribution {
public:
    LinearTestMarginal(double cdf_intercept,
                       double cdf_slope,
                       double ppf_intercept,
                       double ppf_slope,
                       double logpdf_value,
                       double mean_value,
                       double std_value,
                       PDFDiff diff) :
        cdf_intercept_(cdf_intercept),
        cdf_slope_(cdf_slope),
        ppf_intercept_(ppf_intercept),
        ppf_slope_(ppf_slope),
        logpdf_value_(logpdf_value),
        mean_value_(mean_value),
        std_value_(std_value),
        diff_(diff)
    {}

    std::vector<double> rvs(std::size_t n) override {
        return std::vector<double>(n, mean_value_);
    }

    double logpdf(double /*x*/) override {
        return logpdf_value_;
    }

    PDFDiff f_df_ddf(double /*x*/) override {
        return diff_;
    }

    double cdf(double x) override {
        return cdf_intercept_ + cdf_slope_ * x;
    }

    double ppf(double u) override {
        return ppf_intercept_ + ppf_slope_ * u;
    }

    double mean() override {
        return mean_value_;
    }

    double std() override {
        return std_value_;
    }

private:
    double cdf_intercept_;
    double cdf_slope_;
    double ppf_intercept_;
    double ppf_slope_;
    double logpdf_value_;
    double mean_value_;
    double std_value_;
    PDFDiff diff_;
};

class DeterministicCopula final : public ICopula {
public:
    DeterministicCopula(std::vector<double> single_u,
                        std::vector<std::vector<double>> batch_u,
                        double log_density_value,
                        LogDensityDiff diff) :
        single_u_(std::move(single_u)),
        batch_u_(std::move(batch_u)),
        log_density_value_(log_density_value),
        diff_(std::move(diff))
    {}

    std::vector<std::vector<double>> sample_u(std::size_t n) override {
        assert(n == batch_u_.size());
        return batch_u_;
    }

    std::vector<double> sample_u() override {
        return single_u_;
    }

    double log_density(std::vector<double> u) override {
        last_log_density_u = std::move(u);
        ++log_density_calls;
        return log_density_value_;
    }

    RealMatrix dlog_density(std::vector<double> /*u*/) override {
        return diff_.dlog_c;
    }

    RealMatrix ddlog_density(std::vector<double> /*u*/) override {
        return diff_.ddlog_c;
    }

    LogDensityDiff log_c_dc_ddc(std::vector<double> u) override {
        last_diff_u = std::move(u);
        ++diff_calls;
        return diff_;
    }

    std::vector<double> last_log_density_u;
    std::vector<double> last_diff_u;
    std::size_t log_density_calls{0};
    std::size_t diff_calls{0};

private:
    std::vector<double> single_u_;
    std::vector<std::vector<double>> batch_u_;
    double log_density_value_;
    LogDensityDiff diff_;
};

static std::vector<std::unique_ptr<IMarginalDistribution>> make_simple_marginals() {
    std::vector<std::unique_ptr<IMarginalDistribution>> marginals;
    marginals.emplace_back(std::make_unique<LinearTestMarginal>(
        0.0, 1.0,
        10.0, 2.0,
        -1.2,
        10.0, 2.0,
        PDFDiff{2.0, 0.5, 0.25}
    ));
    marginals.emplace_back(std::make_unique<LinearTestMarginal>(
        0.0, 1.0,
        -1.0, 3.0,
        -0.7,
        -1.0, 3.0,
        PDFDiff{3.0, -0.3, 0.12}
    ));
    return marginals;
}

static LogDensityDiff make_zero_diff() {
    LogDensityDiff diff;
    diff.log_c = 0.0;
    diff.dlog_c = RealMatrix(2, 1);
    diff.ddlog_c = RealMatrix(2, 2);

    for (std::size_t i = 0; i < 2; ++i) {
        diff.dlog_c.at(i, 0) = 0.0;
        for (std::size_t j = 0; j < 2; ++j) {
            diff.ddlog_c.at(i, j) = 0.0;
        }
    }

    return diff;
}

static LogDensityDiff make_curvature_diff() {
    LogDensityDiff diff;
    diff.log_c = 0.0;
    diff.dlog_c = RealMatrix(2, 1);
    diff.ddlog_c = RealMatrix(2, 2);

    diff.dlog_c.at(0, 0) = 0.7;
    diff.dlog_c.at(1, 0) = -0.4;

    diff.ddlog_c.at(0, 0) = -1.0;
    diff.ddlog_c.at(0, 1) = 0.2;
    diff.ddlog_c.at(1, 0) = 0.4;
    diff.ddlog_c.at(1, 1) = -2.0;

    return diff;
}

int main() {
    std::cout << "== Running UNIT tests for JointDistribution ==\n";

    {
        auto copula = std::make_unique<DeterministicCopula>(
            std::vector<double>{0.25, 0.80},
            std::vector<std::vector<double>>{},
            0.0,
            make_zero_diff()
        );

        JointDistribution joint(make_simple_marginals(), std::move(copula));
        auto x = joint.sample();

        assert(x.size() == 2);
        assert(approx(x[0], 10.0 + 2.0 * 0.25));
        assert(approx(x[1], -1.0 + 3.0 * 0.80));
    }

    {
        std::vector<std::vector<double>> U{
            {0.10, 0.20},
            {0.70, 0.90},
            {0.40, 0.60}
        };

        auto copula = std::make_unique<DeterministicCopula>(
            std::vector<double>{},
            U,
            0.0,
            make_zero_diff()
        );

        JointDistribution joint(make_simple_marginals(), std::move(copula));
        auto X = joint.sample(3);

        std::vector<std::vector<double>> expected{
            {10.0 + 2.0 * 0.10, -1.0 + 3.0 * 0.20},
            {10.0 + 2.0 * 0.70, -1.0 + 3.0 * 0.90},
            {10.0 + 2.0 * 0.40, -1.0 + 3.0 * 0.60}
        };

        assert(equal_samples(X, expected, 1e-12));
    }

    {
        auto copula = std::make_unique<DeterministicCopula>(
            std::vector<double>{},
            std::vector<std::vector<double>>{},
            0.4,
            make_zero_diff()
        );
        auto* raw_copula = copula.get();

        JointDistribution joint(make_simple_marginals(), std::move(copula));
        double lp = joint.logpdf(std::vector<double>{0.20, 0.80});

        assert(approx(lp, -1.2 - 0.7 + 0.4));
        assert(raw_copula->log_density_calls == 1);
        assert(equal_vec(raw_copula->last_log_density_u,
                         std::vector<double>{0.20, 0.80},
                         1e-15));
    }

    {
        auto copula = std::make_unique<DeterministicCopula>(
            std::vector<double>{},
            std::vector<std::vector<double>>{},
            1.0,
            make_zero_diff()
        );
        auto* raw_copula = copula.get();

        JointDistribution joint(make_simple_marginals(), std::move(copula));
        (void)joint.logpdf(std::vector<double>{-10.0, 2.0});

        assert(raw_copula->last_log_density_u.size() == 2);
        assert(approx(raw_copula->last_log_density_u[0], 1e-13, 1e-18));
        assert(approx(raw_copula->last_log_density_u[1], 1.0 - 1e-13, 1e-18));
    }

    {
        auto copula = std::make_unique<DeterministicCopula>(
            std::vector<double>{},
            std::vector<std::vector<double>>{},
            0.0,
            make_zero_diff()
        );

        JointDistribution joint(make_simple_marginals(), std::move(copula));

        bool threw = false;
        try {
            (void)joint.logpdf(std::vector<double>{1.0});
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        assert(threw);
    }

    {
        auto copula = std::make_unique<DeterministicCopula>(
            std::vector<double>{},
            std::vector<std::vector<double>>{},
            0.0,
            make_curvature_diff()
        );
        auto* raw_copula = copula.get();

        JointDistribution joint(make_simple_marginals(), std::move(copula));
        RealMatrix W = joint.curvature(std::vector<double>{0.25, 0.75});

        assert(raw_copula->diff_calls == 1);
        assert(equal_vec(raw_copula->last_diff_u,
                         std::vector<double>{0.25, 0.75},
                         1e-15));

        assert(approx(W.at(0, 0), 3.5875, 1e-12));
        assert(approx(W.at(1, 1), 17.85, 1e-12));
        assert(approx(W.at(0, 1), -1.8, 1e-12));
        assert(approx(W.at(1, 0), -1.8, 1e-12));
    }

    {
        auto copula = std::make_unique<DeterministicCopula>(
            std::vector<double>{},
            std::vector<std::vector<double>>{},
            0.0,
            make_zero_diff()
        );

        JointDistribution joint(make_simple_marginals(), std::move(copula));

        bool threw = false;
        try {
            (void)joint.curvature(std::vector<double>{1.0});
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        assert(threw);
    }

    {
        auto copula = std::make_unique<DeterministicCopula>(
            std::vector<double>{0.5, 0.5},
            std::vector<std::vector<double>>{},
            0.0,
            make_zero_diff()
        );

        JointDistribution joint(make_simple_marginals(), std::move(copula));

        assert(joint.dim() == 2);

        auto stds = joint.get_stds();
        assert(stds.size() == 2);
        assert(approx(stds[0], 2.0));
        assert(approx(stds[1], 3.0));
    }

    std::cout << "\nAll JointDistribution UNIT tests passed!\n";
    return 0;
}
