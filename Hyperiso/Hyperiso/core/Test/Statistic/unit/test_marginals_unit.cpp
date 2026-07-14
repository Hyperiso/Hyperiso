#include "Include.h"
#include "GaussianMarginal.h"
#include "FlatMarginal.h"
#include "MarginalFactory.h"
#include "MarginalConfigFactory.h"
#include "NuisanceSpec.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>

static bool approx(double a, double b, double eps = 1e-9) {
    return std::fabs(a - b) <= eps;
}

static bool equal_vec(const Vector& a, const Vector& b, double eps = 0.0) {
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

static bool all_in_range(const Vector& v, double lo, double hi) {
    for (double x : v) {
        if (x < lo || x > hi) return false;
    }
    return true;
}

int main() {
    std::cout << "== Running UNIT tests for Marginals ==\n";

    {
        GaussianMarginal g(1.5, 2.0, 1234);

        assert(approx(g.mean(), 1.5));
        assert(approx(g.std(), 2.0));
        assert(approx(g.cdf(1.5), 0.5, 1e-12));
        assert(approx(g.ppf(0.5), 1.5, 1e-12));

        const double x = 2.1;
        const double pi = std::acos(-1.0);
        const double expected_logpdf =
            -0.5 * (std::log(2.0 * pi * 2.0 * 2.0) + std::pow((x - 1.5) / 2.0, 2.0));

        assert(approx(g.logpdf(x), expected_logpdf, 1e-12));
        assert(approx(g.ppf(g.cdf(x)), x, 1e-8));
    }

    {
        GaussianMarginal g1(0.0, 1.0, 42);
        GaussianMarginal g2(0.0, 1.0, 42);

        Vector s1 = g1.rvs(8);
        Vector s2 = g2.rvs(8);

        assert(s1.size() == 8);
        assert(equal_vec(s1, s2));
    }

    {
        FlatMarginal f(-2.0, 4.0, 2024);

        assert(approx(f.mean(), 1.0));
        assert(approx(f.std(), 6.0 / std::sqrt(12.0), 1e-12));

        assert(approx(f.logpdf(0.0), std::log(1.0 / 6.0), 1e-12));
        assert(f.logpdf(-2.0) == -1e100);
        assert(f.logpdf(4.0) == -1e100);
        assert(f.logpdf(100.0) == -1e100);

        assert(approx(f.cdf(-2.0), 0.0, 1e-12));
        assert(approx(f.cdf(4.0), 1.0, 1e-12));
        assert(approx(f.ppf(0.0), -2.0, 1e-12));
        assert(approx(f.ppf(1.0), 4.0, 1e-12));

        const double x = 1.25;
        assert(approx(f.ppf(f.cdf(x)), x, 1e-12));
    }

    {
        FlatMarginal f1(-3.0, 7.0, 7);
        FlatMarginal f2(-3.0, 7.0, 7);

        Vector s1 = f1.rvs(32);
        Vector s2 = f2.rvs(32);

        assert(s1.size() == 32);
        assert(equal_vec(s1, s2));
        assert(all_in_range(s1, -3.0, 7.0));
    }

    {
        MarginalConfig cfg = GaussianMarginalCfg{2.0, 0.75};
        auto dist = MarginalFactory::create(MarginalType::GAUSSIAN, cfg, 11);

        assert(dist != nullptr);
        assert(approx(dist->mean(), 2.0));
        assert(approx(dist->std(), 0.75));
        assert(approx(dist->cdf(2.0), 0.5, 1e-12));
    }

    {
        MarginalConfig cfg = FlatMarginalCfg{-std::sqrt(3.0), std::sqrt(3.0)};
        auto dist = MarginalFactory::create(MarginalType::FLAT, cfg, 17);

        assert(dist != nullptr);
        assert(approx(dist->mean(), 0.0, 1e-12));
        assert(approx(dist->std(), 1.0, 1e-12));
        assert(dist->logpdf(-std::sqrt(3.0)) == -1e100);
    }


    {
        const ParamId runtime_pid{ParameterType::DECAY, BlockName{"B_Xs"}, LhaID{7}};
        const ParamId config_pid{BlockName{"B_Xs"}, LhaID{7}};
        const NuisanceSpec spec{config_pid, {0.9, 4.0}, MarginalType::FLAT};

        MarginalConfig cfg = MarginalConfigFactory{}.create(
            runtime_pid,
            MarginalType::FLAT,
            spec
        );
        const auto& flat_cfg = std::get<FlatMarginalCfg>(cfg);

        assert(approx(flat_cfg.a, 0.9, 1e-12));
        assert(approx(flat_cfg.b, 4.0, 1e-12));
    }

    std::cout << "\nAll marginal UNIT tests passed!\n";
    return 0;
}