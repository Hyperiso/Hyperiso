#include "Include.h"
#include "ICopula.h"
#include "GaussianCopula.h"
#include "StudentTCopula.h"
#include "CopulaFactory.h"

#include <gsl/gsl_cdf.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>
#include <numeric>

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

static bool equal_samples(const std::vector<Vector>& a,
                          const std::vector<Vector>& b,
                          double eps = 0.0) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (!equal_vec(a[i], b[i], eps)) return false;
    }
    return true;
}

static bool in_unit_interval(double x) {
    return x >= 0.0 && x <= 1.0;
}

static bool in_unit_cube(const Vector& u) {
    for (double x : u) {
        if (!in_unit_interval(x)) return false;
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

static double corr(const std::vector<double>& x, const std::vector<double>& y) {
    assert(x.size() == y.size());
    assert(!x.empty());

    double mx = mean(x);
    double my = mean(y);

    double sxx = 0.0, syy = 0.0, sxy = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        double dx = x[i] - mx;
        double dy = y[i] - my;
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
    }

    return sxy / std::sqrt(sxx * syy);
}

static double empirical_gaussian_latent_corr(ICopula& cop, std::size_t n) {
    std::vector<double> x, y;
    x.reserve(n);
    y.reserve(n);

    for (std::size_t i = 0; i < n; ++i) {
        Vector u = cop.sample_u();
        x.push_back(gsl_cdf_ugaussian_Pinv(std::clamp(u[0], 1e-12, 1.0 - 1e-12)));
        y.push_back(gsl_cdf_ugaussian_Pinv(std::clamp(u[1], 1e-12, 1.0 - 1e-12)));
    }

    return corr(x, y);
}

static double empirical_student_latent_corr(ICopula& cop, int nu, std::size_t n) {
    std::vector<double> x, y;
    x.reserve(n);
    y.reserve(n);

    for (std::size_t i = 0; i < n; ++i) {
        Vector u = cop.sample_u();
        x.push_back(gsl_cdf_tdist_Pinv(std::clamp(u[0], 1e-12, 1.0 - 1e-12), nu));
        y.push_back(gsl_cdf_tdist_Pinv(std::clamp(u[1], 1e-12, 1.0 - 1e-12), nu));
    }

    return corr(x, y);
}

int main() {
    std::cout << "== Running INTEGRATION tests for Copulas ==\n";

    {
        GaussianCopulaConfig gcfg;
        gcfg.R = corr2(0.65);

        StudentTCopulaConfig tcfg;
        tcfg.R  = corr2(0.65);
        tcfg.nu = 6;

        std::vector<std::unique_ptr<ICopula>> cops;
        cops.push_back(CopulaFactory::create(CopulaType::GAUSSIAN, gcfg, 111));
        cops.push_back(CopulaFactory::create(CopulaType::STUDENT_T, tcfg, 222));

        assert(cops.size() == 2);

        for (auto& cop : cops) {
            auto u = cop->sample_u();
            assert(u.size() == 2);
            assert(in_unit_cube(u));

            auto U = cop->sample_u(20);
            assert(U.size() == 20);

            for (const auto& ui : U) {
                assert(ui.size() == 2);
                assert(in_unit_cube(ui));
                assert(std::isfinite(cop->log_density(ui)));
            }
        }
    }

    {
        GaussianCopulaConfig gcfg;
        gcfg.R = corr2(0.45);

        auto c1 = CopulaFactory::create(CopulaType::GAUSSIAN, gcfg, 12345);
        auto c2 = CopulaFactory::create(CopulaType::GAUSSIAN, gcfg, 12345);

        auto s1 = c1->sample_u(25);
        auto s2 = c2->sample_u(25);

        assert(equal_samples(s1, s2));

        StudentTCopulaConfig tcfg;
        tcfg.R  = corr2(0.45);
        tcfg.nu = 8;

        auto t1 = CopulaFactory::create(CopulaType::STUDENT_T, tcfg, 54321);
        auto t2 = CopulaFactory::create(CopulaType::STUDENT_T, tcfg, 54321);

        auto r1 = t1->sample_u(25);
        auto r2 = t2->sample_u(25);

        assert(equal_samples(r1, r2));
    }

    {
        GaussianCopulaConfig cfg;
        cfg.R = corr2(0.80);

        auto cop = CopulaFactory::create(CopulaType::GAUSSIAN, cfg, 2026);
        double rho_emp = empirical_gaussian_latent_corr(*cop, 4000);

        assert(rho_emp > 0.55);
    }

    {
        StudentTCopulaConfig cfg;
        cfg.R  = corr2(0.80);
        cfg.nu = 5;

        auto cop = CopulaFactory::create(CopulaType::STUDENT_T, cfg, 2027);
        double rho_emp = empirical_student_latent_corr(*cop, cfg.nu, 4000);

        assert(rho_emp > 0.50);
    }

    {
        GaussianCopulaConfig gcfg;
        gcfg.R = corr2(0.0);

        auto gcop = CopulaFactory::create(CopulaType::GAUSSIAN, gcfg, 3001);
        double grho = empirical_gaussian_latent_corr(*gcop, 4000);
        assert(std::fabs(grho) < 0.08);

        StudentTCopulaConfig tcfg;
        tcfg.R  = corr2(0.0);
        tcfg.nu = 10;

        auto tcop = CopulaFactory::create(CopulaType::STUDENT_T, tcfg, 3002);
        double trho = empirical_student_latent_corr(*tcop, tcfg.nu, 4000);
        assert(std::fabs(trho) < 0.08);
    }

    {
        GaussianCopulaConfig gcfg;
        gcfg.R = corr2(0.0);

        auto gcop = CopulaFactory::create(CopulaType::GAUSSIAN, gcfg, 4001);
        assert(approx(gcop->log_density(Vector{0.3, 0.9}), 0.0, 1e-10));

        StudentTCopulaConfig tcfg;
        tcfg.R  = corr2(0.0);
        tcfg.nu = 7;

        auto tcop = CopulaFactory::create(CopulaType::STUDENT_T, tcfg, 4002);
        assert(approx(tcop->log_density(Vector{0.3, 0.9}), 0.0, 1e-8));
    }

    {
        StudentTCopulaConfig cfg;
        cfg.R  = corr2(0.60);
        cfg.nu = 4;

        auto a = CopulaFactory::create(CopulaType::GAUSSIAN, cfg, 555);
        auto b = CopulaFactory::create(CopulaType::STUDENT_T, cfg, 555);

        auto sa = a->sample_u(30);
        auto sb = b->sample_u(30);

        assert(equal_samples(sa, sb));
    }

    std::cout << "\nCopula INTEGRATION tests passed!\n";
    return 0;
}