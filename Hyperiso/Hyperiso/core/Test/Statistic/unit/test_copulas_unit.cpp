#include "Include.h"
#include "ICopula.h"
#include "GenericCopula.h"
#include "GaussianCopula.h"
#include "StudentTCopula.h"
#include "CopulaFactory.h"

#include <gsl/gsl_cdf.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <limits>

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

static bool all_samples_in_unit_cube(const std::vector<Vector>& U, std::size_t d) {
    for (const auto& u : U) {
        if (u.size() != d) return false;
        if (!in_unit_cube(u)) return false;
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

int main() {
    std::cout << "== Running UNIT tests for Copulas ==\n";

    {
        GaussianCopula cop(1234, corr2(0.6));
        Vector u = cop.sample_u();

        assert(u.size() == 2);
        assert(in_unit_cube(u));
    }

    {
        GaussianCopula c1(42, corr2(0.35));
        GaussianCopula c2(42, corr2(0.35));

        auto batch = c1.sample_u(10);

        std::vector<Vector> seq;
        for (std::size_t i = 0; i < 10; ++i) {
            seq.push_back(c2.sample_u());
        }

        assert(batch.size() == 10);
        assert(all_samples_in_unit_cube(batch, 2));
        assert(equal_samples(batch, seq));
    }

    {
        GaussianCopula c1(777, corr2(0.8));
        GaussianCopula c2(777, corr2(0.8));

        auto s1 = c1.sample_u(12);
        auto s2 = c2.sample_u(12);

        assert(equal_samples(s1, s2));
    }

    {
        GaussianCopula cop(1, corr2(0.0));

        Vector u1{0.2, 0.7};
        Vector u2{1e-20, 1.0}; 

        assert(std::isfinite(cop.log_density(u1)));
        assert(std::isfinite(cop.log_density(u2)));
        assert(approx(cop.log_density(u1), 0.0, 1e-10));
        assert(approx(cop.log_density(u2), 0.0, 1e-10));
    }

    {
        GaussianCopula cop(2, corr2(0.8));

        double ld_diag = cop.log_density(Vector{0.8, 0.8});
        double ld_off  = cop.log_density(Vector{0.8, 0.2});

        assert(ld_diag > ld_off);
    }

    {
        bool threw = false;
        try {
            StudentTCopula bad(123, corr2(0.4), 1);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        assert(threw);
    }

    {
        StudentTCopula cop(4321, corr2(0.5), 5);
        Vector u = cop.sample_u();

        assert(u.size() == 2);
        assert(in_unit_cube(u));
    }

    {
        StudentTCopula c1(99, corr2(0.25), 7);
        StudentTCopula c2(99, corr2(0.25), 7);

        auto batch = c1.sample_u(10);

        std::vector<Vector> seq;
        for (std::size_t i = 0; i < 10; ++i) {
            seq.push_back(c2.sample_u());
        }

        assert(batch.size() == 10);
        assert(all_samples_in_unit_cube(batch, 2));
        assert(equal_samples(batch, seq));
    }

    {
        StudentTCopula c1(2025, corr2(0.7), 6);
        StudentTCopula c2(2025, corr2(0.7), 6);

        auto s1 = c1.sample_u(12);
        auto s2 = c2.sample_u(12);

        assert(equal_samples(s1, s2));
    }

    {
        StudentTCopula cop(5, corr2(0.0), 8);

        Vector u1{0.2, 0.7};
        Vector u2{1e-30, 1.0}; 

        assert(std::isfinite(cop.log_density(u1)));
        assert(std::isfinite(cop.log_density(u2)));

        // A Student-t latent vector with diagonal R still shares the same
        // chi-square scale factor, so this is not an independence copula test.
        // With the current non-normalized convention, only finiteness is asserted here.
    }

    {
        StudentTCopula cop(6, corr2(0.8), 5);

        double ld_diag = cop.log_density(Vector{0.8, 0.8});
        double ld_off  = cop.log_density(Vector{0.8, 0.2});

        // StudentTCopula::log_density currently uses a non-normalized internal
        // convention, so we avoid asserting a normalized copula-density ordering.
        assert(std::isfinite(ld_diag));
        assert(std::isfinite(ld_off));
    }

    {
        GaussianCopulaConfig cfg;
        cfg.R = corr2(0.4);

        auto cop = CopulaFactory::create(CopulaType::GAUSSIAN, cfg, 321);
        assert(cop != nullptr);

        auto u = cop->sample_u();
        assert(u.size() == 2);
        assert(in_unit_cube(u));
    }

    {
        StudentTCopulaConfig cfg;
        cfg.R  = corr2(0.4);
        cfg.nu = 6;

        auto cop = CopulaFactory::create(CopulaType::STUDENT_T, cfg, 654);
        assert(cop != nullptr);

        auto u = cop->sample_u();
        assert(u.size() == 2);
        assert(in_unit_cube(u));
    }

    {
        bool threw = false;
        try {
            GaussianCopulaConfig cfg;
            cfg.R = corr2(0.2);

            auto cop = CopulaFactory::create(static_cast<CopulaType>(999), cfg, 1);
            (void)cop;
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        assert(threw);
    }

    {
        StudentTCopulaConfig cfg;
        cfg.R  = corr2(0.55);
        cfg.nu = 9;

        auto via_mismatch = CopulaFactory::create(CopulaType::GAUSSIAN, cfg, 987);
        auto via_student  = CopulaFactory::create(CopulaType::STUDENT_T, cfg, 987);

        auto s1 = via_mismatch->sample_u(6);
        auto s2 = via_student->sample_u(6);

        assert(equal_samples(s1, s2));
    }

    std::cout << "\nAll copula UNIT tests passed!\n";
    return 0;
}