#include "Include.h"
#include "MarginalFactory.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

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

static double empirical_mean(const Vector& v) {
    double s = 0.0;
    for (double x : v) s += x;
    return v.empty() ? 0.0 : s / static_cast<double>(v.size());
}

int main() {
    std::cout << "== Running INTEGRATION tests for Marginals ==\n";

    {
        std::vector<std::unique_ptr<IMarginalDistribution>> dists;
        dists.push_back(
            MarginalFactory::create(MarginalType::GAUSSIAN, GaussianMarginalCfg{10.0, 2.0}, 123)
        );
        dists.push_back(
            MarginalFactory::create(MarginalType::FLAT, FlatMarginalCfg{-std::sqrt(3.0), std::sqrt(3.0)}, 456)
        );

        assert(dists.size() == 2);

        for (auto& d : dists) {
            for (double p : {0.1, 0.25, 0.5, 0.75, 0.9}) {
                const double q = d->ppf(p);
                const double p_back = d->cdf(q);
                assert(approx(p_back, p, 1e-8));
            }

            Vector s = d->rvs(64);
            assert(s.size() == 64);
        }
    }

    {
        auto d1 = MarginalFactory::create(MarginalType::GAUSSIAN, GaussianMarginalCfg{0.0, 1.0}, 999);
        auto d2 = MarginalFactory::create(MarginalType::GAUSSIAN, GaussianMarginalCfg{0.0, 1.0}, 999);

        Vector s1 = d1->rvs(16);
        Vector s2 = d2->rvs(16);

        assert(equal_vec(s1, s2));
    }

    {
        auto flat = MarginalFactory::create(MarginalType::FLAT, FlatMarginalCfg{-2.0, 2.0}, 77);
        auto gaus = MarginalFactory::create(MarginalType::GAUSSIAN, GaussianMarginalCfg{0.0, 1.0}, 77);

        Vector sf = flat->rvs(2000);
        Vector sg = gaus->rvs(4000);

        assert(all_in_range(sf, -2.0, 2.0));
        assert(std::fabs(empirical_mean(sf) - flat->mean()) < 0.10);
        assert(std::fabs(empirical_mean(sg) - gaus->mean()) < 0.10);
    }

    {
        std::vector<MarginalConfig> cfgs = {
            GaussianMarginalCfg{3.0, 0.5},
            FlatMarginalCfg{1.0, 5.0}
        };

        std::vector<MarginalType> kinds = {
            MarginalType::GAUSSIAN,
            MarginalType::FLAT
        };

        for (std::size_t i = 0; i < cfgs.size(); ++i) {
            auto dist = MarginalFactory::create(kinds[i], cfgs[i], static_cast<unsigned int>(100 + i));
            assert(dist != nullptr);

            const double p = 0.3;
            const double q = dist->ppf(p);
            const double p2 = dist->cdf(q);

            assert(approx(p2, p, 1e-8));

            Vector s = dist->rvs(10);
            assert(s.size() == 10);
        }
    }

    std::cout << "\nMarginal INTEGRATION tests passed!\n";
    return 0;
}