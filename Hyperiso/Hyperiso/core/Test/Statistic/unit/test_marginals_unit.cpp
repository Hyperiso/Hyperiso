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
#include <map>
#include <stdexcept>
#include <string>

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

class StubStatParameterProxy final : public IStatParameterProxy {
public:
    mutable std::size_t scalar_reads {0};

    std::shared_ptr<Parameter> get_param(const ParamId&) const override {
        throw std::runtime_error("Unexpected get_param(ParamId) call");
    }

    std::shared_ptr<Parameter> get_param(
        const std::string&,
        const LhaID&
    ) const override {
        throw std::runtime_error("Unexpected get_param(block, id) call");
    }

    scalar_t operator()(const ParamId&, DataType data_type) const override {
        ++scalar_reads;
        if (data_type == DataType::VALUE) {
            return scalar_t{2.0};
        }
        if (data_type == DataType::STD_COMBINED) {
            return scalar_t{0.75};
        }
        throw std::runtime_error("Unexpected parameter data type");
    }

    std::map<ExperimentObs, double> operator()(
        const ObservableId&,
        DataType
    ) const override {
        throw std::runtime_error("Unexpected observable proxy call");
    }

    std::map<ExperimentObs, double> operator()(
        const BinnedObservableId&,
        DataType
    ) const override {
        throw std::runtime_error("Unexpected binned-observable proxy call");
    }

    scalar_t operator()(
        const std::string&,
        const LhaID&,
        DataType
    ) const override {
        throw std::runtime_error("Unexpected block proxy call");
    }

    std::map<ExperimentObs, std::shared_ptr<Parameter>> get_obs_param(
        const BinnedObservableId&
    ) const override {
        throw std::runtime_error("Unexpected get_obs_param call");
    }
};

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

        const double flat_logpdf = std::log(1.0 / 6.0);
        assert(approx(f.logpdf(0.0), flat_logpdf, 1e-12));
        assert(approx(f.logpdf(-2.0), flat_logpdf, 1e-12));
        assert(approx(f.logpdf(4.0), flat_logpdf, 1e-12));
        assert(std::isinf(f.logpdf(100.0)) && f.logpdf(100.0) < 0.0);

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
        assert(approx(
            dist->logpdf(-std::sqrt(3.0)),
            -std::log(2.0 * std::sqrt(3.0)),
            1e-12
        ));
    }


    {
        auto parameter_proxy = std::make_shared<StubStatParameterProxy>();
        MarginalConfigFactory config_factory(parameter_proxy);

        const ParamId runtime_pid{ParameterType::DECAY, BlockName{"B_Xs"}, LhaID{7}};

        MarginalConfig gaussian_cfg = config_factory.create(
            runtime_pid,
            MarginalType::GAUSSIAN
        );
        const auto& gaussian = std::get<GaussianMarginalCfg>(gaussian_cfg);
        assert(approx(gaussian.mu, 2.0, 1e-12));
        assert(approx(gaussian.sigma, 0.75, 1e-12));
        assert(parameter_proxy->scalar_reads == 2);

        const ParamId config_pid{BlockName{"B_Xs"}, LhaID{7}};
        const NuisanceSpec spec{config_pid, {0.9, 4.0}, MarginalType::FLAT};

        const std::size_t reads_before_explicit_bounds =
            parameter_proxy->scalar_reads;
        MarginalConfig flat_cfg_variant = config_factory.create(
            runtime_pid,
            MarginalType::FLAT,
            spec
        );
        const auto& flat_cfg = std::get<FlatMarginalCfg>(flat_cfg_variant);

        assert(approx(flat_cfg.a, 0.9, 1e-12));
        assert(approx(flat_cfg.b, 4.0, 1e-12));
        assert(parameter_proxy->scalar_reads == reads_before_explicit_bounds);
    }

    {
        bool threw = false;
        try {
            MarginalConfigFactory invalid_factory(nullptr);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        assert(threw && "MarginalConfigFactory must reject a null parameter port");
    }

    std::cout << "\nAll marginal UNIT tests passed!\n";
    return 0;
}