#include "Include.h"
#include "Profiler.h"

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

static fit_app::ParameterDefinition make_def(const std::string& name,
                                             double value,
                                             double step_hint = 1.0) {
    fit_app::ParameterDefinition d;
    d.name = name;
    d.value = value;
    d.step_hint = step_hint;
    return d;
}

class SimpleLikelihood : public ILikelihood {
public:
    std::vector<fit_app::ParameterDefinition> defs = {
        make_def("x", 1.0),
        make_def("y", 2.0),
        make_def("z", 3.0)
    };

    mutable std::vector<double> last_theta;
    mutable int nll_calls = 0;

    double nll(const std::vector<double>& theta) const override {
        nll_calls++;
        last_theta = theta;
        double out = 0.0;
        for (double x : theta) out += x * x;
        return out;
    }

    std::vector<fit_app::ParameterDefinition> get_param_defs() const override {
        return defs;
    }

    std::size_t dim() const override {
        return defs.size();
    }
};

int main() {
    std::cout << "== Profiler UNIT ==\n";

    {
        auto like = std::make_shared<SimpleLikelihood>();
        Profiler profiler(nullptr, ProfilerMode::MINUIT);

        ProfileRequest pr;
        pr.free_params = {};
        pr.fixed_params = {{1, 3.0}};
        pr.start = {1.0, 2.0, 4.0};

        ProfileResult r = profiler.profile(like, pr);

        assert(r.converged);
        assert(approx(r.nll_hat, 1.0 * 1.0 + 3.0 * 3.0 + 4.0 * 4.0));
        assert(like->nll_calls == 1);
        assert((like->last_theta == std::vector<double>{1.0, 3.0, 4.0}));
        assert(r.theta_hat.size() == 3);
        assert(approx(r.theta_hat.at(0), 1.0));
        assert(approx(r.theta_hat.at(1), 3.0));
        assert(approx(r.theta_hat.at(2), 4.0));
    }

    {
        auto like = std::make_shared<SimpleLikelihood>();
        Profiler profiler(nullptr, ProfilerMode::LAPLACE_NUISANCE);

        ProfileRequest pr;
        pr.free_params = {};
        pr.fixed_params = {{0, 2.0}, {2, -1.0}};
        pr.start = {}; // direct path fills a zero vector before applying fixed values.

        ProfileResult r = profiler.profile(like, pr);
        assert(r.converged);
        assert(approx(r.nll_hat, 2.0 * 2.0 + 0.0 * 0.0 + (-1.0) * (-1.0)));
        assert((like->last_theta == std::vector<double>{2.0, 0.0, -1.0}));
    }

    {
        auto like = std::make_shared<SimpleLikelihood>();
        Profiler profiler(nullptr);

        ProfileRequest pr;
        pr.free_params = {};
        pr.fixed_params = {};
        pr.start = {1.0, 2.0};

        bool threw = false;
        try {
            (void)profiler.profile(like, pr);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        assert(threw);
    }

    {
        auto like = std::make_shared<SimpleLikelihood>();
        Profiler profiler(nullptr);

        ProfileRequest pr;
        pr.free_params = {};
        pr.fixed_params = {{99, 1.0}};
        pr.start = {1.0, 2.0, 3.0};

        bool threw = false;
        try {
            (void)profiler.profile(like, pr);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        assert(threw);
    }

    std::cout << "\nProfiler UNIT tests passed!\n";
    return 0;
}
