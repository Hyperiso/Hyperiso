#include "Include.h"
#include "ChiSquaredLikelihood.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

static bool approx(double a, double b, double eps = 1e-12) {
    return std::fabs(a - b) <= eps;
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

static RealMatrix mat2(double a00, double a01, double a10, double a11) {
    RealMatrix M(2, 2);
    M.at(0, 0) = a00;
    M.at(0, 1) = a01;
    M.at(1, 0) = a10;
    M.at(1, 1) = a11;
    return M;
}

static std::shared_ptr<LikelihoodContext> make_context() {
    auto ctx = std::make_shared<LikelihoodContext>();
    ctx->exp_obs_values = {1.0, -2.0};
    ctx->fp_defs = {make_def("p0", 0.5), make_def("p1", -0.5)};
    ctx->nuis_defs = {};
    return ctx;
}

int main() {
    std::cout << "== ChiSquaredLikelihood UNIT ==\n";

    {
        auto ctx = make_context();
        RealMatrix Cinv = mat2(2.0, 0.5, 0.5, 1.0);

        ModelFn model = [](const std::vector<double>& p,
                           const std::vector<double>& eta) {
            assert(eta.empty());
            return std::vector<double>{p[0] + 1.0, 2.0 * p[1]};
        };

        ChiSquaredLikelihood like(model, ctx, 2, Cinv);

        const std::vector<double> theta{2.0, 1.0};
        const std::vector<double> r{2.0, 4.0}; // pred=(3,2), obs=(1,-2)
        const double q = r[0] * (2.0 * r[0] + 0.5 * r[1])
                       + r[1] * (0.5 * r[0] + 1.0 * r[1]);
        const double expected = 0.5 * q;

        assert(like.dim() == 2);
        assert(like.p_dimension() == 2);
        assert(like.eta_dimension() == 0);
        assert(approx(like.nll(theta), expected));
    }

    {
        auto ctx = make_context();
        RealMatrix Cinv = mat2(3.0, 0.0, 0.0, 4.0);
        ModelFn model = [](const std::vector<double>& p,
                           const std::vector<double>&) {
            return std::vector<double>{p[0], p[1]};
        };

        ChiSquaredLikelihood like(model, ctx, 2, Cinv);

        RealMatrix W = like.observable_curvature(std::vector<double>{0.0, 0.0});
        assert(W.rows() == 2 && W.cols() == 2);
        assert(approx(W.at(0, 0), 3.0));
        assert(approx(W.at(1, 1), 4.0));
        assert(approx(W.at(0, 1), 0.0));
        assert(approx(W.at(1, 0), 0.0));

        RealMatrix Weta = like.nuisance_curvature({});
        assert(Weta.rows() == 0 && Weta.cols() == 0);

        bool threw = false;
        try {
            (void)like.nuisance_curvature(std::vector<double>{1.0});
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        assert(threw);
    }

    {
        auto ctx = make_context();
        RealMatrix Cinv = mat2(1.0, 0.0, 0.0, 1.0);
        ModelFn model = [](const std::vector<double>& p,
                           const std::vector<double>&) {
            return std::vector<double>{p[0], p[1]};
        };

        ChiSquaredLikelihood like(model, ctx, 2, Cinv);

        bool threw = false;
        try {
            (void)like.nll(std::vector<double>{1.0});
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        assert(threw);
    }

    {
        auto ctx = make_context();
        RealMatrix Cinv = mat2(1.0, 0.0, 0.0, 1.0);
        ModelFn wrong_size = [](const std::vector<double>&,
                                const std::vector<double>&) {
            return std::vector<double>{1.0};
        };

        ChiSquaredLikelihood like(wrong_size, ctx, 2, Cinv);
        assert(like.nll(std::vector<double>{0.0, 0.0}) == 1e100);
    }

    {
        auto ctx = make_context();
        RealMatrix Cinv = mat2(1.0, 0.0, 0.0, 1.0);
        ModelFn nonfinite_model = [](const std::vector<double>&,
                                     const std::vector<double>&) {
            return std::vector<double>{std::numeric_limits<double>::infinity(), 0.0};
        };

        ChiSquaredLikelihood like(nonfinite_model, ctx, 2, Cinv);
        assert(like.nll(std::vector<double>{0.0, 0.0}) == 1e100);
    }

    {
        auto ctx = make_context();
        RealMatrix bad_cinv = mat2(1.0, std::numeric_limits<double>::infinity(), 0.0, 1.0);
        ModelFn model = [](const std::vector<double>& p,
                           const std::vector<double>&) {
            return std::vector<double>{p[0], p[1]};
        };

        ChiSquaredLikelihood like(model, ctx, 2, bad_cinv);
        assert(like.nll(std::vector<double>{0.0, 0.0}) == 1e100);
    }

    std::cout << "\nChiSquaredLikelihood UNIT tests passed!\n";
    return 0;
}
