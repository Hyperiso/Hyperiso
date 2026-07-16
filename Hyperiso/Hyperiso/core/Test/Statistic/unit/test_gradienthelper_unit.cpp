#include "Include.h"
#include "GradientHelper.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

static bool approx(double a, double b, double eps = 1e-6) {
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

static fit_app::ParameterDefinition make_limited_def(
    const std::string& name,
    double value,
    double step_hint,
    double lo,
    double hi
) {
    auto d = make_def(name, value, step_hint);
    d.limits = std::pair<double, double>{lo, hi};
    return d;
}

class QuadraticProfileableLikelihood : public IProfileableLikelihood {
public:
    double y = 1.5;
    std::vector<double> a = {2.0, -3.0};
    std::vector<fit_app::ParameterDefinition> defs;

    QuadraticProfileableLikelihood() {
        defs = {
            make_def("p0", 1.0, 0.1),
            make_def("eta0", 0.0, 1.0),
            make_def("eta1", 0.0, 1.0)
        };
    }

    double nll(const std::vector<double>& theta) const override {
        return nll_from_split(
            std::vector<double>{theta[0]},
            std::vector<double>{theta[1], theta[2]}
        );
    }

    std::vector<fit_app::ParameterDefinition> get_param_defs() const override {
        return defs;
    }

    std::size_t dim() const override { return 3; }
    std::size_t p_dimension() const override { return 1; }
    std::size_t eta_dimension() const override { return 2; }
    std::vector<double> central_p() const override { return {defs[0].value}; }
    std::vector<double> central_eta() const override { return {defs[1].value, defs[2].value}; }

    std::vector<double> predict(
        const std::vector<double>& p,
        const std::vector<double>& eta
    ) const override {
        return {p[0] + a[0] * eta[0] + a[1] * eta[1]};
    }

    std::vector<double> residuals(
        const std::vector<double>& p,
        const std::vector<double>& eta
    ) const override {
        return {predict(p, eta)[0] - y};
    }

    double nll_from_split(
        const std::vector<double>& p,
        const std::vector<double>& eta
    ) const override {
        const double r = residuals(p, eta)[0];
        return 0.5 * r * r + 0.5 * (eta[0] * eta[0] + eta[1] * eta[1]);
    }

    RealMatrix observable_curvature(const std::vector<double>&) const override {
        RealMatrix W(1, 1);
        W.at(0, 0) = 1.0;
        return W;
    }

    RealMatrix nuisance_curvature(const std::vector<double>&) const override {
        RealMatrix W(2, 2);
        W.at(0, 0) = 1.0;
        W.at(1, 1) = 1.0;
        W.at(0, 1) = 0.0;
        W.at(1, 0) = 0.0;
        return W;
    }
};

int main() {
    std::cout << "== GradientHelper UNIT ==\n";

    {
        assert(approx(fd_step(0.0, 0.0), 1e-5, 1e-12));
        assert(approx(fd_step(100.0, 0.0), 1e-3, 1e-12));
        assert(approx(fd_step(100.0, 0.2), 1e-3, 1e-12));
        assert(approx(fd_step(100.0, 0.01), 1e-4, 1e-12));
    }

    {
        auto d = make_limited_def("eta", 0.0, 1.0, -0.01, 0.01);
        const double h = eta_fd_step_with_limits(d, 0.0);
        assert(h > 0.0);
        assert(h <= 0.0045 + 1e-14);

        // At an exactly closed interval boundary [0,0], the current helper keeps
        // the generic finite-difference step instead of throwing.
        auto degenerate = make_limited_def("eta_degenerate", 0.0, 1.0, 0.0, 0.0);
        const double h_degenerate = eta_fd_step_with_limits(degenerate, 0.0);
        assert(std::isfinite(h_degenerate));
        assert(h_degenerate > 0.0);

        // A genuine invalid step is obtained when both the evaluation point and
        // the step hint cannot produce a finite fallback step.
        auto invalid = make_def("eta_invalid", 0.0, 0.0);
        bool threw = false;
        try {
            (void)eta_fd_step_with_limits(invalid, std::numeric_limits<double>::infinity());
        } catch (const std::runtime_error&) {
            threw = true;
        }
        assert(threw);
    }

    {
        assert((all_eta_indices(4) == std::vector<std::size_t>{0, 1, 2, 3}));
        assert((eta_index_complement(5, std::vector<std::size_t>{1, 3}) ==
                std::vector<std::size_t>{0, 2, 4}));
        assert(eta_index_contains(std::vector<std::size_t>{2, 5}, 5));
        assert(!eta_index_contains(std::vector<std::size_t>{2, 5}, 4));
        assert(eta_index_list_string(std::vector<std::size_t>{1, 4, 7}) == "{1,4,7}");
    }

    {
        RealMatrix M(3, 3);
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                M.at(i, j) = 10.0 * static_cast<double>(i) + static_cast<double>(j);
            }
        }

        RealMatrix S = principal_submatrix_by_indices(M, std::vector<std::size_t>{0, 2});
        assert(S.rows() == 2 && S.cols() == 2);
        assert(approx(S.at(0, 0), 0.0));
        assert(approx(S.at(0, 1), 2.0));
        assert(approx(S.at(1, 0), 20.0));
        assert(approx(S.at(1, 1), 22.0));
    }

    {
        RealMatrix M(2, 2);
        M.at(0, 0) = 1.0;
        M.at(0, 1) = 2.0;
        M.at(1, 0) = 3.0;
        M.at(1, 1) = 4.0;

        auto mv = matvec(M, std::vector<double>{5.0, 6.0});
        assert(mv.size() == 2);
        assert(approx(mv[0], 17.0));
        assert(approx(mv[1], 39.0));
        assert(approx(dot(std::vector<double>{1.0, 2.0}, std::vector<double>{3.0, 4.0}), 11.0));

        bool threw_mv = false;
        try {
            (void)matvec(M, std::vector<double>{1.0});
        } catch (const std::runtime_error&) {
            threw_mv = true;
        }
        assert(threw_mv);

        bool threw_dot = false;
        try {
            (void)dot(std::vector<double>{1.0}, std::vector<double>{1.0, 2.0});
        } catch (const std::runtime_error&) {
            threw_dot = true;
        }
        assert(threw_dot);
    }

    {
        RealMatrix A(2, 2);
        A.at(0, 0) = 2.0;
        A.at(0, 1) = 1.0;
        A.at(1, 0) = 3.0;
        A.at(1, 1) = 4.0;

        RealMatrix S = force_symmetric_checked(A, "A");
        assert(approx(S.at(0, 0), 2.0));
        assert(approx(S.at(1, 1), 4.0));
        assert(approx(S.at(0, 1), 2.0));
        assert(approx(S.at(1, 0), 2.0));

        RealMatrix bad(2, 2);
        bad.at(0, 0) = 1.0;
        bad.at(0, 1) = std::numeric_limits<double>::infinity();
        bad.at(1, 0) = 0.0;
        bad.at(1, 1) = 1.0;
        bool threw = false;
        try {
            (void)force_symmetric_checked(bad, "bad");
        } catch (const std::runtime_error&) {
            threw = true;
        }
        assert(threw);
    }

    {
        RealMatrix H(2, 2);
        H.at(0, 0) = 2.0;
        H.at(0, 1) = 0.0;
        H.at(1, 0) = 0.0;
        H.at(1, 1) = 3.0;

        RealMatrix R = regularize_spd_local(H, 1e-10, "H");
        assert(R.rows() == 2 && R.cols() == 2);
        assert(R.at(0, 0) > 0.0);
        assert(R.at(1, 1) > 0.0);
        assert(approx(R.at(0, 1), R.at(1, 0), 1e-12));
    }

    {
        QuadraticProfileableLikelihood like;
        const std::vector<double> p{1.0};
        const std::vector<double> eta{0.2, -0.4};

        const EtaDerivatives der = compute_eta_derivatives_subset(
            like,
            p,
            eta,
            std::vector<std::size_t>{0, 1}
        );

        assert(der.g_eta.size() == 2);
        assert(der.J_eta.rows() == 1 && der.J_eta.cols() == 2);

        const double r = p[0] + 2.0 * eta[0] - 3.0 * eta[1] - like.y;
        assert(approx(der.J_eta.at(0, 0), 2.0, 1e-5));
        assert(approx(der.J_eta.at(0, 1), -3.0, 1e-5));
        assert(approx(der.g_eta[0], 2.0 * r + eta[0], 1e-4));
        assert(approx(der.g_eta[1], -3.0 * r + eta[1], 1e-4));

        auto g = eta_gradient_nll(like, p, eta);
        assert(g.size() == 2);
        assert(approx(g[0], der.g_eta[0], 1e-6));
        assert(approx(g[1], der.g_eta[1], 1e-6));
    }

    std::cout << "\nGradientHelper UNIT tests passed!\n";
    return 0;
}
