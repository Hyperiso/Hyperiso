#include "Include.h"
#include "IProfilingStrategy.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

static bool approx(double a, double b, double eps = 1e-12) {
    return std::fabs(a - b) <= eps;
}

static FitResult make_fit_result() {
    FitResult fr;
    fr.p_hat = {1.0, 2.0, 3.0};
    fr.eta_hat = {0.1, 0.2};
    fr.p_hat_std = {0.5, 0.6, 0.7};
    fr.p_hat_correlations = RealMatrix(3, 3);
    fr.ell_hat = 42.0;
    return fr;
}

int main() {
    std::cout << "== IProfilingStrategy UNIT ==\n";

    {
        FitResult fr = make_fit_result();
        SliceProfilingStrategy s(0, 2, fr);

        assert(s.get_x_id() == 0);
        assert(s.get_y_id() == 2);

        auto warm = s.init_warm_start();
        assert(warm.size() == 2);
        assert(approx(warm.at(3), 0.1));
        assert(approx(warm.at(4), 0.2));

        std::map<std::size_t, double> current{{3, 0.7}, {4, 0.8}};
        ProfileRequest pr = s.build_request(10.0, 20.0, current);

        assert(pr.start.size() == 5);
        assert(pr.fixed_params.size() == 3);
        assert(pr.free_params.size() == 2);

        assert(approx(pr.fixed_params.at(0), 10.0));
        assert(approx(pr.fixed_params.at(1), 2.0));
        assert(approx(pr.fixed_params.at(2), 20.0));

        assert((pr.free_params == std::vector<std::size_t>{3, 4}));

        assert(approx(pr.start[0], 10.0));
        assert(approx(pr.start[1], 2.0));
        assert(approx(pr.start[2], 20.0));
        assert(approx(pr.start[3], 0.7));
        assert(approx(pr.start[4], 0.8));
    }

    {
        FitResult fr = make_fit_result();
        SliceProfilingStrategy s(0, 2, fr);

        bool threw = false;
        try {
            (void)s.build_request(1.0, 2.0, std::map<std::size_t, double>{{3, 0.0}});
        } catch (const std::runtime_error&) {
            threw = true;
        }
        assert(threw);
    }

    {
        FitResult fr = make_fit_result();
        ProjectionProfilingStrategy s(0, 2, fr);

        assert(s.get_x_id() == 0);
        assert(s.get_y_id() == 2);

        auto warm = s.init_warm_start();
        assert(warm.size() == 3);
        assert(approx(warm.at(1), 2.0));
        assert(approx(warm.at(3), 0.1));
        assert(approx(warm.at(4), 0.2));
        assert(!warm.contains(0));
        assert(!warm.contains(2));

        std::map<std::size_t, double> current{{1, 11.0}, {3, 33.0}, {4, 44.0}};
        ProfileRequest pr = s.build_request(10.0, 20.0, current);

        assert(pr.start.size() == 5);
        assert(pr.fixed_params.size() == 2);
        assert(pr.free_params.size() == 3);

        assert(approx(pr.fixed_params.at(0), 10.0));
        assert(approx(pr.fixed_params.at(2), 20.0));
        assert((pr.free_params == std::vector<std::size_t>{1, 3, 4}));

        assert(approx(pr.start[0], 10.0));
        assert(approx(pr.start[1], 11.0));
        assert(approx(pr.start[2], 20.0));
        assert(approx(pr.start[3], 33.0));
        assert(approx(pr.start[4], 44.0));
    }

    {
        FitResult fr = make_fit_result();
        ProjectionProfilingStrategy s(0, 2, fr);

        bool threw = false;
        try {
            (void)s.build_request(1.0, 2.0, std::map<std::size_t, double>{{1, 0.0}, {3, 0.0}});
        } catch (const std::runtime_error&) {
            threw = true;
        }
        assert(threw);
    }

    std::cout << "\nIProfilingStrategy UNIT tests passed!\n";
    return 0;
}
