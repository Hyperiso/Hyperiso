#include "Fit.h"

FitResult MLEstimator::fit(const Vector& p0) const {
    std::size_t p_dim = p0.size();

    auto f = [this, p_dim] (Vector p_and_eta) -> double {
        auto p_first = p_and_eta.begin();
        auto p_last = p_first + p_dim;
        std::vector<double> p (p_first, p_last);
        std::vector<double> eta (p_last, p_and_eta.end());
        double v = this->like_.nll(p, eta);
        return std::isfinite(v) ? v : 1e300;
    };

    MinimizationContext min_ctx;
    // min_ctx.bfgs_max_iter = max_iter;
    // min_ctx.final_tol = tol;
    min_ctx.final_tol = 1e-8;
    min_ctx.switch_tol = 1e-3;
    min_ctx.simplex_initial_step_size = 0.2;
    min_ctx.simplex_max_iter = 1000;

    Vector start = p0;
    Vector eta_start = like_.get_eta_central_values();
    start.insert(start.end(), eta_start.begin(), eta_start.end());

    Vector p_scales (p0);
    std::transform(p_scales.begin(), p_scales.end(), p_scales.begin(), [] (double x) { return std::abs(x); });
    Vector eta_scales = like_.get_eta_standard_devs();
    Vector scales;
    scales.insert(scales.end(), p_scales.begin(), p_scales.end());
    scales.insert(scales.end(), eta_scales.begin(), eta_scales.end());

    LOG_INFO("First minimization in fit");
    MinimizationResult min_res = minimize_combined(f, start, scales, min_ctx);

    for (auto v : p0) std::cout << v << " ";
    for (auto v : eta_start) std::cout << v << " ";
    std::cout << std::endl;
    for (auto v : min_res.argmin) std::cout << v << " ";
    std::cout << std::endl;

    auto p_and_eta_hat = min_res.argmin;
    auto first = p_and_eta_hat.begin();
    auto last = first + p_dim;
    std::vector<double> p_hat (first, last);
    std::vector<double> eta_hat (last, p_and_eta_hat.end());

    // std::cout << "m_t =" << p_hat[0] << std::endl;
    // std::cout << "f_Bd =" << eta_hat[0] << std::endl;

    // double f_start = f(start);
    // double f_end   = min_res.min;
    // std::cout << "[DBG] f_start=" << f_start << " f_end=" << f_end << std::endl;

    auto f_pr = [this, p_dim] (Vector p) -> double {
        double v = this->like_.nll_profiled(p);
        return std::isfinite(v) ? v : 1e300;
    };

    // double x_min = min_res.argmin[0] - 0.1 * p0[0];
    // double y_min = min_res.argmin[1] - 0.1 * p0[1];
    // size_t n_x = 10;
    // size_t n_y = 10;
    // double step_x = 0.2 * p0[0] / n_x;
    // double step_y = 0.2 * p0[1] / n_y;

    // std::ofstream fs;
    // fs.open("profiled_likelihood.csv");
    // fs << "x,y,f\n";

    // LOG_INFO("Starting loop");
    // double x = x_min;
    // double y = y_min;
    // for (size_t i = 0; i < n_x; i++) {
    //     y = y_min;
    //     for (size_t j = 0; j < n_y; j++){
    //         LOG_INFO("x =", x, ", y =", y);
    //         fs << x << "," << y << "," << f_pr({x, y}) << '\n';
    //         y += step_y;
    //     }
    //     x += step_x;
    // }

    // fs.close();

    // RealMatrix H = hessian(f_pr, p_hat);

    // std::cout << H << std::endl;

    RealMatrix H_inv = inverse_hessian(f_pr, p_hat, p_scales);

    // debug
    // std::cout << H_inv << std::endl;

    double delta_chi2 = gsl_cdf_chisq_Pinv(0.682689492137, p_dim);
    // std::cout << delta_chi2 << std::endl;
    RealMatrix cov = 2 * delta_chi2 * H_inv;
    std::cout << cov << std::endl;

    Vector p_std (p_dim, 0.0);
    for (size_t i = 0; i < p_dim; i++) {
        p_std[i] = std::sqrt(cov.at(i, i));
    }

    RealMatrix corr (p_dim, p_dim);
    for (size_t i = 0; i < p_dim; i++) {
        for (size_t j = 0; j < p_dim; j++) {
            corr.at(i, j) = cov.at(i, j) / (p_std[i] * p_std[j]);
        }
    }

    std::cout << corr << std::endl;
        
    return FitResult {p_hat, eta_hat, p_std, corr, min_res.min}; 
}

std::set<std::vector<std::pair<double, double>>> MLEstimator::contour(double z, std::array<double, 4> bounds, double ell_hat) const {
    // TODO : Make dimension check available
    // if (this->dim != 2)
        // throw std::runtime_error("Contour routines are only available for 2-dimensional fitting systems.")

    auto f = [this, ell_hat] (Vector p) {
        return this->wilks_T(p, ell_hat);
    };

    MarchingSquaresExtractor extractor(f, bounds, 4);
    double delta_T = gsl_cdf_chisq_Pinv(gsl_sf_erf(z / RT2), 2);
    return extractor.find_iso_contour(delta_T);
}
