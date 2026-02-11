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
    min_ctx.max_iter = max_iter;
    min_ctx.tol = tol;

    Vector start = p0;
    Vector eta_start = like_.get_eta_central_values();
    start.insert(start.end(), eta_start.begin(), eta_start.end());

    auto p = p0;
    auto eta = like_.get_eta_central_values();
    double a = like_.nll(p, eta);
    double b = like_.nll(p, eta);
    std::cout << std::setprecision(17)
            << "nll repeat diff = " << (a-b) << "\n";

    // {
    //     std::vector<double> p(p0.begin(), p0.end());
    //     std::vector<double> eta = like_.get_eta_central_values();

    //     double n1 = like_.nll(p, eta);
    //     double n2 = like_.nll(p, eta);

    //     std::cout << std::setprecision(17)
    //             << "[DBG] fixed-point nll1=" << n1
    //             << " nll2=" << n2
    //             << " diff=" << (n1 - n2) << "\n";
    // }

    // MinimizationResult min_res = minimize(f, start, min_ctx);
    MinimizationResult min_res = minimize_BFGS(f, start, min_ctx);

    auto p_and_eta_hat = min_res.argmin;
    auto first = p_and_eta_hat.begin();
    auto last = first + p_dim;
    std::vector<double> p_hat (first, last);
    std::vector<double> eta_hat (last, p_and_eta_hat.end());

    double f_start = f(start);
    double f_end   = min_res.min;
    std::cout << "[DBG] f_start=" << f_start << " f_end=" << f_end << std::endl;

    auto f_pr = [this, p_dim] (Vector p) -> double {
        double v = this->like_.nll_profiled(p);
        return std::isfinite(v) ? v : 1e300;
    };

    RealMatrix H = hessian(f_pr, p_hat);

    std::cout << H << std::endl;

    RealMatrix H_inv = inverse_hessian(f_pr, p_hat);

    // debug
    std::cout << H_inv << std::endl;

    double delta_chi2 = gsl_cdf_chisq_Pinv(0.682689492137, p_dim);
    std::cout << delta_chi2 << std::endl;
    RealMatrix cov = 2 * delta_chi2 * H_inv;

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
