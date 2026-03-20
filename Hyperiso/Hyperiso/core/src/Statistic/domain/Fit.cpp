#include "Fit.h"
// #include "Fit.h"

// FitResult MLEstimator::fit(const std::vector<double>& p0) const {
//     std::size_t p_dim = p0.size();

//     auto f = [this, p_dim] (std::vector<double>p_and_eta) -> double {
//         auto p_first = p_and_eta.begin();
//         auto p_last = p_first + p_dim;
//         std::vector<double> p (p_first, p_last);
//         std::vector<double> eta (p_last, p_and_eta.end());
//         double v = this->like_.nll(p, eta);
//         return std::isfinite(v) ? v : 1e300;
//     };

//     MinimizationContext min_ctx;
//     // min_ctx.bfgs_max_iter = max_iter;
//     // min_ctx.final_tol = tol;
//     min_ctx.final_tol = 1e-8;
//     min_ctx.switch_tol = 1e-3;
//     min_ctx.simplex_initial_step_size = 0.1;
//     min_ctx.simplex_max_iter = this->max_iter;

//     std::vector<double>start = p0;
//     std::vector<double>eta_start = like_.get_eta_central_values();
//     start.insert(start.end(), eta_start.begin(), eta_start.end());

//     std::vector<double>p_scales (p0);
//     std::transform(p_scales.begin(), p_scales.end(), p_scales.begin(), [] (double x) { return std::abs(x); });
//     std::vector<double>eta_scales = like_.get_eta_standard_devs();
//     std::vector<double>scales;
//     scales.insert(scales.end(), p_scales.begin(), p_scales.end());
//     scales.insert(scales.end(), eta_scales.begin(), eta_scales.end());
//     std::cout << std::setprecision(17);
//     LOG_INFO("First minimization in fit");
//     std::cout << "start="; for (auto v : start) std::cout << v << " ";
//     std::cout << "\nscales="; for (auto v : scales) std::cout << v << " ";
//     std::cout << "\n";
//     MinimizationResult min_res = minimize_combined(f, start, scales, min_ctx);

//     for (auto v : p0) std::cout << v << " ";
//     for (auto v : eta_start) std::cout << v << " ";
//     std::cout << std::endl;
//     for (auto v : min_res.argmin) std::cout << v << " ";
//     std::cout << std::endl;

//     auto p_and_eta_hat = min_res.argmin;
//     auto first = p_and_eta_hat.begin();
//     auto last = first + p_dim;
//     std::vector<double> p_hat (first, last);
//     std::vector<double> eta_hat (last, p_and_eta_hat.end());

//     // std::cout << "m_t =" << p_hat[0] << std::endl;
//     // std::cout << "f_Bd =" << eta_hat[0] << std::endl;

//     // double f_start = f(start);
//     // double f_end   = min_res.min;
//     // std::cout << "[DBG] f_start=" << f_start << " f_end=" << f_end << std::endl;

//     auto f_pr = [this, p_dim] (std::vector<double>p) -> double {
//         double v = this->like_.nll_profiled(p);
//         return std::isfinite(v) ? v : 1e300;
//     };

//     // double x_min = min_res.argmin[0] - 0.1 * p0[0];
//     // double y_min = min_res.argmin[1] - 0.1 * p0[1];
//     // size_t n_x = 10;
//     // size_t n_y = 10;
//     // double step_x = 0.2 * p0[0] / n_x;
//     // double step_y = 0.2 * p0[1] / n_y;

//     // std::ofstream fs;
//     // fs.open("profiled_likelihood.csv");
//     // fs << "x,y,f\n";

//     // LOG_INFO("Starting loop");
//     // double x = x_min;
//     // double y = y_min;
//     // for (size_t i = 0; i < n_x; i++) {
//     //     y = y_min;
//     //     for (size_t j = 0; j < n_y; j++){
//     //         LOG_INFO("x =", x, ", y =", y);
//     //         fs << x << "," << y << "," << f_pr({x, y}) << '\n';
//     //         y += step_y;
//     //     }
//     //     x += step_x;
//     // }

//     // fs.close();

//     // std::cout << "Hessian : " << std::endl;

//     // RealMatrix H = hessian(f_pr, p_hat, p_scales);

//     // std::cout << H << std::endl;

//     RealMatrix H_inv = inverse_hessian(f_pr, p_hat, p_scales);

//     std::cout << "H_inv : " << std::endl;
//     // debug
//     std::cout << H_inv << std::endl;

//     // std::cout << "product : " << std::endl;

//     // std::cout << H_inv * H << std::endl;

//     double delta_chi2 = gsl_cdf_chisq_Pinv(0.682689492137, p_dim);
//     // std::cout << delta_chi2 << std::endl;
//     RealMatrix cov = 2 * delta_chi2 * H_inv;
//     std::cout << "Covariance matrix : " << std::endl;
//     std::cout << cov << std::endl;

//     std::vector<double>p_std (p_dim, 0.0);
//     for (size_t i = 0; i < p_dim; i++) {
//         p_std[i] = std::sqrt(cov.at(i, i));
//     }

//     RealMatrix corr (p_dim, p_dim);
//     for (size_t i = 0; i < p_dim; i++) {
//         for (size_t j = 0; j < p_dim; j++) {
//             corr.at(i, j) = cov.at(i, j) / (p_std[i] * p_std[j]);
//         }
//     }

//     std::cout << corr << std::endl;
        
//     double nll_joint = like_.nll(p_hat, eta_hat);
// double nll_prof  = like_.nll_profiled(p_hat);

// // std::cout << std::setprecision(17);
// // std::cout << "[TEST C] nll(p_hat, eta_hat) = " << nll_joint << "\n";
// // std::cout << "[TEST C] nll_profiled(p_hat) = " << nll_prof << "\n";
// // std::cout << "[TEST C] diff = " << (nll_joint - nll_prof) << "\n";

//     return FitResult {p_hat, eta_hat, p_std, corr, min_res.min}; 
// }

// std::set<std::vector<std::pair<double, double>>> MLEstimator::contour(double z, std::array<double, 4> bounds, double ell_hat) const {
//     // TODO : Make dimension check available
//     // if (this->dim != 2)
//         // throw std::runtime_error("Contour routines are only available for 2-dimensional fitting systems.")

//     auto f = [this, ell_hat] (std::vector<double>p) {
//         return this->wilks_T(p, ell_hat);
//     };

//     MarchingSquaresExtractor extractor(f, bounds, 4);
//     double delta_T = gsl_cdf_chisq_Pinv(gsl_sf_erf(z / RT2), 2);
//     return extractor.find_iso_contour(delta_T);
// }

MLFitter::MLFitter(std::shared_ptr<LikelihoodContext> ctx, const ModelFn &model) {
    this->like_ = std::make_shared<BaseLikelihood>(model, std::move(ctx), ctx->fp_defs.size());
}

FitResult MLFitter::maximum_likelihood_fit(const std::vector<double>&p0) {
    std::shared_ptr<Profiler> profiler = std::make_shared<Profiler>(fit_app::make_minuit_backend());
    ProfileRequest pr_model;
    for (size_t i = p0.size(); i < like_->dim(); i++) {
        pr_model.free_params.push_back(i);
        pr_model.start.push_back(this->like_->get_param_defs().at(i).value);
    }

    auto f = fit_app::LambdaObjectiveFunction(
        [this, profiler, pr_model] (const std::vector<double>& p) { 
            ProfileRequest pr;
            pr.free_params = pr_model.free_params;
            pr.start = pr_model.start;

            for (size_t i = 0; i < p.size(); i++)
                pr.fixed_params[i] = p[i];
            
            return profiler->profile(this->like_, pr).nll_hat;
        },
        0.5
    );
    
    fit_app::FitOptions opt;
    opt.run_hesse = false;
    opt.verbose = false;

    std::vector<fit_app::ParameterDefinition> theta_0 = std::vector(
        this->like_->get_param_defs().begin(),
        this->like_->get_param_defs().begin() + p0.size()
    );

    for (size_t i = 0; i < p0.size(); i++)
        theta_0[i].value = p0[i];

    std::unique_ptr<fit_app::IFitBackend> minimizer = fit_app::make_minuit_backend();
    fit_app::BackendFitResult master_fit_res = minimizer->minimize(f, theta_0, opt);

    if (!master_fit_res.diagnostics.ok)
        LOG_WARN("Initial ML fit failed to converge. Fit result is probably wrong.");

    FitResult fr;
    fr.ell_hat = master_fit_res.diagnostics.fmin;
    fr.p_hat = master_fit_res.values;
    
    RealMatrix cov = master_fit_res.covariance;
    for (size_t i = 0; i < cov.rows(); i++)
        fr.p_hat_std.emplace_back(std::sqrt(cov.at(i, i)));

    fr.p_hat_correlations = RealMatrix(cov);
    for (size_t i = 0; i < cov.rows(); i++) {
        for (size_t j = 0; j < cov.cols(); j++) {
            fr.p_hat_correlations.at(i, j) /= (fr.p_hat_std.at(i) * fr.p_hat_std.at(j));
        }
    }
    
    this->master_fit_success = master_fit_res.diagnostics.ok;
    this->master_fit_result = fr;
    return fr;
}

Contour MLFitter::contour(std::size_t x_id, std::size_t y_id, double z, std::array<double, 4> bounds, ProfilingMethod method = ProfilingMethod::SLICE) const {
    if (this->master_fit_success)
        LOG_ERROR("InvalidState", "ML fi must have converged before contour computation is available.");

    // TODO : Make resolution and contouring algorithms available as options to the user.
    ContourConfig cc;
    cc.fr = this->master_fit_result;
    cc.x_id = x_id;
    cc.y_id = y_id;
    cc.primary_contour_method = ContourAlgorithm::MINUIT;
    cc.fallback_contour_method = ContourAlgorithm::AMS;
    cc.profiling_method = method;

    ContourEngine ce(this->like_, cc);
    Contour cl = ce.compute_contour(z * z / 2, bounds, 200);

    return cl;
}
