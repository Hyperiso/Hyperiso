// // #include "Fit.h"

// // MLFitter::MLFitter(std::shared_ptr<LikelihoodContext> ctx, const ModelFn& model) {
// //     const std::size_t p_dim = ctx->fp_defs.size();
// //     this->like_ = std::make_shared<BaseLikelihood>(model, ctx, p_dim);
// // }

// // FitResult MLFitter::maximum_likelihood_fit(const std::vector<double>& p0) {
// //     const auto defs = like_->get_param_defs();
// //     const std::size_t p_dim = p0.size();
// //     const std::size_t dim = defs.size();

// //     std::vector<fit_app::ParameterDefinition> theta0 = defs;
// //     for (std::size_t i = 0; i < p_dim; ++i) {
// //         theta0[i].value = p0[i];
// //     }

// //     auto f = fit_app::LambdaObjectiveFunction(
// //         [this](const std::vector<double>& theta) {
// //             return like_->nll(theta);
// //         },
// //         0.5
// //     );

// //     fit_app::FitOptions opt;
// //     opt.run_hesse = true;
// //     opt.verbose = false;

// //     std::unique_ptr<fit_app::IFitBackend> minimizer = fit_app::make_minuit_backend();
// //     fit_app::BackendFitResult res = minimizer->minimize(f, theta0, opt);

// //     FitResult fr;
// //     fr.ell_hat = res.diagnostics.fmin;
// //     fr.p_hat.assign(res.values.begin(), res.values.begin() + p_dim);
// //     fr.eta_hat.assign(res.values.begin() + p_dim, res.values.end());

// //     // bloc p×p de la covariance complète
// //     RealMatrix H = res.covariance.inv();
// //     RealMatrix H_p_p(p_dim, p_dim);
// //     RealMatrix H_p_eta(p_dim, dim - p_dim);
// //     RealMatrix H_eta_eta(dim - p_dim, dim - p_dim);

// //     for (std::size_t i = 0; i < p_dim; ++i) {
// //         for (std::size_t j = 0; j < p_dim; ++j) {
// //             H_p_p.at(i, j) = H.at(i, j);
// //         }
// //     }

// //     for (std::size_t i = 0; i < p_dim; ++i) {
// //         for (std::size_t j = p_dim; j < dim; ++j) {
// //             H_p_eta.at(i, j - p_dim) = H.at(i, j);
// //         }
// //     }

// //     for (std::size_t i = p_dim; i < dim; ++i) {
// //         for (std::size_t j = p_dim; j < dim; ++j) {
// //             H_eta_eta.at(i - p_dim, j - p_dim) = H.at(i, j);
// //         }
// //     }

// //     RealMatrix H_prof = H_p_p - H_p_eta * H_eta_eta.inv() * H_p_eta.transpose();
// //     RealMatrix cov_prof = H_prof.inv();

// //     fr.p_hat_std.resize(p_dim);
// //     fr.p_hat_correlations = RealMatrix(p_dim, p_dim);

// //     for (std::size_t i = 0; i < p_dim; ++i) {
// //         fr.p_hat_std[i] = std::sqrt(std::max(0.0, cov_prof.at(i, i)));
// //     }
// //     for (std::size_t i = 0; i < p_dim; ++i) {
// //         for (std::size_t j = 0; j < p_dim; ++j) {
// //             const double si = fr.p_hat_std[i];
// //             const double sj = fr.p_hat_std[j];
// //             fr.p_hat_correlations.at(i, j) =
// //                 (si > 0.0 && sj > 0.0) ? cov_prof.at(i, j) / (si * sj) : 0.0;
// //         }
// //     }

// //     this->master_fit_success = res.diagnostics.ok;
// //     this->master_fit_result = fr;
// //     return fr;
// // }

// // Contour MLFitter::contour(std::size_t x_id, std::size_t y_id, double z, std::array<double, 4> bounds, ContourOptions options) const {
// //     if (!this->master_fit_success)
// //         LOG_ERROR("InvalidState", "ML fit must have converged before contour computation is available.");

// //     ContourConfig cc;
// //     cc.fr = this->master_fit_result;
// //     cc.x_id = x_id;
// //     cc.y_id = y_id;
// //     cc.primary_contour_method = options.primary_contour_method;
// //     cc.fallback_contour_method = options.fallback_contour_method;
// //     cc.profiling_method = options.profiling_method;

// //     ContourEngine ce(this->like_, cc);
// //     Contour cl = ce.compute_contour(z, bounds, options.resolution);

// //     return cl;
// // }

// #include "Fit.h"

// #include <algorithm>
// #include <iostream>
// #include <limits>
// #include <sstream>
// #include <string>

// namespace {

// void log_fit_diagnostics(const fit_app::FitDiagnostics& d) {
//     std::cout << "[FIT] Minuit diagnostics: ok=" << d.ok
//               << ", has_valid_parameters=" << d.has_valid_parameters
//               << ", has_valid_covar=" << d.has_valid_covar
//               << ", has_posdef_covar=" << d.has_posdef_covar
//               << ", has_accurate_covar=" << d.has_accurate_covar
//               << ", made_posdef=" << d.made_posdef
//               << ", hesse_failed=" << d.hesse_failed
//               << ", reached_call_limit=" << d.reached_call_limit
//               << ", above_max_edm=" << d.above_max_edm
//               << ", fmin=" << d.fmin
//               << ", edm=" << d.edm
//               << ", nfcn=" << d.nfcn
//               << std::endl;

//     if (!d.cov_eigs.empty()) {
//         double min_eig = std::numeric_limits<double>::infinity();
//         double max_eig = -std::numeric_limits<double>::infinity();
//         std::size_t non_pos = 0;

//         for (double eig : d.cov_eigs) {
//             min_eig = std::min(min_eig, eig);
//             max_eig = std::max(max_eig, eig);
//             if (!(eig > 0.0)) ++non_pos;
//         }

//         std::cout << "[FIT] Minuit covariance eigenvalues: min=" << min_eig
//                   << ", max=" << max_eig
//                   << ", non_positive=" << non_pos
//                   << ", cond=" << d.cond_number
//                   << std::endl;
//     }
// }

// void log_matrix_diagnostics(const std::string& label, const RealMatrix& M) {
//     std::cout << "[FIT] Matrix " << label
//               << " shape=" << M.rows() << "x" << M.cols()
//               << std::endl;

//     if (M.rows() == 0 || M.cols() == 0) {
//         std::cout << "[FIT] Matrix " << label << " is empty" << std::endl;
//         return;
//     }

//     if (M.rows() != M.cols()) {
//         std::cout << "[FIT] Matrix " << label << " is not square" << std::endl;
//         return;
//     }

//     try {
//         RealMatrix sym = 0.5 * (M + M.transpose());
//         EigenSystem eig = sym.eig();

//         double min_eig = std::numeric_limits<double>::infinity();
//         double max_eig = -std::numeric_limits<double>::infinity();
//         std::size_t non_pos = 0;
//         std::size_t tiny = 0;

//         for (std::size_t i = 0; i < eig.D.rows(); ++i) {
//             const double lambda = eig.D.at(i, i);
//             min_eig = std::min(min_eig, lambda);
//             max_eig = std::max(max_eig, lambda);
//             if (!(lambda > 0.0)) ++non_pos;
//             if (std::abs(lambda) < 1e-12) ++tiny;
//         }

//         std::cout << "[FIT] Matrix " << label
//                   << " eig_min=" << min_eig
//                   << ", eig_max=" << max_eig
//                   << ", non_positive=" << non_pos
//                   << ", tiny(|eig|<1e-12)=" << tiny
//                   << std::endl;
//     } catch (const std::exception& e) {
//         std::cout << "[FIT] Matrix " << label
//                   << " diagnostics failed: " << e.what()
//                   << std::endl;
//     }
// }

// RealMatrix invert_or_throw(const RealMatrix& M, const std::string& label) {
//     log_matrix_diagnostics(label, M);
//     try {
//         return M.inv();
//     } catch (const std::exception& e) {
//         std::ostringstream oss;
//         oss << "Failed to invert " << label << ": " << e.what();
//         throw std::runtime_error(oss.str());
//     }
// }

// } // namespace

// MLFitter::MLFitter(std::shared_ptr<LikelihoodContext> ctx, const ModelFn& model) {
//     const std::size_t p_dim = ctx->fp_defs.size();
//     this->like_ = std::make_shared<BaseLikelihood>(model, ctx, p_dim);
// }

// FitResult MLFitter::maximum_likelihood_fit(const std::vector<double>& p0) {
//     const auto defs = like_->get_param_defs();
//     const std::size_t p_dim = p0.size();
//     const std::size_t dim = defs.size();

//     std::vector<fit_app::ParameterDefinition> theta0 = defs;
//     for (std::size_t i = 0; i < p_dim; ++i) {
//         theta0[i].value = p0[i];
//     }

//     auto f = fit_app::LambdaObjectiveFunction(
//         [this](const std::vector<double>& theta) {
//             return like_->nll(theta);
//         },
//         0.5
//     );

//     fit_app::FitOptions opt;
//     opt.run_hesse = true;
//     opt.verbose = false;

//     std::unique_ptr<fit_app::IFitBackend> minimizer = fit_app::make_minuit_backend();
//     fit_app::BackendFitResult res = minimizer->minimize(f, theta0, opt);

//     log_fit_diagnostics(res.diagnostics);
//     log_matrix_diagnostics("Minuit covariance", res.covariance);

//     if (!res.diagnostics.has_valid_covar) {
//         throw std::runtime_error(
//             "Minuit returned no valid covariance matrix. "
//             "This usually means the likelihood has flat or redundant directions."
//         );
//     }

//     FitResult fr;
//     fr.ell_hat = res.diagnostics.fmin;
//     fr.p_hat.assign(res.values.begin(), res.values.begin() + p_dim);
//     fr.eta_hat.assign(res.values.begin() + p_dim, res.values.end());

//     // bloc p×p de la covariance complète
//     RealMatrix H = invert_or_throw(res.covariance, "Minuit covariance");
//     RealMatrix H_p_p(p_dim, p_dim);
//     RealMatrix H_p_eta(p_dim, dim - p_dim);
//     RealMatrix H_eta_eta(dim - p_dim, dim - p_dim);

//     for (std::size_t i = 0; i < p_dim; ++i) {
//         for (std::size_t j = 0; j < p_dim; ++j) {
//             H_p_p.at(i, j) = H.at(i, j);
//         }
//     }

//     for (std::size_t i = 0; i < p_dim; ++i) {
//         for (std::size_t j = p_dim; j < dim; ++j) {
//             H_p_eta.at(i, j - p_dim) = H.at(i, j);
//         }
//     }

//     for (std::size_t i = p_dim; i < dim; ++i) {
//         for (std::size_t j = p_dim; j < dim; ++j) {
//             H_eta_eta.at(i - p_dim, j - p_dim) = H.at(i, j);
//         }
//     }

//     RealMatrix cov_prof;
//     if (dim == p_dim) {
//         cov_prof = invert_or_throw(H_p_p, "H_p_p (no nuisance block)");
//     } else {
//         RealMatrix H_eta_eta_inv = invert_or_throw(H_eta_eta, "H_eta_eta");
//         RealMatrix H_prof = H_p_p - H_p_eta * H_eta_eta_inv * H_p_eta.transpose();
//         log_matrix_diagnostics("H_prof", H_prof);
//         cov_prof = invert_or_throw(H_prof, "H_prof");
//     }

//     fr.p_hat_std.resize(p_dim);
//     fr.p_hat_correlations = RealMatrix(p_dim, p_dim);

//     for (std::size_t i = 0; i < p_dim; ++i) {
//         fr.p_hat_std[i] = std::sqrt(std::max(0.0, cov_prof.at(i, i)));
//     }
//     for (std::size_t i = 0; i < p_dim; ++i) {
//         for (std::size_t j = 0; j < p_dim; ++j) {
//             const double si = fr.p_hat_std[i];
//             const double sj = fr.p_hat_std[j];
//             fr.p_hat_correlations.at(i, j) =
//                 (si > 0.0 && sj > 0.0) ? cov_prof.at(i, j) / (si * sj) : 0.0;
//         }
//     }

//     this->master_fit_success = res.diagnostics.ok;
//     this->master_fit_result = fr;
//     return fr;
// }

// Contour MLFitter::contour(std::size_t x_id, std::size_t y_id, double z, std::array<double, 4> bounds, ContourOptions options) const {
//     if (!this->master_fit_success)
//         LOG_ERROR("InvalidState", "ML fit must have converged before contour computation is available.");

//     ContourConfig cc;
//     cc.fr = this->master_fit_result;
//     cc.x_id = x_id;
//     cc.y_id = y_id;
//     cc.primary_contour_method = options.primary_contour_method;
//     cc.fallback_contour_method = options.fallback_contour_method;
//     cc.profiling_method = options.profiling_method;

//     ContourEngine ce(this->like_, cc);
//     Contour cl = ce.compute_contour(z, bounds, options.resolution);

//     return cl;
// }


#include "Fit.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <cmath>
#include <numeric>

namespace {

void log_fit_diagnostics(const fit_app::FitDiagnostics& d) {
    std::cout << "[FIT] Minuit diagnostics: ok=" << d.ok
              << ", has_valid_parameters=" << d.has_valid_parameters
              << ", has_valid_covar=" << d.has_valid_covar
              << ", has_posdef_covar=" << d.has_posdef_covar
              << ", has_accurate_covar=" << d.has_accurate_covar
              << ", made_posdef=" << d.made_posdef
              << ", hesse_failed=" << d.hesse_failed
              << ", reached_call_limit=" << d.reached_call_limit
              << ", above_max_edm=" << d.above_max_edm
              << ", fmin=" << d.fmin
              << ", edm=" << d.edm
              << ", nfcn=" << d.nfcn
              << std::endl;

    if (!d.cov_eigs.empty()) {
        double min_eig = std::numeric_limits<double>::infinity();
        double max_eig = -std::numeric_limits<double>::infinity();
        std::size_t non_pos = 0;

        for (double eig : d.cov_eigs) {
            min_eig = std::min(min_eig, eig);
            max_eig = std::max(max_eig, eig);
            if (!(eig > 0.0)) ++non_pos;
        }

        std::cout << "[FIT] Minuit covariance eigenvalues: min=" << min_eig
                  << ", max=" << max_eig
                  << ", non_positive=" << non_pos
                  << ", cond=" << d.cond_number
                  << std::endl;
    }
}

void log_matrix_diagnostics(const std::string& label, const RealMatrix& M) {
    std::cout << "[FIT] Matrix " << label
              << " shape=" << M.rows() << "x" << M.cols()
              << std::endl;

    if (M.rows() == 0 || M.cols() == 0) {
        std::cout << "[FIT] Matrix " << label << " is empty" << std::endl;
        return;
    }

    if (M.rows() != M.cols()) {
        std::cout << "[FIT] Matrix " << label << " is not square" << std::endl;
        return;
    }

    try {
        RealMatrix sym = 0.5 * (M + M.transpose());
        EigenSystem eig = sym.eig();

        double min_eig = std::numeric_limits<double>::infinity();
        double max_eig = -std::numeric_limits<double>::infinity();
        std::size_t non_pos = 0;
        std::size_t tiny = 0;

        for (std::size_t i = 0; i < eig.D.rows(); ++i) {
            const double lambda = eig.D.at(i, i);
            min_eig = std::min(min_eig, lambda);
            max_eig = std::max(max_eig, lambda);
            if (!(lambda > 0.0)) ++non_pos;
            if (std::abs(lambda) < 1e-12) ++tiny;
        }

        std::cout << "[FIT] Matrix " << label
                  << " eig_min=" << min_eig
                  << ", eig_max=" << max_eig
                  << ", non_positive=" << non_pos
                  << ", tiny(|eig|<1e-12)=" << tiny
                  << std::endl;
    } catch (const std::exception& e) {
        std::cout << "[FIT] Matrix " << label
                  << " diagnostics failed: " << e.what()
                  << std::endl;
    }
}

RealMatrix invert_or_throw(const RealMatrix& M, const std::string& label) {
    log_matrix_diagnostics(label, M);
    try {
        return M.inv();
    } catch (const std::exception& e) {
        std::ostringstream oss;
        oss << "Failed to invert " << label << ": " << e.what();
        throw std::runtime_error(oss.str());
    }
}

} // namespace

namespace {

std::vector<std::size_t> make_fixed_p_indices(std::size_t p_dim) {
    std::vector<std::size_t> idx(p_dim);
    std::iota(idx.begin(), idx.end(), 0);
    return idx;
}

void fill_nan_profile_errors(FitResult& fr, std::size_t p_dim) {
    const double qnan = std::numeric_limits<double>::quiet_NaN();
    fr.p_hat_std.assign(p_dim, qnan);
    fr.p_hat_correlations = RealMatrix(p_dim, p_dim);
    for (std::size_t i = 0; i < p_dim; ++i) {
        for (std::size_t j = 0; j < p_dim; ++j) {
            fr.p_hat_correlations.at(i, j) = qnan;
        }
    }
}

void fill_profile_result_from_cov(const RealMatrix& cov_prof, FitResult& fr, std::size_t p_dim) {
    fr.p_hat_std.resize(p_dim);
    fr.p_hat_correlations = RealMatrix(p_dim, p_dim);

    for (std::size_t i = 0; i < p_dim; ++i) {
        fr.p_hat_std[i] = std::sqrt(std::max(0.0, cov_prof.at(i, i)));
    }

    for (std::size_t i = 0; i < p_dim; ++i) {
        for (std::size_t j = 0; j < p_dim; ++j) {
            const double si = fr.p_hat_std[i];
            const double sj = fr.p_hat_std[j];
            fr.p_hat_correlations.at(i, j) =
                (si > 0.0 && sj > 0.0) ? cov_prof.at(i, j) / (si * sj) : 0.0;
        }
    }
}

double choose_profile_fd_step(const fit_app::ParameterDefinition& def,
                              double x,
                              double step_scale)
{
    double h = step_scale * fit_app::safe_step(x, def.step_hint);

    if (def.limits.has_value()) {
        const auto [lo, hi] = def.limits.value();
        const double dist_lo = x - lo;
        const double dist_hi = hi - x;
        const double max_sym = 0.45 * std::max(0.0, std::min(dist_lo, dist_hi));
        if (max_sym > 0.0) {
            h = std::min(h, max_sym);
        }
    }

    if (!std::isfinite(h) || h <= 0.0) {
        h = 0.0;
    }
    return h;
}

RealMatrix regularize_spd(const RealMatrix& H, double rel_floor) {
    RealMatrix sym = 0.5 * (H + H.transpose());
    EigenSystem eig = sym.eig();

    double max_abs_eig = 0.0;
    for (std::size_t i = 0; i < eig.D.rows(); ++i) {
        max_abs_eig = std::max(max_abs_eig, std::abs(eig.D.at(i, i)));
    }

    const double floor = std::max(1e-10, rel_floor * std::max(1.0, max_abs_eig));

    RealMatrix Dreg(eig.D.rows(), eig.D.cols());
    for (std::size_t i = 0; i < eig.D.rows(); ++i) {
        Dreg.at(i, i) = std::max(eig.D.at(i, i), floor);
    }

    return eig.P * Dreg * eig.P.transpose();
}

double profiled_nll_at(const fit_app::IFitBackend& minimizer,
                       const fit_app::IObjectiveFunction& objective,
                       const std::vector<fit_app::ParameterDefinition>& defs,
                       const std::vector<double>& theta_anchor,
                       std::size_t p_dim,
                       const std::vector<double>& p_fixed,
                       const fit_app::FitOptions& profile_opt)
{
    std::vector<fit_app::ParameterDefinition> local_defs = defs;
    for (std::size_t i = 0; i < local_defs.size(); ++i) {
        local_defs[i].value = theta_anchor[i];
    }
    for (std::size_t i = 0; i < p_dim; ++i) {
        local_defs[i].value = p_fixed[i];
    }

    const std::vector<std::size_t> fixed_idx = make_fixed_p_indices(p_dim);

    fit_app::BackendFitResult prof =
        minimizer.minimize_with_fixed(objective, local_defs, profile_opt, fixed_idx, p_fixed);

    if (!prof.diagnostics.has_valid_parameters) {
        throw std::runtime_error("Profile minimization failed while building fallback Hessian.");
    }

    return prof.diagnostics.fmin;
}

RealMatrix numerical_profile_hessian(const fit_app::IFitBackend& minimizer,
                                     const fit_app::IObjectiveFunction& objective,
                                     const std::vector<fit_app::ParameterDefinition>& defs,
                                     const std::vector<double>& theta_hat,
                                     std::size_t p_dim,
                                     const fit_app::FitOptions& profile_opt,
                                     double step_scale,
                                     double f0)
{
    RealMatrix H(p_dim, p_dim);

    std::vector<double> p_hat(theta_hat.begin(), theta_hat.begin() + p_dim);
    std::vector<double> h(p_dim, 0.0);

    for (std::size_t i = 0; i < p_dim; ++i) {
        h[i] = choose_profile_fd_step(defs[i], p_hat[i], step_scale);
        if (!(h[i] > 0.0)) {
            std::ostringstream oss;
            oss << "Cannot build profile Hessian: finite-difference step collapsed for parameter "
                << defs[i].name;
            throw std::runtime_error(oss.str());
        }
    }

    for (std::size_t i = 0; i < p_dim; ++i) {
        auto p_plus = p_hat;
        auto p_minus = p_hat;
        p_plus[i] += h[i];
        p_minus[i] -= h[i];

        const double f_plus =
            profiled_nll_at(minimizer, objective, defs, theta_hat, p_dim, p_plus, profile_opt);
        const double f_minus =
            profiled_nll_at(minimizer, objective, defs, theta_hat, p_dim, p_minus, profile_opt);

        H.at(i, i) = (f_plus - 2.0 * f0 + f_minus) / (h[i] * h[i]);

        for (std::size_t j = i + 1; j < p_dim; ++j) {
            auto p_pp = p_hat;
            auto p_pm = p_hat;
            auto p_mp = p_hat;
            auto p_mm = p_hat;

            p_pp[i] += h[i]; p_pp[j] += h[j];
            p_pm[i] += h[i]; p_pm[j] -= h[j];
            p_mp[i] -= h[i]; p_mp[j] += h[j];
            p_mm[i] -= h[i]; p_mm[j] -= h[j];

            const double f_pp =
                profiled_nll_at(minimizer, objective, defs, theta_hat, p_dim, p_pp, profile_opt);
            const double f_pm =
                profiled_nll_at(minimizer, objective, defs, theta_hat, p_dim, p_pm, profile_opt);
            const double f_mp =
                profiled_nll_at(minimizer, objective, defs, theta_hat, p_dim, p_mp, profile_opt);
            const double f_mm =
                profiled_nll_at(minimizer, objective, defs, theta_hat, p_dim, p_mm, profile_opt);

            const double hij = (f_pp - f_pm - f_mp + f_mm) / (4.0 * h[i] * h[j]);
            H.at(i, j) = hij;
            H.at(j, i) = hij;
        }
    }

    return 0.5 * (H + H.transpose());
}

} // namespace

MLFitter::MLFitter(std::shared_ptr<LikelihoodContext> ctx, const ModelFn& model, MLFitOptions options)
    : fit_options_(options)
{
    const std::size_t p_dim = ctx->fp_defs.size();
    this->like_ = std::make_shared<BaseLikelihood>(model, ctx, p_dim);
}

FitResult MLFitter::maximum_likelihood_fit(const std::vector<double>& p0) {
    const auto defs = like_->get_param_defs();
    const std::size_t p_dim = p0.size();
    const std::size_t dim = defs.size();

    std::vector<fit_app::ParameterDefinition> theta0 = defs;
    for (std::size_t i = 0; i < p_dim; ++i) {
        theta0[i].value = p0[i];
    }

    auto f = fit_app::LambdaObjectiveFunction(
        [this](const std::vector<double>& theta) {
            return like_->nll(theta);
        },
        0.5
    );

    fit_app::FitOptions opt;
    opt.run_hesse = fit_options_.run_hesse;
    opt.verbose = fit_options_.verbose;
    if (fit_options_.strategy > 0) {
        opt.strategy = fit_options_.strategy;
    }
    if (fit_options_.max_fcn > 0) {
        opt.max_fcn = fit_options_.max_fcn;
    }
    if (fit_options_.tolerance > 0.0) {
        opt.tolerance = fit_options_.tolerance;
    }

#ifdef FIT_APP_HAS_RUN_MINOS
    opt.run_minos = fit_options_.request_minos;
#else
    if (fit_options_.request_minos) {
        std::cout << "[FIT] MINOS was requested, but fit_app::FitOptions has no run_minos flag in this build.\n"
                  << "[FIT] Define FIT_APP_HAS_RUN_MINOS only if your backend exposes opt.run_minos.\n";
    }
#endif

    std::unique_ptr<fit_app::IFitBackend> minimizer = fit_app::make_minuit_backend();
    fit_app::BackendFitResult res = minimizer->minimize(f, theta0, opt);

    log_fit_diagnostics(res.diagnostics);
    log_matrix_diagnostics("Minuit covariance", res.covariance);

    FitResult fr;
    fr.ell_hat = res.diagnostics.fmin;
    fr.p_hat.assign(res.values.begin(), res.values.begin() + p_dim);
    fr.eta_hat.assign(res.values.begin() + p_dim, res.values.end());

    bool have_profile_covariance = false;

    // 1) Chemin nominal : covariance Minuit valide
    if (res.diagnostics.has_valid_covar) {
        try {
            RealMatrix H = invert_or_throw(res.covariance, "Minuit covariance");

            RealMatrix H_p_p(p_dim, p_dim);
            RealMatrix H_p_eta(p_dim, dim - p_dim);
            RealMatrix H_eta_eta(dim - p_dim, dim - p_dim);

            for (std::size_t i = 0; i < p_dim; ++i) {
                for (std::size_t j = 0; j < p_dim; ++j) {
                    H_p_p.at(i, j) = H.at(i, j);
                }
            }

            for (std::size_t i = 0; i < p_dim; ++i) {
                for (std::size_t j = p_dim; j < dim; ++j) {
                    H_p_eta.at(i, j - p_dim) = H.at(i, j);
                }
            }

            for (std::size_t i = p_dim; i < dim; ++i) {
                for (std::size_t j = p_dim; j < dim; ++j) {
                    H_eta_eta.at(i - p_dim, j - p_dim) = H.at(i, j);
                }
            }

            RealMatrix cov_prof;
            if (dim == p_dim) {
                cov_prof = invert_or_throw(H_p_p, "H_p_p (no nuisance block)");
            } else {
                RealMatrix H_eta_eta_inv = invert_or_throw(H_eta_eta, "H_eta_eta");
                RealMatrix H_prof = H_p_p - H_p_eta * H_eta_eta_inv * H_p_eta.transpose();
                log_matrix_diagnostics("H_prof", H_prof);
                cov_prof = invert_or_throw(H_prof, "H_prof");
            }

            fill_profile_result_from_cov(cov_prof, fr, p_dim);
            have_profile_covariance = true;
        } catch (const std::exception& e) {
            std::cout << "[FIT] Failed to build profiled covariance from Minuit covariance: "
                      << e.what() << std::endl;
        }
    }

    // 2) Fallback propre : Hessienne numérique de la profile likelihood sur p
    if (!have_profile_covariance && fit_options_.allow_profile_hessian_fallback) {
        try {
            std::cout << "[FIT] Falling back to numerical profile Hessian on fit parameters.\n";

            fit_app::FitOptions profile_opt = opt;
            profile_opt.run_hesse = false;
            profile_opt.verbose = false;

            RealMatrix H_prof_num = numerical_profile_hessian(
                *minimizer,
                f,
                defs,
                res.values,
                p_dim,
                profile_opt,
                fit_options_.profile_hessian_step_scale,
                res.diagnostics.fmin
            );

            log_matrix_diagnostics("Numerical profile Hessian", H_prof_num);

            RealMatrix H_prof_reg = regularize_spd(
                H_prof_num,
                fit_options_.profile_hessian_eig_floor_rel
            );

            log_matrix_diagnostics("Regularized numerical profile Hessian", H_prof_reg);

            RealMatrix cov_prof = invert_or_throw(
                H_prof_reg,
                "Regularized numerical profile Hessian"
            );

            fill_profile_result_from_cov(cov_prof, fr, p_dim);
            have_profile_covariance = true;
        } catch (const std::exception& e) {
            std::cout << "[FIT] Numerical profile-Hessian fallback failed: "
                      << e.what() << std::endl;
        }
    }

    // 3) Si tout échoue, on garde le MLE mais on marque les erreurs comme indisponibles
    if (!have_profile_covariance) {
        std::cout << "[FIT] No covariance could be constructed. Returning MLE with NaN errors.\n";
        fill_nan_profile_errors(fr, p_dim);
    }

    // Pour la suite du workflow, un minimum avec paramètres valides est déjà exploitable
    this->master_fit_success = res.diagnostics.has_valid_parameters;
    this->master_fit_result = fr;
    return fr;
}

// FitResult MLFitter::maximum_likelihood_fit(const std::vector<double>& p0) {
//     const auto defs = like_->get_param_defs();
//     const std::size_t p_dim = p0.size();
//     const std::size_t dim = defs.size();

//     std::vector<fit_app::ParameterDefinition> theta0 = defs;
//     for (std::size_t i = 0; i < p_dim; ++i) {
//         theta0[i].value = p0[i];
//     }

//     auto f = fit_app::LambdaObjectiveFunction(
//         [this](const std::vector<double>& theta) {
//             return like_->nll(theta);
//         },
//         0.5
//     );

//     fit_app::FitOptions opt;
//     opt.run_hesse = fit_options_.run_hesse;
//     opt.verbose = fit_options_.verbose;
//     if (fit_options_.strategy > 0) {
//         opt.strategy = fit_options_.strategy;
//     }
//     if (fit_options_.max_fcn > 0) {
//         opt.max_fcn = fit_options_.max_fcn;
//     }
//     if (fit_options_.tolerance > 0.0) {
//         opt.tolerance = fit_options_.tolerance;
//     }

// #ifdef FIT_APP_HAS_RUN_MINOS
//     opt.run_minos = fit_options_.request_minos;
// #else
//     if (fit_options_.request_minos) {
//         std::cout << "[FIT] MINOS was requested, but fit_app::FitOptions has no run_minos flag in this build.\n"
//                   << "[FIT] Define FIT_APP_HAS_RUN_MINOS only if your backend exposes opt.run_minos.\n";
//     }
// #endif

//     if (fit_options_.request_minos) {
//         std::cout << "[FIT] Requesting MINOS for the master Minuit fit.\n";
//     }
//     if (!fit_options_.run_hesse) {
//         std::cout << "[FIT] Warning: maximum_likelihood_fit() still uses the covariance matrix to build profiled errors.\n"
//                   << "[FIT] Disabling HESSE may therefore make the fit fail if the backend does not return a valid covariance anyway.\n";
//     }

//     std::unique_ptr<fit_app::IFitBackend> minimizer = fit_app::make_minuit_backend();
//     fit_app::BackendFitResult res = minimizer->minimize(f, theta0, opt);

//     log_fit_diagnostics(res.diagnostics);
//     log_matrix_diagnostics("Minuit covariance", res.covariance);

//     if (!res.diagnostics.has_valid_covar) {
//         throw std::runtime_error(
//             "Minuit returned no valid covariance matrix. "
//             "This usually means the likelihood has flat or redundant directions."
//         );
//     }

//     FitResult fr;
//     fr.ell_hat = res.diagnostics.fmin;
//     fr.p_hat.assign(res.values.begin(), res.values.begin() + p_dim);
//     fr.eta_hat.assign(res.values.begin() + p_dim, res.values.end());

//     // bloc p×p de la covariance complète
//     RealMatrix H = invert_or_throw(res.covariance, "Minuit covariance");
//     RealMatrix H_p_p(p_dim, p_dim);
//     RealMatrix H_p_eta(p_dim, dim - p_dim);
//     RealMatrix H_eta_eta(dim - p_dim, dim - p_dim);

//     for (std::size_t i = 0; i < p_dim; ++i) {
//         for (std::size_t j = 0; j < p_dim; ++j) {
//             H_p_p.at(i, j) = H.at(i, j);
//         }
//     }

//     for (std::size_t i = 0; i < p_dim; ++i) {
//         for (std::size_t j = p_dim; j < dim; ++j) {
//             H_p_eta.at(i, j - p_dim) = H.at(i, j);
//         }
//     }

//     for (std::size_t i = p_dim; i < dim; ++i) {
//         for (std::size_t j = p_dim; j < dim; ++j) {
//             H_eta_eta.at(i - p_dim, j - p_dim) = H.at(i, j);
//         }
//     }

//     RealMatrix cov_prof;
//     if (dim == p_dim) {
//         cov_prof = invert_or_throw(H_p_p, "H_p_p (no nuisance block)");
//     } else {
//         RealMatrix H_eta_eta_inv = invert_or_throw(H_eta_eta, "H_eta_eta");
//         RealMatrix H_prof = H_p_p - H_p_eta * H_eta_eta_inv * H_p_eta.transpose();
//         log_matrix_diagnostics("H_prof", H_prof);
//         cov_prof = invert_or_throw(H_prof, "H_prof");
//     }

//     fr.p_hat_std.resize(p_dim);
//     fr.p_hat_correlations = RealMatrix(p_dim, p_dim);

//     for (std::size_t i = 0; i < p_dim; ++i) {
//         fr.p_hat_std[i] = std::sqrt(std::max(0.0, cov_prof.at(i, i)));
//     }
//     for (std::size_t i = 0; i < p_dim; ++i) {
//         for (std::size_t j = 0; j < p_dim; ++j) {
//             const double si = fr.p_hat_std[i];
//             const double sj = fr.p_hat_std[j];
//             fr.p_hat_correlations.at(i, j) =
//                 (si > 0.0 && sj > 0.0) ? cov_prof.at(i, j) / (si * sj) : 0.0;
//         }
//     }

//     this->master_fit_success = res.diagnostics.ok;
//     this->master_fit_result = fr;
//     return fr;
// }

Contour MLFitter::contour(std::size_t x_id, std::size_t y_id, double z, std::array<double, 4> bounds, ContourOptions options) const {
    if (!this->master_fit_success)
        LOG_ERROR("InvalidState", "ML fit must have converged before contour computation is available.");

    ContourConfig cc;
    cc.fr = this->master_fit_result;
    cc.x_id = x_id;
    cc.y_id = y_id;
    cc.primary_contour_method = options.primary_contour_method;
    cc.fallback_contour_method = options.fallback_contour_method;
    cc.profiling_method = options.profiling_method;

    ContourEngine ce(this->like_, cc);
    Contour cl = ce.compute_contour(z, bounds, options.resolution);

    return cl;
}