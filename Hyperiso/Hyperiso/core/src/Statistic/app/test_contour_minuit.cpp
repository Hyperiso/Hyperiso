#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "StatisticManager.h"
#include "ObservableInterfaceAdapter2.h"
#include "StatCorrelationProxy.h"
#include "StatParameterProxy.h"
#include "ObservableInterface.h"
#include "StatParamSourcesProxy.h"
#include "Fit.h"
#include "Likelihood.h"

#include "minuit-cpp/FCNBase.hh"
#include "minuit-cpp/FunctionMinimum.hh"
#include "minuit-cpp/MnContours.hh"
#include "minuit-cpp/MnEigen.hh"
#include "minuit-cpp/MnHesse.hh"
#include "minuit-cpp/MnMigrad.hh"
#include "minuit-cpp/MnMinos.hh"
#include "minuit-cpp/MnStrategy.hh"
#include "minuit-cpp/MnUserCovariance.hh"
#include "minuit-cpp/MnUserParameters.hh"
#include "minuit-cpp/MnUserParameterState.hh"
#include "minuit-cpp/MinosError.hh"

namespace M2 = MinuitCpp;

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------

template <class T>
static std::string stream_str(const T& x) {
    std::ostringstream oss;
    oss << x;
    return oss.str();
}

static void print_vec(const std::vector<double>& vec) {
    std::cout << "[ ";
    for (std::size_t i = 0; i < vec.size(); ++i) {
        std::cout << std::setprecision(17) << vec[i]
                  << (i + 1 == vec.size() ? " " : ", ");
    }
    std::cout << "]\n";
}

static std::vector<double> linspace(double a, double b, std::size_t n) {
    std::vector<double> out(n);
    if (n == 0) return out;
    if (n == 1) {
        out[0] = a;
        return out;
    }
    for (std::size_t i = 0; i < n; ++i) {
        out[i] = a + (b - a) * double(i) / double(n - 1);
    }
    return out;
}

static double safe_step(double value, double scale_hint) {
    double a = std::fabs(value);
    double s = std::fabs(scale_hint);

    double step = 0.0;
    if (std::isfinite(s) && s > 0.0) step = 0.05 * s;
    if (std::isfinite(a) && a > 0.0) step = std::max(step, 0.01 * a);
    if (!std::isfinite(step) || step <= 0.0) step = 1e-3;

    return step;
}

static void save_bestfit_csv(const std::string& path,
                             const std::vector<std::string>& names,
                             const Vector& vals,
                             const Vector& errs) {
    std::ofstream out(path);
    out << "name,value,error\n";
    for (std::size_t i = 0; i < vals.size(); ++i) {
        out << names[i] << ","
            << std::setprecision(17) << vals[i] << ","
            << errs[i] << "\n";
    }
}

static void save_contours_csv(const std::string& path,
                              const std::string& xname,
                              const std::string& yname,
                              const std::vector<std::pair<double, double>>& c68,
                              const std::vector<std::pair<double, double>>& c95) {
    std::ofstream out(path);
    out << "# x=" << xname << "\n";
    out << "# y=" << yname << "\n";
    out << "cl,x,y\n";
    for (const auto& p : c68) {
        out << "0.683," << std::setprecision(17) << p.first << "," << p.second << "\n";
    }
    for (const auto& p : c95) {
        out << "0.95," << std::setprecision(17) << p.first << "," << p.second << "\n";
    }
}

static void save_grid_csv(const std::string& path,
                          const std::string& xname,
                          const std::string& yname,
                          const std::vector<double>& xs,
                          const std::vector<double>& ys,
                          const std::vector<double>& z) {
    const std::size_t nx = xs.size();
    const std::size_t ny = ys.size();

    std::ofstream out(path);
    out << "# x=" << xname << "\n";
    out << "# y=" << yname << "\n";
    out << "x,y,delta_nll\n";
    out << std::setprecision(17);

    for (std::size_t iy = 0; iy < ny; ++iy) {
        for (std::size_t ix = 0; ix < nx; ++ix) {
            out << xs[ix] << "," << ys[iy] << "," << z[iy * nx + ix] << "\n";
        }
    }
}

// -----------------------------------------------------------------------------
// Minuit wrappers
// -----------------------------------------------------------------------------

struct ParamLimit {
    std::size_t idx;
    double low;
    double high;
};

struct MinuitFitOptions {
    double up = 0.5;               // NLL => 0.5 for 1D 1σ
    unsigned strategy = 2;         // robust
    unsigned max_fcn = 100000;
    double tolerance = 0.2;        // Minuit "toler", not raw EDM
    bool run_hesse = true;
    unsigned hesse_maxcalls = 0;
    bool run_minos_checks = false;
    unsigned minos_maxcalls = 200000;
    bool verbose = true;
};

struct MinuitJointFit {
    Vector x_hat;                       // full vector [p, eta]
    std::vector<double> x_err;          // HESSE errors
    std::vector<double> cov_eigs;       // covariance eigenvalues
    double cond_number = std::numeric_limits<double>::infinity();

    double fmin = std::numeric_limits<double>::quiet_NaN();
    double edm  = std::numeric_limits<double>::quiet_NaN();
    int nfcn    = -1;

    bool ok = false;
    bool has_valid_covar = false;
    bool has_posdef_covar = false;
    bool has_accurate_covar = false;
    bool made_posdef = false;

    RealMatrix cov;
    std::shared_ptr<M2::FunctionMinimum> min;
};

class GenericFCN final : public M2::FCNBase {
public:
    GenericFCN(std::function<double(const std::vector<double>&)> f, double up)
        : f_(std::move(f)), up_(up) {}

    double operator()(const std::vector<double>& x) const override {
        try {
            const double v = f_(x);
            return std::isfinite(v) ? v : 1e300;
        } catch (...) {
            return 1e300;
        }
    }

    double Up() const override { return up_; }

private:
    std::function<double(const std::vector<double>&)> f_;
    double up_;
};

static void log_minuit_summary(const std::string& tag, const M2::FunctionMinimum& min) {
    std::cout << "\n=== [" << tag << "] FunctionMinimum ===\n";
    std::cout << "IsValid            = " << min.IsValid() << "\n";
    std::cout << "HasValidParameters = " << min.HasValidParameters() << "\n";
    std::cout << "HasValidCovariance = " << min.HasValidCovariance() << "\n";
    std::cout << "HasAccurateCovar   = " << min.HasAccurateCovar() << "\n";
    std::cout << "HasPosDefCovar     = " << min.HasPosDefCovar() << "\n";
    std::cout << "HasMadePosDefCovar = " << min.HasMadePosDefCovar() << "\n";
    std::cout << "HesseFailed        = " << min.HesseFailed() << "\n";
    std::cout << "Fval               = " << std::setprecision(17) << min.Fval() << "\n";
    std::cout << "EDM                = " << std::setprecision(17) << min.Edm() << "\n";
    std::cout << "NFcn               = " << min.NFcn() << "\n";
    std::cout << "Up                 = " << min.Up() << "\n";
}

static MinuitJointFit minuit_migrad_hesse(
    const std::function<double(const std::vector<double>&)>& f,
    const std::vector<std::string>& names,
    const std::vector<double>& x0,
    const std::vector<double>& scale_hints,
    const std::vector<ParamLimit>& limits,
    const MinuitFitOptions& opt,
    const std::vector<unsigned>& minos_check_indices = {}
) {
    if (x0.size() != names.size() || x0.size() != scale_hints.size()) {
        throw std::invalid_argument("minuit_migrad_hesse: names/x0/scale_hints size mismatch");
    }

    GenericFCN fcn(f, opt.up);

    M2::MnUserParameters upar;
    for (std::size_t i = 0; i < x0.size(); ++i) {
        const double step = safe_step(x0[i], scale_hints[i]);
        upar.Add(names[i].c_str(), x0[i], step);
    }

    for (const auto& lim : limits) {
        if (lim.idx < names.size()) {
            upar.SetLimits(names[lim.idx].c_str(), lim.low, lim.high);
        }
    }

    M2::MnMigrad migrad(fcn, upar, opt.strategy);
    M2::FunctionMinimum min = migrad(opt.max_fcn, opt.tolerance);

    if (opt.run_hesse) {
        M2::MnHesse hesse(opt.strategy);
        hesse(fcn, min, opt.hesse_maxcalls);
    }

    if (opt.verbose) {
        log_minuit_summary("MIGRAD+HESSE", min);
    }

    std::vector<double> eigs;
    double cond = std::numeric_limits<double>::infinity();

    if (min.HasValidCovariance()) {
        M2::MnEigen eigen;
        eigs = eigen(min.UserState().Covariance());

        if (!eigs.empty() && eigs.front() > 0.0) {
            cond = eigs.back() / eigs.front();
        }

        if (opt.verbose) {
            if (!eigs.empty()) {
                std::cout << "Cov eigen min/max = "
                          << std::setprecision(6) << eigs.front()
                          << " / " << eigs.back()
                          << "   (cond ~ " << cond << ")\n";
            } else {
                std::cout << "Cov eigenvalues unavailable\n";
            }
        }
    }

    if (opt.run_minos_checks && min.IsValid() && !minos_check_indices.empty()) {
        GenericFCN fcn_minos(f, opt.up);
        M2::MnMinos minos(fcn_minos, min, opt.strategy);

        for (unsigned idx : minos_check_indices) {
            if (idx >= names.size()) continue;
            auto me = minos.Minos(idx, opt.minos_maxcalls);
            auto e = me();
            std::cout << "MINOS check [" << names[idx] << "] valid=" << me.IsValid()
                      << " err=(" << e.first << ", " << e.second << ")\n";
        }
    }

    MinuitJointFit out;
    out.fmin = min.Fval();
    out.edm  = min.Edm();
    out.nfcn = min.NFcn();
    out.ok   = min.IsValid();

    out.has_valid_covar    = min.HasValidCovariance();
    out.has_posdef_covar   = min.HasPosDefCovar();
    out.has_accurate_covar = min.HasAccurateCovar();
    out.made_posdef        = min.HasMadePosDefCovar();

    out.cov_eigs = std::move(eigs);
    out.cond_number = cond;
    out.min = std::make_shared<M2::FunctionMinimum>(min);

    const auto& st = min.UserState();
    const std::size_t n = x0.size();

    out.x_hat.resize(n);
    out.x_err.assign(n, 0.0);
    out.cov = RealMatrix(n, n);

    for (std::size_t i = 0; i < n; ++i) {
        out.x_hat[i] = st.Value(names[i].c_str());
        out.x_err[i] = st.Error(names[i].c_str());
    }

    if (min.HasValidCovariance()) {
        const auto& cov = st.Covariance();
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                out.cov.at(i, j) = cov(i, j);
            }
        }
    }

    return out;
}

static bool try_mncontours_with_refit(
    const std::function<double(const std::vector<double>&)>& f_joint,
    const std::vector<std::string>& names,
    const MinuitJointFit& mj_base,
    const std::vector<double>& scale_hints,
    const std::vector<ParamLimit>& limits,
    unsigned px,
    unsigned py,
    double up_contour,
    unsigned npoints,
    std::vector<std::pair<double, double>>& out_points
) {
    if (!mj_base.ok || !mj_base.min) return false;

    // Important: refit court avec Up cohérent pour éviter le mismatch fcn.Up() / min.Up()
    MinuitFitOptions opt2;
    opt2.up = up_contour;
    opt2.strategy = 2;
    opt2.max_fcn = 30000;
    opt2.tolerance = 0.2;
    opt2.run_hesse = true;
    opt2.verbose = false;

    std::vector<double> refit_scales = scale_hints;
    const double scale_factor = std::sqrt(up_contour / 0.5);

    for (std::size_t i = 0; i < refit_scales.size() && i < mj_base.x_err.size(); ++i) {
        const double err = std::fabs(mj_base.x_err[i]);
        if (std::isfinite(err) && err > 0.0) {
            refit_scales[i] = std::max(refit_scales[i], err * scale_factor);
        }
    }

    MinuitJointFit mj_up = minuit_migrad_hesse(
        f_joint, names, mj_base.x_hat, refit_scales, limits, opt2
    );

    if (!mj_up.ok || !mj_up.min) return false;

    GenericFCN fcn(f_joint, up_contour);
    M2::MnContours contours(fcn, *mj_up.min, 2);
    out_points = contours(px, py, npoints);

    return out_points.size() >= 4;
}

static double profiled_nll_at_fixed_xy(
    const std::function<double(const std::vector<double>&)>& f_joint,
    const std::vector<std::string>& names,
    const std::vector<double>& x_start,
    const std::vector<double>& scale_hints,
    const std::vector<ParamLimit>& limits,
    unsigned px,
    unsigned py,
    double xval,
    double yval,
    unsigned strategy,
    unsigned max_fcn
) {
    GenericFCN fcn(f_joint, 0.5);

    M2::MnUserParameters upar;
    for (std::size_t i = 0; i < x_start.size(); ++i) {
        upar.Add(names[i].c_str(), x_start[i], safe_step(x_start[i], scale_hints[i]));
    }

    for (const auto& lim : limits) {
        if (lim.idx < names.size()) {
            upar.SetLimits(names[lim.idx].c_str(), lim.low, lim.high);
        }
    }

    M2::MnMigrad migrad(fcn, upar, strategy);
    migrad.SetValue(px, xval);
    migrad.SetValue(py, yval);
    migrad.Fix(px);
    migrad.Fix(py);

    M2::FunctionMinimum min = migrad(max_fcn, 0.5);

    if (!min.IsValid()) {
        std::vector<double> x = x_start;
        x[px] = xval;
        x[py] = yval;
        return f_joint(x);
    }

    return min.Fval();
}

// -----------------------------------------------------------------------------
// Local estimator using Minuit
// -----------------------------------------------------------------------------

struct JointFitOutput {
    FitResult fr;
    MinuitJointFit mj;
    std::vector<std::string> names;       // full order [p, eta]
    std::vector<double> scale_hints;      // full order [p, eta]
    std::vector<ParamLimit> limits;       // full order indices
};

class MinuitMLEstimatorLocal {
public:
    using ModelFn = ProfiledLikelihood::ModelFn;

    MinuitMLEstimatorLocal(LikelihoodContext ctx,
                           ModelFn model,
                           std::size_t max_fcn,
                           double tolerance,
                           unsigned strategy)
        : like_(std::move(ctx))
        , model_(std::move(model))
        , max_fcn_(max_fcn)
        , tolerance_(tolerance)
        , strategy_(strategy) {}

    std::function<double(const std::vector<double>&)> make_joint_f(std::size_t p_dim) const {
        return [this, p_dim](const std::vector<double>& x) -> double {
            Vector p(x.begin(), x.begin() + p_dim);
            Vector eta(x.begin() + p_dim, x.end());
            return nll(p, eta);
        };
    }

    JointFitOutput fit_joint_with_minuit(
        const std::vector<ParamId>& p_ids,
        const std::vector<ParamId>& eta_ids,
        const Vector& p0
    ) const {
        const std::size_t p_dim = p0.size();

        Vector eta0 = like_.nuisance_central_values;
        Vector eta_scales = like_.nuisance_dist->get_stds();

        if (eta0.size() != eta_scales.size()) {
            throw std::runtime_error("eta central values and eta stds do not have same size");
        }

        std::vector<double> x0;
        x0.reserve(p_dim + eta0.size());
        x0.insert(x0.end(), p0.begin(), p0.end());
        x0.insert(x0.end(), eta0.begin(), eta0.end());

        std::vector<double> scale_hints;
        scale_hints.reserve(x0.size());

        for (std::size_t i = 0; i < p_dim; ++i) {
            double hint = std::fabs(p0[i]);
            if (hint < 1e-3) hint = 0.01;
            scale_hints.push_back(hint);
        }
        for (double s : eta_scales) {
            scale_hints.push_back(std::max(1e-12, std::fabs(s)));
        }

        std::vector<std::string> names;
        names.reserve(x0.size());
        for (const auto& pid : p_ids) names.push_back(stream_str(pid));
        for (const auto& pid : eta_ids) names.push_back(stream_str(pid));

        std::vector<ParamLimit> limits;
        for (std::size_t i = 0; i < p_dim; ++i) {
            if (names[i].find("FCONST") != std::string::npos) {
                limits.push_back(ParamLimit{i, 0.05, 0.35});
            }
        }
        for (std::size_t i = p_dim; i < names.size(); ++i) {
        const std::string& nm = names[i];
        const double c = x0[i];
        const double s = std::max(1e-12, std::fabs(scale_hints[i]));

        if (nm.find("SMINPUTS:3") != std::string::npos) {
            limits.push_back(ParamLimit{i, 0.05, 0.30});
        } else if (nm.find("MASS:") != std::string::npos ||
                nm.find("FLIFE:") != std::string::npos ||
                nm.find("FCONST:") != std::string::npos ||
                nm.find("FMASS:") != std::string::npos ||
                nm.find("SMINPUTS:5") != std::string::npos ||
                nm.find("SMINPUTS:6") != std::string::npos) {
            limits.push_back(ParamLimit{i, std::max(1e-12, c - 5.0 * s), c + 5.0 * s});
        }
    }
        auto f_joint = make_joint_f(p_dim);

        MinuitFitOptions opt;
        opt.up = 0.5;
        opt.strategy = strategy_;
        opt.max_fcn = static_cast<unsigned>(max_fcn_);
        opt.tolerance = tolerance_;
        opt.run_hesse = true;
        opt.run_minos_checks = (p_dim >= 2);
        opt.minos_maxcalls = std::max<unsigned>(200000, opt.max_fcn);
        opt.verbose = true;

        std::vector<unsigned> minos_check;
        if (p_dim >= 2) minos_check = {0u, 1u};

        MinuitJointFit mj = minuit_migrad_hesse(
            f_joint, names, x0, scale_hints, limits, opt, minos_check
        );

        FitResult fr;
        fr.ell_hat = mj.fmin;

        fr.p_hat.assign(mj.x_hat.begin(), mj.x_hat.begin() + p_dim);
        fr.eta_hat.assign(mj.x_hat.begin() + p_dim, mj.x_hat.end());

        fr.p_hat_std.assign(p_dim, 0.0);
        fr.p_hat_correlations = RealMatrix(p_dim, p_dim);

        if (mj.min && mj.min->HasValidCovariance()) {
            const auto& cov = mj.min->UserState().Covariance();

            for (std::size_t i = 0; i < p_dim; ++i) {
                fr.p_hat_std[i] = std::sqrt(std::max(0.0, cov(i, i)));
            }

            for (std::size_t i = 0; i < p_dim; ++i) {
                for (std::size_t j = 0; j < p_dim; ++j) {
                    const double di = std::sqrt(std::max(0.0, cov(i, i)));
                    const double dj = std::sqrt(std::max(0.0, cov(j, j)));
                    fr.p_hat_correlations.at(i, j) =
                        (di > 0.0 && dj > 0.0) ? (cov(i, j) / (di * dj)) : 0.0;
                }
            }
        } else {
            for (std::size_t i = 0; i < p_dim && i < mj.x_err.size(); ++i) {
                fr.p_hat_std[i] = mj.x_err[i];
                fr.p_hat_correlations.at(i, i) = 1.0;
            }
        }

        return JointFitOutput{fr, mj, names, scale_hints, limits};
    }

private:
    double nll(const Vector& p, const Vector& eta) const {
        Vector pred = model_(p, eta);

        Vector r(pred.size());
        for (std::size_t i = 0; i < pred.size(); ++i) {
            r[i] = pred[i] - like_.exp_obs_values[i];
        }

        const double ell_obs = like_.exp_obs_dist->logpdf(r);
        const double ell_eta = like_.nuisance_dist->logpdf(eta);
        return -(ell_obs + ell_eta);
    }

    LikelihoodContext like_;
    ModelFn model_;
    std::size_t max_fcn_;
    double tolerance_;
    unsigned strategy_;
};

// -----------------------------------------------------------------------------
// main
// -----------------------------------------------------------------------------

int main(int argc, char** argv) {
    HyperisoMaster hyp;
    HyperisoConfig config_hyp;
    config_hyp.model = Model::SM;
    hyp.init("lha/si_input.flha", config_hyp);

    auto oint = std::make_shared<ObservableInterface>();
    oint->add_observable(ObservableMapper::to_id(Observables::BR_BS_MUMU_UNTAG), QCDOrder::LO, true)
        .add_observable(ObservableMapper::to_id(Observables::BR_BD_MUMU), QCDOrder::LO, true);

    StatisticConfig config;
    config.MC_draws = 100;

    // IMPORTANT:
    // ici MLE_max_iter sert maintenant de max_fcn Minuit,
    // et MLE_tol sert de "tolerance" Minuit.
    config.MLE_max_iter = 150000;
    config.MLE_tol = 0.2;

    config.p_specs = {
        ParamId{ParameterType::FLAVOR, "FCONST", {511, 1}},
        ParamId{ParameterType::FLAVOR, "FCONST", {531, 1}}
    };

    auto model = std::make_shared<ObservableInterfaceAdapterObs>(oint);

    StatisticManager stat(
        config,
        model,
        std::make_shared<StatCorrelationProxy>(),
        std::make_shared<StatParameterProxy>(),
        std::make_shared<StatParamSourcesProxy>()
    );

    LOG_INFO("fill_cache #1");
    stat.fill_cache();

    auto start_u = std::chrono::steady_clock::now();
    stat.compute_uncertainties();
    auto stop_u = std::chrono::steady_clock::now();

    auto us_u = std::chrono::duration_cast<std::chrono::microseconds>(stop_u - start_u).count();
    std::cout << "Uncertainty estimation time: " << us_u << " us\n";

    // refresh cache after uncertainties
    LOG_INFO("fill_cache #2");
    stat.fill_cache();

    // Build context like compute_MLE()
    auto p_specs_map = stat.get_p_specs();
    auto eta_specs_real = stat.get_all_obss_deps();
    for (const auto& [pid, _] : p_specs_map) eta_specs_real.erase(pid);
    auto exp_obs_map = stat.get_obs_exp();

    auto unz_p   = unzip(p_specs_map);
    auto unz_eta = unzip(eta_specs_real);
    auto unz_obs = unzip(exp_obs_map);

    std::vector<ParamId> p_ids = unz_p.ids;
    std::vector<ParamId> eta_ids = unz_eta.ids;
    std::vector<BinnedObservableId> obs_ids = unz_obs.ids;

    auto nuisance_dist = stat.build_nuisance_distribution();
    auto exp_obs_dist  = stat.build_exp_data_distribution();

    if (nuisance_dist->get_stds().size() != unz_eta.vals.size()) {
        std::cerr << "[ERROR] nuisance std size = " << nuisance_dist->get_stds().size()
                  << " but eta central size = " << unz_eta.vals.size() << "\n";
        return 3;
    }

    if (exp_obs_dist->dim() != unz_obs.vals.size()) {
        std::cerr << "[ERROR] exp obs dim = " << exp_obs_dist->dim()
                  << " but values size = " << unz_obs.vals.size() << "\n";
        return 4;
    }

    LikelihoodContext ctx;
    ctx.nuisance_dist = std::move(nuisance_dist);
    ctx.exp_obs_dist  = std::move(exp_obs_dist);
    ctx.nuisance_central_values = unz_eta.vals;
    ctx.exp_obs_values = unz_obs.vals;

    auto model_fn = [model, obs_ids, p_ids, eta_ids](const Vec& p_vec, const Vec& eta_vec) -> Vec {
        auto pred_map = model->predict_optimized(zip(p_ids, p_vec), zip(eta_ids, eta_vec));

        Vec out;
        out.reserve(obs_ids.size());

        for (const auto& bid : obs_ids) {
            const auto& vec = pred_map.at(bid.s);

            auto it = std::find_if(vec.begin(), vec.end(), [&](const ObservableValue& ov) {
                auto bin = ov.bin.value_or(std::pair<double, double>{0., 0.});
                return bin == bid.p;
            });

            if (it == vec.end()) {
                throw std::runtime_error("Missing predicted observable/bin");
            }

            out.push_back(it->value);
        }

        return out;
    };

    auto start_m = std::chrono::steady_clock::now();
    MinuitMLEstimatorLocal est(std::move(ctx), model_fn, config.MLE_max_iter, config.MLE_tol, 2);
    JointFitOutput fit_out = est.fit_joint_with_minuit(p_ids, eta_ids, unz_p.vals);
    auto stop_m = std::chrono::steady_clock::now();

    const auto& fr = fit_out.fr;
    const auto& mj = fit_out.mj;

    auto us_m = std::chrono::duration_cast<std::chrono::microseconds>(stop_m - start_m).count();
    std::cout << "\nMLE (Minuit) fitting time: " << us_m << " us\n";

    std::cout << "ell_hat = " << std::setprecision(17) << fr.ell_hat << "\n";
    std::cout << "p_hat = "; print_vec(fr.p_hat);
    std::cout << "p_hat_std = "; print_vec(fr.p_hat_std);
    std::cout << "p_hat_correlations:\n" << fr.p_hat_correlations << "\n";

    if (!mj.ok) {
        std::cerr << "[ERROR] Minuit fit invalid.\n";
        return 5;
    }

    if (!mj.has_valid_covar || !mj.has_posdef_covar) {
        std::cerr << "[WARN] Covariance is not fully healthy."
                  << " valid=" << mj.has_valid_covar
                  << " posdef=" << mj.has_posdef_covar
                  << " accurate=" << mj.has_accurate_covar
                  << " cond=" << mj.cond_number << "\n";
    }

    // Save best-fit always
    {
        std::vector<std::string> p_names;
        for (const auto& pid : p_ids) p_names.push_back(stream_str(pid));
        save_bestfit_csv("bestfit.csv", p_names, fr.p_hat, fr.p_hat_std);
        std::cout << "[INFO] Wrote bestfit.csv\n";
    }

    // Contours only if exactly 2 fit parameters
    if (p_ids.size() == 2) {
        const double up68 = 2.30 / 2.0;   // 1.15
        const double up95 = 5.99 / 2.0;   // 2.995

        const unsigned px = 0;
        const unsigned py = 1;
        const unsigned npoints = 80;

        auto f_joint = est.make_joint_f(/*p_dim=*/2);

        std::vector<std::pair<double, double>> c68;
        std::vector<std::pair<double, double>> c95;

        bool ok68 = try_mncontours_with_refit(
            f_joint, fit_out.names, mj, fit_out.scale_hints, fit_out.limits,
            px, py, up68, npoints, c68
        );

        bool ok95 = try_mncontours_with_refit(
            f_joint, fit_out.names, mj, fit_out.scale_hints, fit_out.limits,
            px, py, up95, npoints, c95
        );

        const std::string xname = stream_str(p_ids[0]);
        const std::string yname = stream_str(p_ids[1]);

        if (ok68 && ok95) {
            save_contours_csv("contours.csv", xname, yname, c68, c95);
            std::cout << "[INFO] Wrote contours.csv\n";
        } else {
            std::cerr << "[WARN] MnContours failed; fallback to grid scan.\n";

            double x0 = fr.p_hat[0];
            double y0 = fr.p_hat[1];
            double sx = std::max(0.01, fr.p_hat_std[0]);
            double sy = std::max(0.01, fr.p_hat_std[1]);

            double xlo = std::max(0.05, x0 - 5.0 * sx);
            double xhi = std::min(0.35, x0 + 5.0 * sx);
            double ylo = std::max(0.05, y0 - 5.0 * sy);
            double yhi = std::min(0.35, y0 + 5.0 * sy);

            if (!(xhi > xlo)) { xlo = 0.10; xhi = 0.30; }
            if (!(yhi > ylo)) { ylo = 0.10; yhi = 0.30; }

            std::vector<double> xs = linspace(xlo, xhi, 81);
            std::vector<double> ys = linspace(ylo, yhi, 81);
            std::vector<double> z;
            z.reserve(xs.size() * ys.size());

            for (double y : ys) {
                for (double x : xs) {
                    double fprof = profiled_nll_at_fixed_xy(
                        f_joint,
                        fit_out.names,
                        mj.x_hat,
                        fit_out.scale_hints,
                        fit_out.limits,
                        px, py,
                        x, y,
                        /*strategy=*/2,
                        /*max_fcn=*/40000
                    );

                    z.push_back(fprof - fr.ell_hat);
                }
            }

            save_grid_csv("grid.csv", xname, yname, xs, ys, z);
            std::cout << "[INFO] Wrote grid.csv\n";
        }
    }

    return 0;
}