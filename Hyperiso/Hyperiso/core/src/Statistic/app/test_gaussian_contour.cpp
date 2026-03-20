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
#include <utility>
#include <vector>

#include "StatisticManager.h"
#include "ObservableInterfaceAdapter2.h"
#include "StatCorrelationProxy.h"
#include "StatParameterProxy.h"
#include "ObservableInterface.h"
#include "StatParamSourcesProxy.h"
#include "Fit.h"
#include "BaseLikelihood.h"

#include "minuit-cpp/FCNBase.hh"
#include "minuit-cpp/FunctionMinimum.hh"
#include "minuit-cpp/MnContours.hh"
#include "minuit-cpp/MnEigen.hh"
#include "minuit-cpp/MnHesse.hh"
#include "minuit-cpp/MnMigrad.hh"
#include "minuit-cpp/MnUserCovariance.hh"
#include "minuit-cpp/MnUserParameters.hh"
#include "minuit-cpp/MnUserParameterState.hh"

namespace M2 = MinuitCpp;

namespace contour_app {

// -----------------------------------------------------------------------------
// Small utilities
// -----------------------------------------------------------------------------

template <class T>
std::string to_string_any(const T& x) {
    std::ostringstream oss;
    oss << x;
    return oss.str();
}

void print_vec(const std::vector<double>& vec) {
    std::cout << "[ ";
    for (std::size_t i = 0; i < vec.size(); ++i) {
        std::cout << std::setprecision(17) << vec[i]
                  << (i + 1 == vec.size() ? " " : ", ");
    }
    std::cout << "]\n";
}

std::vector<double> linspace(double a, double b, std::size_t n) {
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

double safe_step(double value, double scale_hint) {
    double a = std::fabs(value);
    double s = std::fabs(scale_hint);

    double step = 0.0;
    if (std::isfinite(s) && s > 0.0) step = 0.05 * s;
    if (std::isfinite(a) && a > 0.0) step = std::max(step, 0.01 * a);
    if (!std::isfinite(step) || step <= 0.0) step = 1e-3;
    return step;
}

// -----------------------------------------------------------------------------
// CSV export
// -----------------------------------------------------------------------------

class CsvExporter {
public:
    static void save_bestfit(const std::string& path,
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

    static void save_contours(const std::string& path,
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

    static void save_grid(const std::string& path,
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
};

// -----------------------------------------------------------------------------
// Minuit core
// -----------------------------------------------------------------------------

struct ParamLimit {
    std::size_t idx;
    double low;
    double high;
};

struct MinuitFitOptions {
    double up = 0.5;               // NLL => 0.5 for 1D 1σ
    unsigned strategy = 2;
    unsigned max_fcn = 100000;
    double tolerance = 0.2;
    bool run_hesse = true;
    unsigned hesse_maxcalls = 0;
    bool verbose = true;
};

struct MinuitJointFit {
    Vector x_hat;
    std::vector<double> x_err;
    std::vector<double> cov_eigs;
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

class MinuitRunner {
public:
    static MinuitJointFit fit(const std::function<double(const std::vector<double>&)>& f,
                              const std::vector<std::string>& names,
                              const std::vector<double>& x0,
                              const std::vector<double>& scale_hints,
                              const std::vector<ParamLimit>& limits,
                              const MinuitFitOptions& opt) {
        if (x0.size() != names.size() || x0.size() != scale_hints.size()) {
            throw std::invalid_argument("MinuitRunner::fit: names/x0/scale_hints size mismatch");
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
            log_summary("MIGRAD+HESSE", min);
        }

        return extract_result(min, names, x0.size(), opt.verbose);
    }

    static void log_summary(const std::string& tag, const M2::FunctionMinimum& min) {
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

private:
    static MinuitJointFit extract_result(const M2::FunctionMinimum& min,
                                         const std::vector<std::string>& names,
                                         std::size_t n,
                                         bool verbose) {
        std::vector<double> eigs;
        double cond = std::numeric_limits<double>::infinity();

        if (min.HasValidCovariance()) {
            M2::MnEigen eigen;
            eigs = eigen(min.UserState().Covariance());

            double min_pos = std::numeric_limits<double>::infinity();
            double max_pos = 0.0;
            for (double e : eigs) {
                if (std::isfinite(e) && e > 0.0) {
                    min_pos = std::min(min_pos, e);
                    max_pos = std::max(max_pos, e);
                }
            }
            if (min_pos < std::numeric_limits<double>::infinity() && max_pos > 0.0) {
                cond = max_pos / min_pos;
            }

            if (verbose && !eigs.empty()) {
                std::cout << "Cov eigen min/max = "
                          << std::setprecision(6) << eigs.front()
                          << " / " << eigs.back()
                          << "   (cond ~ " << cond << ")\n";
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
};

// -----------------------------------------------------------------------------
// Problem / fit domain objects
// -----------------------------------------------------------------------------

struct JointLikelihoodProblem {
    std::vector<std::string> names;
    std::vector<double> x0;
    std::vector<double> scale_hints;
    std::vector<ParamLimit> limits;
    std::function<double(const std::vector<double>&)> f_joint;
};

struct JointFitOutput {
    FitResult fr;
    MinuitJointFit mj;
    JointLikelihoodProblem problem;
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

    JointFitOutput fit_joint_with_minuit(const std::vector<ParamId>& p_ids,
                                         const std::vector<ParamId>& eta_ids,
                                         const Vector& p0) const {
        const std::size_t p_dim = p0.size();
        Vector eta0 = like_.nuisance_central_values;
        Vector eta_scales = like_.nuisance_dist->get_stds();

        if (eta0.size() != eta_scales.size()) {
            throw std::runtime_error("eta central values and eta stds do not have same size");
        }

        JointLikelihoodProblem problem;
        problem.x0.reserve(p_dim + eta0.size());
        problem.x0.insert(problem.x0.end(), p0.begin(), p0.end());
        problem.x0.insert(problem.x0.end(), eta0.begin(), eta0.end());

        problem.scale_hints.reserve(problem.x0.size());
        for (std::size_t i = 0; i < p_dim; ++i) {
            double hint = std::fabs(p0[i]);
            if (hint < 1e-3) hint = 0.01;
            problem.scale_hints.push_back(hint);
        }
        for (double s : eta_scales) {
            problem.scale_hints.push_back(std::max(1e-12, std::fabs(s)));
        }

        problem.names.reserve(problem.x0.size());
        for (const auto& pid : p_ids) problem.names.push_back(to_string_any(pid));
        for (const auto& pid : eta_ids) problem.names.push_back(to_string_any(pid));

        for (std::size_t i = 0; i < p_dim; ++i) {
            if (problem.names[i].find("FCONST") != std::string::npos) {
                problem.limits.push_back(ParamLimit{i, 0.05, 0.35});
            }
        }

        for (std::size_t i = p_dim; i < problem.names.size(); ++i) {
            const std::string& nm = problem.names[i];
            const double c = problem.x0[i];
            const double s = std::max(1e-12, std::fabs(problem.scale_hints[i]));

            if (nm.find("SMINPUTS:3") != std::string::npos) {
                problem.limits.push_back(ParamLimit{i, 0.05, 0.30});
            } else if (nm.find("MASS:") != std::string::npos ||
                       nm.find("FLIFE:") != std::string::npos ||
                       nm.find("FCONST:") != std::string::npos ||
                       nm.find("FMASS:") != std::string::npos ||
                       nm.find("SMINPUTS:5") != std::string::npos ||
                       nm.find("SMINPUTS:6") != std::string::npos) {
                problem.limits.push_back(ParamLimit{i, std::max(1e-12, c - 5.0 * s), c + 5.0 * s});
            }
        }

        problem.f_joint = make_joint_f(p_dim);

        MinuitFitOptions opt;
        opt.up = 0.5;
        opt.strategy = strategy_;
        opt.max_fcn = static_cast<unsigned>(max_fcn_);
        opt.tolerance = tolerance_;
        opt.run_hesse = true;
        opt.verbose = true;

        MinuitJointFit mj = MinuitRunner::fit(problem.f_joint,
                                              problem.names,
                                              problem.x0,
                                              problem.scale_hints,
                                              problem.limits,
                                              opt);

        FitResult fr;
        fr.ell_hat = mj.fmin;
        fr.p_hat.assign(mj.x_hat.begin(), mj.x_hat.begin() + p_dim);
        fr.eta_hat.assign(mj.x_hat.begin() + p_dim, mj.x_hat.end());
        fr.p_hat_std.assign(p_dim, 0.0);
        fr.p_hat_correlations = RealMatrix(p_dim, p_dim);

        if (mj.has_valid_covar) {
            for (std::size_t i = 0; i < p_dim; ++i) {
                fr.p_hat_std[i] = std::sqrt(std::max(0.0, mj.cov.at(i, i)));
            }

            for (std::size_t i = 0; i < p_dim; ++i) {
                for (std::size_t j = 0; j < p_dim; ++j) {
                    const double di = std::sqrt(std::max(0.0, mj.cov.at(i, i)));
                    const double dj = std::sqrt(std::max(0.0, mj.cov.at(j, j)));
                    fr.p_hat_correlations.at(i, j) =
                        (di > 0.0 && dj > 0.0) ? (mj.cov.at(i, j) / (di * dj)) : 0.0;
                }
            }
        } else {
            for (std::size_t i = 0; i < p_dim && i < mj.x_err.size(); ++i) {
                fr.p_hat_std[i] = mj.x_err[i];
                fr.p_hat_correlations.at(i, i) = 1.0;
            }
        }

        return JointFitOutput{fr, mj, problem};
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
// Contour strategies
// -----------------------------------------------------------------------------

struct ContourComputationInput {
    std::string xname;
    std::string yname;
    unsigned px = 0;
    unsigned py = 1;
    double best_fval = 0.0;
    Vector p_hat;
    Vector p_std;
    JointLikelihoodProblem problem;
    MinuitJointFit best_fit;
};

struct ContourComputationResult {
    bool success = false;
    bool used_grid = false;
    std::vector<std::pair<double, double>> c68;
    std::vector<std::pair<double, double>> c95;
    std::vector<double> xs;
    std::vector<double> ys;
    std::vector<double> z;
};

struct ProfileXYResult {
    double fmin = 1e300;
    std::vector<double> x_hat;
    bool ok = false;
};

static ProfileXYResult profiled_fit_at_fixed_xy(
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
    unsigned max_fcn,
    double tolerance
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

    M2::FunctionMinimum min = migrad(max_fcn, tolerance);

    ProfileXYResult out;
    out.ok = min.IsValid();
    out.fmin = out.ok ? min.Fval() : 1e300;
    out.x_hat = x_start;

    if (out.ok) {
        const auto& st = min.UserState();
        for (std::size_t i = 0; i < x_start.size(); ++i) {
            out.x_hat[i] = st.Value(names[i].c_str());
        }
    } else {
        out.x_hat[px] = xval;
        out.x_hat[py] = yval;
        out.fmin = f_joint(out.x_hat);
    }

    return out;
}

class IContourStrategy {
public:
    virtual ~IContourStrategy() = default;
    virtual ContourComputationResult compute(const ContourComputationInput& input) const = 0;
};

class MnContoursStrategy final : public IContourStrategy {
public:
    struct Options {
        unsigned npoints = 80;
        unsigned refit_max_fcn = 30000;
        double tolerance = 0.2;
        unsigned strategy = 2;
    };

    MnContoursStrategy() = default;
    explicit MnContoursStrategy(const Options& options) : opt_(options) {}

    ContourComputationResult compute(const ContourComputationInput& input) const override {
        ContourComputationResult result;
        if (!input.best_fit.ok || !input.best_fit.min) {
            return result;
        }

        const bool ok68 = compute_single(input, 2.30 / 2.0, result.c68);
        const bool ok95 = compute_single(input, 5.99 / 2.0, result.c95);

        result.success = ok68 && ok95;
        result.used_grid = false;
        return result;
    }

private:
    bool compute_single(const ContourComputationInput& input,
                        double up_contour,
                        std::vector<std::pair<double, double>>& out_points) const {
        std::vector<double> refit_scales = input.problem.scale_hints;
        const double scale_factor = std::sqrt(up_contour / 0.5);

        for (std::size_t i = 0; i < refit_scales.size() && i < input.best_fit.x_err.size(); ++i) {
            const double err = std::fabs(input.best_fit.x_err[i]);
            if (std::isfinite(err) && err > 0.0) {
                refit_scales[i] = std::max(refit_scales[i], err * scale_factor);
            }
        }

        MinuitFitOptions refit_opt;
        refit_opt.up = up_contour;
        refit_opt.strategy = opt_.strategy;
        refit_opt.max_fcn = opt_.refit_max_fcn;
        refit_opt.tolerance = opt_.tolerance;
        refit_opt.run_hesse = true;
        refit_opt.verbose = false;

        MinuitJointFit refit = MinuitRunner::fit(input.problem.f_joint,
                                                 input.problem.names,
                                                 input.best_fit.x_hat,
                                                 refit_scales,
                                                 input.problem.limits,
                                                 refit_opt);

        if (!refit.ok || !refit.min) {
            return false;
        }

        try {
            GenericFCN fcn(input.problem.f_joint, up_contour);
            M2::MnContours contours(fcn, *refit.min, 2);
            out_points = contours(input.px, input.py, opt_.npoints);
            return out_points.size() >= 4;
        } catch (...) {
            out_points.clear();
            return false;
        }
    }

    Options opt_{};
};

class GridProfileContourStrategy final : public IContourStrategy {
public:
    struct Options {
        std::size_t nx = 31;
        std::size_t ny = 31;
        unsigned strategy = 1;
        unsigned max_fcn = 1200;
        double tolerance = 0.5;
        double n_sigma_window = 4.0;
        double hard_low = 0.05;
        double hard_high = 0.35;
    };

    GridProfileContourStrategy() = default;
    explicit GridProfileContourStrategy(const Options& options) : opt_(options) {}

    ContourComputationResult compute(const ContourComputationInput& input) const override {
        ContourComputationResult result;
        result.used_grid = true;

        double x0 = input.p_hat.at(0);
        double y0 = input.p_hat.at(1);
        double sx = std::max(0.01, input.p_std.at(0));
        double sy = std::max(0.01, input.p_std.at(1));

        double xlo = std::max(opt_.hard_low, x0 - opt_.n_sigma_window * sx);
        double xhi = std::min(opt_.hard_high, x0 + opt_.n_sigma_window * sx);
        double ylo = std::max(opt_.hard_low, y0 - opt_.n_sigma_window * sy);
        double yhi = std::min(opt_.hard_high, y0 + opt_.n_sigma_window * sy);

        if (!(xhi > xlo)) { xlo = 0.10; xhi = 0.30; }
        if (!(yhi > ylo)) { ylo = 0.10; yhi = 0.30; }

        result.xs = linspace(xlo, xhi, opt_.nx);
        result.ys = linspace(ylo, yhi, opt_.ny);
        result.z.assign(opt_.nx * opt_.ny, 1e300);

        std::vector<double> seed = input.best_fit.x_hat;

        for (std::size_t iy = 0; iy < opt_.ny; ++iy) {
            const bool reverse = (iy % 2 == 1);

            if (!reverse) {
                for (std::size_t ix = 0; ix < opt_.nx; ++ix) {
                    auto pr = profile_at_fixed_xy(input,
                                                  seed,
                                                  result.xs[ix],
                                                  result.ys[iy]);
                    result.z[iy * opt_.nx + ix] = pr.fmin - input.best_fval;
                    seed = pr.x_hat;
                }
            } else {
                for (std::size_t k = 0; k < opt_.nx; ++k) {
                    std::size_t ix = opt_.nx - 1 - k;
                    auto pr = profile_at_fixed_xy(input,
                                                  seed,
                                                  result.xs[ix],
                                                  result.ys[iy]);
                    result.z[iy * opt_.nx + ix] = pr.fmin - input.best_fval;
                    seed = pr.x_hat;
                }
            }
        }

        result.success = true;
        return result;
    }

private:
    ProfileXYResult profile_at_fixed_xy(const ContourComputationInput& input,
                                        const std::vector<double>& seed,
                                        double xval,
                                        double yval) const {
        return profiled_fit_at_fixed_xy(input.problem.f_joint,
                                        input.problem.names,
                                        seed,
                                        input.problem.scale_hints,
                                        input.problem.limits,
                                        input.px,
                                        input.py,
                                        xval,
                                        yval,
                                        opt_.strategy,
                                        opt_.max_fcn,
                                        opt_.tolerance);
    }

    Options opt_{};
};

class FallbackContourStrategy final : public IContourStrategy {
public:
    FallbackContourStrategy(std::unique_ptr<IContourStrategy> primary,
                            std::unique_ptr<IContourStrategy> fallback)
        : primary_(std::move(primary)), fallback_(std::move(fallback)) {}

    ContourComputationResult compute(const ContourComputationInput& input) const override {
        ContourComputationResult primary_result = primary_->compute(input);
        if (primary_result.success) {
            return primary_result;
        }

        std::cerr << "[WARN] MnContours failed; fallback to grid scan.\n";
        return fallback_->compute(input);
    }

private:
    std::unique_ptr<IContourStrategy> primary_;
    std::unique_ptr<IContourStrategy> fallback_;
};

// -----------------------------------------------------------------------------
// Application bootstrap
// -----------------------------------------------------------------------------

struct BuiltProblem {
    std::vector<ParamId> p_ids;
    std::vector<ParamId> eta_ids;
    std::vector<ExperimentObs> obs_ids;
    LikelihoodContext ctx;
    std::shared_ptr<ObservableInterfaceAdapterObs> model;
    Vector p0;
};

BuiltProblem build_problem(StatisticManager& stat,
                           const StatisticConfig& config,
                           const std::shared_ptr<ObservableInterfaceAdapterObs>& model) {
    LOG_INFO("fill_cache #1");
    // stat.fill_cache();

    auto start_u = std::chrono::steady_clock::now();
    stat.compute_uncertainties();
    auto stop_u = std::chrono::steady_clock::now();
    auto us_u = std::chrono::duration_cast<std::chrono::microseconds>(stop_u - start_u).count();
    std::cout << "Uncertainty estimation time: " << us_u << " us\n";

    LOG_INFO("fill_cache #2");
    // stat.fill_cache();

    auto p_specs_map = stat.get_p_specs(config.p_specs);
    auto eta_specs_real = stat.get_all_obss_deps();
    for (const auto& [pid, _] : p_specs_map) eta_specs_real.erase(pid);
    auto exp_obs_map = stat.get_obs_exp();

    auto unz_p   = unzip(p_specs_map);
    auto unz_eta = unzip(eta_specs_real);
    auto unz_obs = unzip(exp_obs_map);

    auto nuisance_dist = stat.build_nuisance_distribution();
    auto exp_obs_dist  = stat.build_exp_data_distribution();

    if (nuisance_dist->get_stds().size() != unz_eta.vals.size()) {
        throw std::runtime_error("nuisance std size and eta central size mismatch");
    }
    if (exp_obs_dist->dim() != unz_obs.vals.size()) {
        throw std::runtime_error("exp obs dim and exp obs values size mismatch");
    }

    LikelihoodContext ctx;
    ctx.nuisance_dist = std::move(nuisance_dist);
    ctx.exp_obs_dist  = std::move(exp_obs_dist);
    ctx.nuisance_central_values = unz_eta.vals;
    ctx.exp_obs_values = unz_obs.vals;

    return BuiltProblem{
        unz_p.ids,
        unz_eta.ids,
        unz_obs.ids,
        std::move(ctx),
        model,
        unz_p.vals
    };
}


struct GaussianToyConfig {
    double mu_x = 0.20;
    double mu_y = 0.20;
    double sigma_x = 0.04;
    double sigma_y = 0.04;
    double rho = 0.50;
    double hard_low = 0.05;
    double hard_high = 0.35;
};

double gaussian_toy_nll(const std::vector<double>& x, const GaussianToyConfig& cfg) {
    if (x.size() != 2) {
        throw std::invalid_argument("gaussian_toy_nll expects exactly 2 parameters");
    }

    const double dx = x[0] - cfg.mu_x;
    const double dy = x[1] - cfg.mu_y;

    const double vx = cfg.sigma_x * cfg.sigma_x;
    const double vy = cfg.sigma_y * cfg.sigma_y;
    const double cxy = cfg.rho * cfg.sigma_x * cfg.sigma_y;

    const double det = vx * vy - cxy * cxy;
    if (!(det > 0.0)) {
        throw std::runtime_error("Gaussian covariance is not positive definite");
    }

    const double inv00 =  vy / det;
    const double inv11 =  vx / det;
    const double inv01 = -cxy / det;

    return 0.5 * (inv00 * dx * dx + 2.0 * inv01 * dx * dy + inv11 * dy * dy);
}

RealMatrix gaussian_toy_covariance(const GaussianToyConfig& cfg) {
    RealMatrix cov(2, 2);
    cov.at(0, 0) = cfg.sigma_x * cfg.sigma_x;
    cov.at(0, 1) = cfg.rho * cfg.sigma_x * cfg.sigma_y;
    cov.at(1, 0) = cov.at(0, 1);
    cov.at(1, 1) = cfg.sigma_y * cfg.sigma_y;
    return cov;
}

RealMatrix gaussian_toy_hessian(const GaussianToyConfig& cfg) {
    const double vx = cfg.sigma_x * cfg.sigma_x;
    const double vy = cfg.sigma_y * cfg.sigma_y;
    const double cxy = cfg.rho * cfg.sigma_x * cfg.sigma_y;
    const double det = vx * vy - cxy * cxy;

    RealMatrix h(2, 2);
    h.at(0, 0) =  vy / det;
    h.at(0, 1) = -cxy / det;
    h.at(1, 0) = h.at(0, 1);
    h.at(1, 1) =  vx / det;
    return h;
}

std::vector<std::pair<double, double>>
gaussian_toy_contour_points(const GaussianToyConfig& cfg,
                            double delta_nll,
                            std::size_t npoints = 240) {
    std::vector<std::pair<double, double>> pts;
    pts.reserve(npoints);

    const double scale = std::sqrt(2.0 * delta_nll);
    const double sx = cfg.sigma_x;
    const double sy = cfg.sigma_y;
    const double rho = cfg.rho;
    const double rho_perp = std::sqrt(std::max(0.0, 1.0 - rho * rho));

    constexpr double two_pi = 6.2831853071795864769;

    for (std::size_t i = 0; i < npoints; ++i) {
        const double t = two_pi * double(i) / double(npoints);
        const double u0 = std::cos(t);
        const double u1 = std::sin(t);

        const double x = cfg.mu_x + scale * sx * u0;
        const double y = cfg.mu_y + scale * sy * (rho * u0 + rho_perp * u1);
        pts.emplace_back(x, y);
    }

    return pts;
}

std::vector<double> gaussian_toy_grid(const GaussianToyConfig& cfg,
                                      const std::vector<double>& xs,
                                      const std::vector<double>& ys) {
    std::vector<double> z(xs.size() * ys.size(), 0.0);
    for (std::size_t iy = 0; iy < ys.size(); ++iy) {
        for (std::size_t ix = 0; ix < xs.size(); ++ix) {
            z[iy * xs.size() + ix] = gaussian_toy_nll({xs[ix], ys[iy]}, cfg);
        }
    }
    return z;
}

void save_gaussian_reference(const std::string& path, const GaussianToyConfig& cfg) {
    std::ofstream out(path);
    out << std::setprecision(17);
    out << "key,value\n";
    out << "mu_x," << cfg.mu_x << "\n";
    out << "mu_y," << cfg.mu_y << "\n";
    out << "sigma_x," << cfg.sigma_x << "\n";
    out << "sigma_y," << cfg.sigma_y << "\n";
    out << "rho," << cfg.rho << "\n";
    out << "delta68," << 2.30 / 2.0 << "\n";
    out << "delta95," << 5.99 / 2.0 << "\n";
}

} // namespace contour_app

int main(int argc, char** argv) {
    using namespace contour_app;

    (void)argc;
    (void)argv;

    const GaussianToyConfig cfg;

    std::cout << "=== Gaussian toy validation ===\n";
    std::cout << "Expected best-fit = (" << cfg.mu_x << ", " << cfg.mu_y << ")\n";
    std::cout << "Expected std      = (" << cfg.sigma_x << ", " << cfg.sigma_y << ")\n";
    std::cout << "Expected rho      = " << cfg.rho << "\n";
    std::cout << "Expected covariance:\n" << gaussian_toy_covariance(cfg) << "\n";
    std::cout << "Expected Hessian = covariance^{-1}:\n" << gaussian_toy_hessian(cfg) << "\n";

    JointLikelihoodProblem problem;
    problem.names = {"x", "y"};
    problem.x0 = {0.12, 0.28};
    problem.scale_hints = {cfg.sigma_x, cfg.sigma_y};
    problem.limits = {
        ParamLimit{0, cfg.hard_low, cfg.hard_high},
        ParamLimit{1, cfg.hard_low, cfg.hard_high}
    };
    problem.f_joint = [cfg](const std::vector<double>& x) {
        return gaussian_toy_nll(x, cfg);
    };

    MinuitFitOptions opt;
    opt.up = 0.5;
    opt.strategy = 2;
    opt.max_fcn = 5000;
    opt.tolerance = 0.1;
    opt.run_hesse = true;
    opt.verbose = true;

    auto start_m = std::chrono::steady_clock::now();
    MinuitJointFit mj = MinuitRunner::fit(problem.f_joint,
                                          problem.names,
                                          problem.x0,
                                          problem.scale_hints,
                                          problem.limits,
                                          opt);
    auto stop_m = std::chrono::steady_clock::now();

    auto us_m = std::chrono::duration_cast<std::chrono::microseconds>(stop_m - start_m).count();
    std::cout << "\nMLE (Minuit) fitting time: " << us_m << " us\n";

    Vector p_hat = mj.x_hat;
    Vector p_std(2);
    RealMatrix p_corr(2, 2);

    if (mj.has_valid_covar) {
        for (std::size_t i = 0; i < 2; ++i) {
            p_std[i] = std::sqrt(std::max(0.0, mj.cov.at(i, i)));
        }
        for (std::size_t i = 0; i < 2; ++i) {
            for (std::size_t j = 0; j < 2; ++j) {
                const double di = std::sqrt(std::max(0.0, mj.cov.at(i, i)));
                const double dj = std::sqrt(std::max(0.0, mj.cov.at(j, j)));
                p_corr.at(i, j) = (di > 0.0 && dj > 0.0) ? mj.cov.at(i, j) / (di * dj) : 0.0;
            }
        }
    } else {
        p_std[0] = mj.x_err.size() > 0 ? mj.x_err[0] : 0.0;
        p_std[1] = mj.x_err.size() > 1 ? mj.x_err[1] : 0.0;
        p_corr.at(0, 0) = 1.0;
        p_corr.at(1, 1) = 1.0;
        p_corr.at(0, 1) = 0.0;
        p_corr.at(1, 0) = 0.0;
    }

    std::cout << "ell_hat = " << std::setprecision(17) << mj.fmin << "\n";
    std::cout << "p_hat = "; print_vec(p_hat);
    std::cout << "p_hat_std = "; print_vec(p_std);
    std::cout << "p_hat_correlations:\n" << p_corr << "\n";

    std::cout << "\nAbsolute differences wrt expected:\n";
    std::cout << "|x_hat - mu_x| = " << std::fabs(p_hat[0] - cfg.mu_x) << "\n";
    std::cout << "|y_hat - mu_y| = " << std::fabs(p_hat[1] - cfg.mu_y) << "\n";
    std::cout << "|sigma_x_fit - sigma_x| = " << std::fabs(p_std[0] - cfg.sigma_x) << "\n";
    std::cout << "|sigma_y_fit - sigma_y| = " << std::fabs(p_std[1] - cfg.sigma_y) << "\n";
    std::cout << "|rho_fit - rho| = " << std::fabs(p_corr.at(0, 1) - cfg.rho) << "\n";
    std::cout << "|sigma_x_fit - sigma_y_fit| = " << std::fabs(p_std[0] - p_std[1]) << "\n";
    std::cout << "|x_hat - y_hat| = " << std::fabs(p_hat[0] - p_hat[1]) << "\n";

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

    CsvExporter::save_bestfit("bestfit.csv", problem.names, p_hat, p_std);
    save_gaussian_reference("gaussian_reference.csv", cfg);
    std::cout << "[INFO] Wrote bestfit.csv\n";
    std::cout << "[INFO] Wrote gaussian_reference.csv\n";

    ContourComputationInput contour_input;
    contour_input.xname = problem.names[0];
    contour_input.yname = problem.names[1];
    contour_input.px = 0;
    contour_input.py = 1;
    contour_input.best_fval = mj.fmin;
    contour_input.p_hat = p_hat;
    contour_input.p_std = p_std;
    contour_input.problem = problem;
    contour_input.best_fit = mj;

    MnContoursStrategy::Options native_opt;
    native_opt.npoints = 120;
    native_opt.refit_max_fcn = 10000;
    native_opt.tolerance = 0.1;
    native_opt.strategy = 2;

    MnContoursStrategy native_strategy(native_opt);
    ContourComputationResult native_result = native_strategy.compute(contour_input);

    if (native_result.success) {
        CsvExporter::save_contours("contours.csv",
                                   contour_input.xname,
                                   contour_input.yname,
                                   native_result.c68,
                                   native_result.c95);
        CsvExporter::save_contours("contours_native.csv",
                                   contour_input.xname,
                                   contour_input.yname,
                                   native_result.c68,
                                   native_result.c95);
        std::cout << "[INFO] Wrote contours.csv and contours_native.csv\n";
    } else {
        std::cerr << "[WARN] Native MnContours failed on Gaussian toy.\n";
    }

    GridProfileContourStrategy::Options grid_opt;
    grid_opt.nx = 61;
    grid_opt.ny = 61;
    grid_opt.strategy = 1;
    grid_opt.max_fcn = 1500;
    grid_opt.tolerance = 0.3;
    grid_opt.n_sigma_window = 4.0;
    grid_opt.hard_low = cfg.hard_low;
    grid_opt.hard_high = cfg.hard_high;

    GridProfileContourStrategy grid_strategy(grid_opt);
    ContourComputationResult grid_result = grid_strategy.compute(contour_input);

    if (!grid_result.success) {
        std::cerr << "[ERROR] Grid computation failed.\n";
        return 6;
    }

    CsvExporter::save_grid("grid.csv",
                           contour_input.xname,
                           contour_input.yname,
                           grid_result.xs,
                           grid_result.ys,
                           grid_result.z);
    std::cout << "[INFO] Wrote grid.csv\n";

    const std::vector<double> z_expected = gaussian_toy_grid(cfg, grid_result.xs, grid_result.ys);
    CsvExporter::save_grid("grid_expected.csv",
                           contour_input.xname,
                           contour_input.yname,
                           grid_result.xs,
                           grid_result.ys,
                           z_expected);
    std::cout << "[INFO] Wrote grid_expected.csv\n";

    const auto expected_c68 = gaussian_toy_contour_points(cfg, 2.30 / 2.0);
    const auto expected_c95 = gaussian_toy_contour_points(cfg, 5.99 / 2.0);
    CsvExporter::save_contours("contours_expected.csv",
                               contour_input.xname,
                               contour_input.yname,
                               expected_c68,
                               expected_c95);
    std::cout << "[INFO] Wrote contours_expected.csv\n";

    double max_abs_grid_diff = 0.0;
    for (std::size_t i = 0; i < grid_result.z.size(); ++i) {
        max_abs_grid_diff = std::max(max_abs_grid_diff, std::fabs(grid_result.z[i] - z_expected[i]));
    }
    std::cout << "max |grid_observed - grid_expected| = " << max_abs_grid_diff << "\n";

    return 0;
}
