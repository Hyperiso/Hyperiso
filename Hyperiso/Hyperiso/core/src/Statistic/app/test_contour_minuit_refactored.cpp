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
#include "Likelihood.h"

#include "FitAbstraction.h"


namespace contour_app {

// -----------------------------------------------------------------------------
// Small utilities
// -----------------------------------------------------------------------------

void print_vec(const std::vector<double>& vec) {
    std::cout << "[ ";
    for (std::size_t i = 0; i < vec.size(); ++i) {
        std::cout << std::setprecision(17) << vec[i]
                  << (i + 1 == vec.size() ? " " : ", ");
    }
    std::cout << "]\n";
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
// Backend adapter helpers
// -----------------------------------------------------------------------------

struct ParamLimit {
    std::size_t idx;
    double low;
    double high;
};

struct MinuitFitOptions {
    double up = 0.5;
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
    BackendFitResult backend_fit;
};

static std::vector<ParameterDefinition> make_parameter_definitions(
    const std::vector<std::string>& names,
    const std::vector<double>& values,
    const std::vector<double>& scale_hints,
    const std::vector<ParamLimit>& limits
) {
    if (names.size() != values.size() || names.size() != scale_hints.size()) {
        throw std::invalid_argument("make_parameter_definitions: names/values/scale_hints size mismatch");
    }

    std::vector<ParameterDefinition> parameters;
    parameters.reserve(names.size());

    for (std::size_t i = 0; i < names.size(); ++i) {
        ParameterDefinition parameter;
        parameter.name = names[i];
        parameter.value = values[i];
        parameter.step_hint = scale_hints[i];
        parameters.push_back(std::move(parameter));
    }

    for (const auto& lim : limits) {
        if (lim.idx < parameters.size()) {
            parameters[lim.idx].limits = std::make_pair(lim.low, lim.high);
        }
    }

    return parameters;
}

static FitOptions to_backend_fit_options(const MinuitFitOptions& opt) {
    FitOptions out;
    out.up = opt.up;
    out.strategy = opt.strategy;
    out.max_fcn = opt.max_fcn;
    out.tolerance = opt.tolerance;
    out.run_hesse = opt.run_hesse;
    out.hesse_maxcalls = opt.hesse_maxcalls;
    out.verbose = opt.verbose;
    return out;
}

class MinuitRunner {
public:
    explicit MinuitRunner(const IFitBackend& backend) : backend_(backend) {}

    MinuitJointFit fit(const std::function<double(const std::vector<double>&)>& f,
                       const std::vector<std::string>& names,
                       const std::vector<double>& x0,
                       const std::vector<double>& scale_hints,
                       const std::vector<ParamLimit>& limits,
                       const MinuitFitOptions& opt) const {
        const auto parameters = make_parameter_definitions(names, x0, scale_hints, limits);
        const LambdaObjectiveFunction objective(f, opt.up);
        const BackendFitResult backend_fit = backend_.minimize(objective, parameters, to_backend_fit_options(opt));
        return extract_result(backend_fit, opt.verbose);
    }

private:
    static MinuitJointFit extract_result(const BackendFitResult& backend_fit, bool verbose) {
        MinuitJointFit out;
        out.backend_fit = backend_fit;
        out.x_hat = backend_fit.values;
        out.x_err = backend_fit.errors;
        out.cov = backend_fit.covariance;

        out.fmin = backend_fit.diagnostics.fmin;
        out.edm  = backend_fit.diagnostics.edm;
        out.nfcn = backend_fit.diagnostics.nfcn;
        out.ok   = backend_fit.diagnostics.ok;

        out.has_valid_covar    = backend_fit.diagnostics.has_valid_covar;
        out.has_posdef_covar   = backend_fit.diagnostics.has_posdef_covar;
        out.has_accurate_covar = backend_fit.diagnostics.has_accurate_covar;
        out.made_posdef        = backend_fit.diagnostics.made_posdef;
        out.cov_eigs           = backend_fit.diagnostics.cov_eigs;
        out.cond_number        = backend_fit.diagnostics.cond_number;

        if (verbose && !out.cov_eigs.empty()) {
            std::cout << "Cov eigen min/max = "
                      << std::setprecision(6) << out.cov_eigs.front()
                      << " / " << out.cov_eigs.back()
                      << "   (cond ~ " << out.cond_number << ")\n";
        }

        return out;
    }

    const IFitBackend& backend_;
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

    MinuitMLEstimatorLocal(const IFitBackend& backend,
                           LikelihoodContext ctx,
                           ModelFn model,
                           std::size_t max_fcn,
                           double tolerance,
                           unsigned strategy)
        : backend_(backend)
        , like_(std::move(ctx))
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

        MinuitRunner runner(backend_);
        MinuitJointFit mj = runner.fit(problem.f_joint,
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

    const IFitBackend& backend_;
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
    const IFitBackend* backend = nullptr;
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
    const IFitBackend& backend,
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
    const auto parameters = make_parameter_definitions(names, x_start, scale_hints, limits);
    const LambdaObjectiveFunction objective(f_joint, 0.5);

    FitOptions opt;
    opt.up = 0.5;
    opt.strategy = strategy;
    opt.max_fcn = max_fcn;
    opt.tolerance = tolerance;
    opt.run_hesse = false;
    opt.verbose = false;

    const BackendFitResult fit = backend.minimize_with_fixed(objective,
                                                             parameters,
                                                             opt,
                                                             {px, py},
                                                             {xval, yval});

    ProfileXYResult out;
    out.ok = fit.diagnostics.ok;
    out.fmin = out.ok ? fit.diagnostics.fmin : 1e300;
    out.x_hat = x_start;

    if (out.ok && fit.values.size() == x_start.size()) {
        out.x_hat = fit.values;
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
        if (input.backend == nullptr || !input.best_fit.ok || !input.best_fit.backend_fit.state) {
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

        const auto refit_parameters = make_parameter_definitions(input.problem.names,
                                                                 input.best_fit.x_hat,
                                                                 refit_scales,
                                                                 input.problem.limits);

        const LambdaObjectiveFunction objective(input.problem.f_joint, up_contour);

        FitOptions refit_opt;
        refit_opt.up = up_contour;
        refit_opt.strategy = opt_.strategy;
        refit_opt.max_fcn = opt_.refit_max_fcn;
        refit_opt.tolerance = opt_.tolerance;
        refit_opt.run_hesse = true;
        refit_opt.verbose = false;

        const BackendFitResult refit = input.backend->minimize(objective, refit_parameters, refit_opt);
        if (!refit.diagnostics.ok || !refit.state) {
            return false;
        }

        ContourOptions contour_opt;
        contour_opt.up = up_contour;
        contour_opt.npoints = opt_.npoints;
        contour_opt.strategy = opt_.strategy;
        contour_opt.max_fcn = opt_.refit_max_fcn;
        contour_opt.tolerance = opt_.tolerance;

        const BackendContourResult contour = input.backend->contour(objective,
                                                                    refit,
                                                                    input.px,
                                                                    input.py,
                                                                    contour_opt);
        out_points = contour.points;
        return contour.success;
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
        return profiled_fit_at_fixed_xy(*input.backend,
                                        input.problem.f_joint,
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
    std::vector<BinnedObservableId> obs_ids;
    LikelihoodContext ctx;
    std::shared_ptr<ObservableInterfaceAdapterObs> model;
    Vector p0;
};

BuiltProblem build_problem(StatisticManager& stat,
                           const StatisticConfig& config,
                           const std::shared_ptr<ObservableInterfaceAdapterObs>& model) {
    LOG_INFO("fill_cache #1");
    stat.fill_cache();

    auto start_u = std::chrono::steady_clock::now();
    stat.compute_uncertainties();
    auto stop_u = std::chrono::steady_clock::now();
    auto us_u = std::chrono::duration_cast<std::chrono::microseconds>(stop_u - start_u).count();
    std::cout << "Uncertainty estimation time: " << us_u << " us\n";

    LOG_INFO("fill_cache #2");
    stat.fill_cache();

    auto p_specs_map = stat.get_p_specs();
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

} // namespace contour_app

int main(int argc, char** argv) {
    using namespace contour_app;

    HyperisoMaster hyp;
    HyperisoConfig config_hyp;
    config_hyp.model = Model::SM;
    hyp.init("lha/si_input.flha", config_hyp);

    auto oint = std::make_shared<ObservableInterface>();
    oint->add_observable(ObservableMapper::to_id(Observables::BR_BS_MUMU_UNTAG), QCDOrder::LO, true)
        .add_observable(ObservableMapper::to_id(Observables::BR_BD_MUMU), QCDOrder::LO, true);

    StatisticConfig config;
    config.MC_draws = 100;
    config.MLE_max_iter = 120000;
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

    BuiltProblem built = build_problem(stat, config, model);
    std::unique_ptr<IFitBackend> backend = make_minuit_backend();

    auto model_fn = [model, obs_ids = built.obs_ids, p_ids = built.p_ids, eta_ids = built.eta_ids]
                    (const Vec& p_vec, const Vec& eta_vec) -> Vec {
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
    MinuitMLEstimatorLocal est(*backend, std::move(built.ctx), model_fn, config.MLE_max_iter, config.MLE_tol, 2);
    JointFitOutput fit_out = est.fit_joint_with_minuit(built.p_ids, built.eta_ids, built.p0);
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

    std::vector<std::string> p_names;
    for (const auto& pid : built.p_ids) p_names.push_back(to_string_any(pid));
    CsvExporter::save_bestfit("bestfit.csv", p_names, fr.p_hat, fr.p_hat_std);
    std::cout << "[INFO] Wrote bestfit.csv\n";

    if (built.p_ids.size() == 2) {
        ContourComputationInput contour_input;
        contour_input.xname = to_string_any(built.p_ids[0]);
        contour_input.yname = to_string_any(built.p_ids[1]);
        contour_input.px = 0;
        contour_input.py = 1;
        contour_input.best_fval = fr.ell_hat;
        contour_input.p_hat = fr.p_hat;
        contour_input.p_std = fr.p_hat_std;
        contour_input.problem = fit_out.problem;
        contour_input.best_fit = fit_out.mj;
        contour_input.backend = backend.get();

        FallbackContourStrategy contour_strategy(
            std::make_unique<MnContoursStrategy>(),
            std::make_unique<GridProfileContourStrategy>()
        );

        ContourComputationResult contour_result = contour_strategy.compute(contour_input);

        if (!contour_result.success) {
            std::cerr << "[ERROR] Contour computation failed.\n";
            return 6;
        }

        if (contour_result.used_grid) {
            CsvExporter::save_grid("grid.csv",
                                   contour_input.xname,
                                   contour_input.yname,
                                   contour_result.xs,
                                   contour_result.ys,
                                   contour_result.z);
            std::cout << "[INFO] Wrote grid.csv\n";

            double best_grid = 1e300;
            std::size_t best_ix = 0;
            std::size_t best_iy = 0;
            const std::size_t nx = contour_result.xs.size();
            const std::size_t ny = contour_result.ys.size();

            for (std::size_t iy = 0; iy < ny; ++iy) {
                for (std::size_t ix = 0; ix < nx; ++ix) {
                    double val = contour_result.z[iy * nx + ix];
                    if (val < best_grid) {
                        best_grid = val;
                        best_ix = ix;
                        best_iy = iy;
                    }
                }
            }

            std::cout << "[INFO] grid min delta_nll = " << best_grid
                      << " at (" << contour_result.xs[best_ix] << ", "
                      << contour_result.ys[best_iy] << ")\n";
            std::cout << "[INFO] best-fit           = (" << fr.p_hat[0] << ", " << fr.p_hat[1] << ")\n";
        } else {
            CsvExporter::save_contours("contours.csv",
                                       contour_input.xname,
                                       contour_input.yname,
                                       contour_result.c68,
                                       contour_result.c95);
            std::cout << "[INFO] Wrote contours.csv\n";
        }
    }

    return 0;
}
