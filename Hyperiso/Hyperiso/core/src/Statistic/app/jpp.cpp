#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "StatisticManager.h"
#include "ObservableInterface.h"
#include "ObservableInterfaceAdapter2.h"
#include "StatCorrelationProxy.h"
#include "StatParameterProxy.h"
#include "StatParamSourcesProxy.h"
#include "StatDependencyPruner.h"
#include "Fit.h"
#include "BaseLikelihood.h"
#include "Profiler.h"
#include "IProfilingStrategy.h"
#include "FitAbstraction.h"
#include "NuisanceReader.h"
#include "DefaultNuisancePathsProvider.h"

namespace fit_app {

struct ParamLimit {
    std::size_t idx;
    double low;
    double high;
};

struct BuiltProblem {
    std::vector<ParamId> p_ids;
    std::vector<ParamId> eta_ids;
    std::vector<ExperimentObs> obs_ids;

    std::shared_ptr<LikelihoodContext> ctx;
    std::shared_ptr<BaseLikelihood> like;

    std::vector<std::string> names;
    std::vector<double> x0;
    std::vector<double> scale_hints;
    std::vector<ParamLimit> limits;
};

static std::string param_name(const ParamId& pid) {
    std::ostringstream oss;
    oss << pid;
    return oss.str();
}

static fit_app::ParameterDefinition make_fit_param_def_local(
    const ParamId& pid, double value, double sigma_hint)
{
    fit_app::ParameterDefinition out;
    out.name = param_name(pid);
    out.value = value;
    out.step_hint = (std::isfinite(sigma_hint) && sigma_hint > 0.0)
        ? sigma_hint
        : std::max(1e-3, 0.01 * std::abs(value));

    if (out.name.find("FCONST") != std::string::npos) {
        out.limits = std::make_pair(0.05, 0.35);
    }

    return out;
}

static fit_app::ParameterDefinition make_nuisance_param_def_local(
    const ParamId& pid, double value, double sigma_hint)
{
    fit_app::ParameterDefinition out;
    out.name = param_name(pid);
    out.value = value;

    const double s = (std::isfinite(sigma_hint) && sigma_hint > 0.0)
        ? std::abs(sigma_hint)
        : std::max(1e-3, 0.01 * std::abs(value));

    out.step_hint = s;

    const std::string& nm = out.name;
    if (nm.find("SMINPUTS:3") != std::string::npos) {
        out.limits = std::make_pair(0.05, 0.30);
    } else if (nm.find("MASS:") != std::string::npos ||
               nm.find("FLIFE:") != std::string::npos ||
               nm.find("FCONST:") != std::string::npos ||
               nm.find("FMASS:") != std::string::npos ||
               nm.find("SMINPUTS:5") != std::string::npos ||
               nm.find("SMINPUTS:6") != std::string::npos) {
        out.limits = std::make_pair(std::max(1e-12, value - 5.0 * s), value + 5.0 * s);
    }

    return out;
}

static std::vector<double> ordered_prediction_vector(
    const std::vector<ExperimentObs>& obs_ids,
    const std::map<ObservableId, std::vector<ObservableValue>>& pred_map)
{
    std::vector<double> out;
    out.reserve(obs_ids.size());

    for (const auto& bid : obs_ids) {
        const auto& vec = pred_map.at(bid.obs.s);

        auto it = std::find_if(vec.begin(), vec.end(), [&](const ObservableValue& ov) {
            auto bin = ov.bin.value_or(std::pair<double, double>{0.0, 0.0});
            return bin == bid.obs.p;
        });

        if (it == vec.end()) {
            throw std::runtime_error("Missing predicted observable/bin.");
        }

        out.push_back(it->value);
    }

    return out;
}

static std::vector<fit_app::ParameterDefinition> make_parameter_definitions(
    const std::vector<std::string>& names,
    const std::vector<double>& values,
    const std::vector<double>& scale_hints,
    const std::vector<ParamLimit>& limits)
{
    if (names.size() != values.size() || names.size() != scale_hints.size()) {
        throw std::invalid_argument("names/values/scale_hints size mismatch");
    }

    std::vector<fit_app::ParameterDefinition> parameters;
    parameters.reserve(names.size());

    for (std::size_t i = 0; i < names.size(); ++i) {
        fit_app::ParameterDefinition p;
        p.name = names[i];
        p.value = values[i];
        p.step_hint = scale_hints[i];
        parameters.push_back(std::move(p));
    }

    for (const auto& lim : limits) {
        if (lim.idx < parameters.size()) {
            parameters[lim.idx].limits = std::make_pair(lim.low, lim.high);
        }
    }

    return parameters;
}

static BuiltProblem build_problem(
    StatisticManager& stat,
    const StatisticConfig& config,
    const std::shared_ptr<ObservableInterfaceAdapterObs>& model, std::vector<ParamId> p_specs)
{
    stat.compute_uncertainties();
    stat.update_cache(p_specs);

    auto p_specs_map = stat.get_p_specs(p_specs);
    auto eta_specs_real = stat.get_all_obss_deps();
    for (const auto& [pid, _] : p_specs_map) eta_specs_real.erase(pid);
    auto exp_obs_map = stat.get_obs_exp();

    auto unz_p   = unzip(p_specs_map);
    auto unz_eta = unzip(eta_specs_real);
    auto unz_obs = unzip(exp_obs_map);

    auto ctx = std::make_shared<LikelihoodContext>();
    ctx->nuisance_dist = stat.build_nuisance_distribution();
    ctx->exp_obs_dist  = stat.build_exp_data_distribution();
    ctx->exp_obs_values = unz_obs.vals;

    // Définitions paramètres: mêmes limites que le manager actuel
    ctx->fp_defs.reserve(unz_p.ids.size());
    for (std::size_t i = 0; i < unz_p.ids.size(); ++i) {
        double sigma_hint = std::max(1e-3, std::abs(unz_p.vals[i]));
        ctx->fp_defs.emplace_back(make_fit_param_def_local(unz_p.ids[i], unz_p.vals[i], sigma_hint));
    }

    const auto eta_stds = ctx->nuisance_dist->get_stds();
    ctx->nuis_defs.reserve(unz_eta.ids.size());
    for (std::size_t i = 0; i < unz_eta.ids.size(); ++i) {
        double sigma_hint = (i < eta_stds.size()) ? eta_stds[i] : 1e-3;
        ctx->nuis_defs.emplace_back(make_nuisance_param_def_local(unz_eta.ids[i], unz_eta.vals[i], sigma_hint));
    }

    auto model_fn = [model, obs_ids = unz_obs.ids, p_ids = unz_p.ids, eta_ids = unz_eta.ids]
                    (const std::vector<double>& p_vec, const std::vector<double>& eta_vec) -> std::vector<double> {
        auto pred_map = model->predict_optimized(zip(p_ids, p_vec), zip(eta_ids, eta_vec));
        return ordered_prediction_vector(obs_ids, pred_map);
    };

    auto like = std::make_shared<BaseLikelihood>(model_fn, ctx, unz_p.ids.size());

    BuiltProblem bp;
    bp.p_ids = unz_p.ids;
    bp.eta_ids = unz_eta.ids;
    bp.obs_ids = unz_obs.ids;
    bp.ctx = ctx;
    bp.like = like;

    // vieux style: noms / x0 / scale_hints / limits
    bp.names.reserve(unz_p.ids.size() + unz_eta.ids.size());
    bp.x0.reserve(unz_p.vals.size() + unz_eta.vals.size());
    bp.scale_hints.reserve(unz_p.vals.size() + unz_eta.vals.size());

    for (std::size_t i = 0; i < unz_p.ids.size(); ++i) {
        bp.names.push_back(to_string_any(unz_p.ids[i]));
        bp.x0.push_back(unz_p.vals[i]);
        bp.scale_hints.push_back(std::max(1e-3, std::abs(unz_p.vals[i])));
        if (bp.names.back().find("FCONST") != std::string::npos) {
            bp.limits.push_back({i, 0.05, 0.35});
        }
    }

    for (std::size_t i = 0; i < unz_eta.ids.size(); ++i) {
        const std::size_t idx = unz_p.ids.size() + i;
        const double c = unz_eta.vals[i];
        const double s = (i < eta_stds.size()) ? std::max(1e-12, std::abs(eta_stds[i])) : 1e-3;

        bp.names.push_back(to_string_any(unz_eta.ids[i]));
        bp.x0.push_back(c);
        bp.scale_hints.push_back(s);

        const std::string& nm = bp.names.back();
        if (nm.find("SMINPUTS:3") != std::string::npos) {
            bp.limits.push_back({idx, 0.05, 0.30});
        } else if (nm.find("MASS:") != std::string::npos ||
                   nm.find("FLIFE:") != std::string::npos ||
                   nm.find("FCONST:") != std::string::npos ||
                   nm.find("FMASS:") != std::string::npos ||
                   nm.find("SMINPUTS:5") != std::string::npos ||
                   nm.find("SMINPUTS:6") != std::string::npos) {
            bp.limits.push_back({idx, std::max(1e-12, c - 5.0 * s), c + 5.0 * s});
        }
    }

    return bp;
}

struct OldProfileResult {
    bool ok = false;
    double fmin = 1e300;
    std::vector<double> x_hat;
    fit_app::BackendFitResult raw;
};

static OldProfileResult old_profile_at_fixed_xy(
    const fit_app::IFitBackend& backend,
    const std::shared_ptr<BaseLikelihood>& like,
    const std::vector<std::string>& names,
    const std::vector<double>& x_start,
    const std::vector<double>& scale_hints,
    const std::vector<ParamLimit>& limits,
    unsigned px,
    unsigned py,
    double xval,
    double yval)
{
    const auto parameters = make_parameter_definitions(names, x_start, scale_hints, limits);

    auto objective = fit_app::LambdaObjectiveFunction(
        [like](const std::vector<double>& theta) {
            return like->nll(theta);
        },
        0.5
    );

    fit_app::FitOptions opt;
    opt.up = 0.5;
    opt.strategy = 2;
    opt.max_fcn = 30000;
    opt.tolerance = 0.2;
    opt.run_hesse = false;
    opt.verbose = true;

    auto raw = backend.minimize_with_fixed(objective, parameters, opt, {px, py}, {xval, yval});

    OldProfileResult out;
    out.ok = raw.diagnostics.ok;
    out.raw = raw;
    out.x_hat = x_start;

    if (raw.diagnostics.ok && raw.values.size() == x_start.size()) {
        out.fmin = raw.diagnostics.fmin;
        out.x_hat = raw.values;
    } else {
        out.x_hat[px] = xval;
        out.x_hat[py] = yval;
        out.fmin = like->nll(out.x_hat); // fallback fini comme dans l'ancien esprit
    }

    return out;
}

static void print_theta(const std::string& label, const std::vector<double>& x) {
    std::cout << label << " = (";
    for (std::size_t i = 0; i < x.size(); ++i) {
        std::cout << std::setprecision(17) << x[i];
        if (i + 1 != x.size()) std::cout << ", ";
    }
    std::cout << ")\n";
}

} // namespace fit_app

int main() {
    using namespace fit_app;

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
    std::vector<ParamId> p_specs = {
        ParamId{ParameterType::FLAVOR, "FCONST", {511, 1}},
        ParamId{ParameterType::FLAVOR, "FCONST", {531, 1}}
    };

    auto model = std::make_shared<ObservableInterfaceAdapterObs>(oint);

    std::shared_ptr<INuisancePathsProvider> npp = std::make_shared<DefaultNuisancePathsProvider>();

    StatisticManager stat(
        config,
        model,
        std::make_shared<StatCorrelationProxy>(),
        std::make_shared<StatParameterProxy>(),
        std::make_shared<StatParamSourcesProxy>(),
        std::make_shared<StatDependencyPruner>(),
        std::make_shared<NuisanceReader>(npp)
    );

    BuiltProblem bp = build_problem(stat, config, model, p_specs);

    // MLE jointe actuelle
    auto fitter = std::make_shared<MLFitter>(bp.ctx, [model, obs_ids = bp.obs_ids, p_ids = bp.p_ids, eta_ids = bp.eta_ids]
        (const std::vector<double>& p_vec, const std::vector<double>& eta_vec) -> std::vector<double> {
            auto pred_map = model->predict_optimized(zip(p_ids, p_vec), zip(eta_ids, eta_vec));
            return ordered_prediction_vector(obs_ids, pred_map);
        });

    auto fit = fitter->maximum_likelihood_fit(std::vector<double>(bp.x0.begin(), bp.x0.begin() + bp.p_ids.size()));

    std::vector<double> theta_hat = fit.p_hat;
    theta_hat.insert(theta_hat.end(), fit.eta_hat.begin(), fit.eta_hat.end());

    std::cout << "\n=== JOINT MLE ===\n";
    print_theta("theta_hat", theta_hat);
    std::cout << "ell_hat = " << std::setprecision(17) << fit.ell_hat << "\n";

    // Points à comparer
    std::vector<std::pair<double,double>> test_points = {
        {fit.p_hat[0], fit.p_hat[1]},
        {0.20, 0.20},
        {0.125, 0.125},
        {0.05, 0.20},
        {0.35, 0.35}
    };

    

    auto backend = make_minuit_backend();
    Profiler profiler(make_minuit_backend());

    
    SliceProfilingStrategy slice(0, 1, fit);
    
    std::map<std::size_t, double> current = slice.init_warm_start();

    for (const auto& [x, y] : test_points) {
        std::cout << "\n====================================================\n";
        std::cout << "Testing point (x, y) = (" << x << ", " << y << ")\n";

        std::map<std::size_t, double> fresh = slice.init_warm_start();
        ProfileRequest pr_fresh = slice.build_request(x, y, fresh);
        auto new_res_fresh = profiler.profile(bp.like, pr_fresh);

        std::cout << "\n[NEW Profiler::profile / fresh global seed]\n";
        std::cout << "ok/converged = " << new_res_fresh.converged << "\n";
        std::cout << "fmin         = " << std::setprecision(17) << new_res_fresh.nll_hat << "\n";

        std::vector<double> theta_prof_fresh(theta_hat.size(), std::numeric_limits<double>::quiet_NaN());
        theta_prof_fresh[0] = x;
        theta_prof_fresh[1] = y;
        for (const auto& [idx, val] : new_res_fresh.theta_hat) {
            theta_prof_fresh[idx] = val;
        }
        print_theta("theta_hat_prof_fresh", theta_prof_fresh);
        
        // Old-style direct profiling
        auto old_res = old_profile_at_fixed_xy(
            *backend,
            bp.like,
            bp.names,
            theta_hat,
            bp.scale_hints,
            bp.limits,
            0, 1,
            x, y
        );

        std::cout << "\n[OLD DIRECT minimize_with_fixed]\n";
        std::cout << "ok    = " << old_res.ok << "\n";
        std::cout << "fmin  = " << std::setprecision(17) << old_res.fmin << "\n";
        std::cout << "nfcn  = " << old_res.raw.diagnostics.nfcn << "\n";
        std::cout << "edm   = " << old_res.raw.diagnostics.edm << "\n";
        print_theta("x_hat", old_res.x_hat);

        // New-style current profiler
        ProfileRequest pr = slice.build_request(x, y, current);
        auto new_res = profiler.profile(bp.like, pr);

        std::cout << "\n[NEW Profiler::profile]\n";
        std::cout << "ok/converged = " << new_res.converged << "\n";
        std::cout << "fmin         = " << std::setprecision(17) << new_res.nll_hat << "\n";

        std::vector<double> theta_prof(theta_hat.size(), std::numeric_limits<double>::quiet_NaN());
        theta_prof[0] = x;
        theta_prof[1] = y;
        for (const auto& [idx, val] : new_res.theta_hat) {
            theta_prof[idx] = val;
        }
        print_theta("theta_hat_prof", theta_prof);

        if (new_res.converged) {
            current = new_res.theta_hat;
        }

        std::cout << "\n[COMPARE]\n";
        std::cout << "delta_f = " << std::setprecision(17) << (new_res.nll_hat - old_res.fmin) << "\n";

        std::vector<double> theta_seed = theta_hat;
        theta_seed[0] = x;
        theta_seed[1] = y;
        double f_seed = bp.like->nll(theta_seed);

        std::cout << "\n[SEED ONLY]\n";
        std::cout << "f_seed = " << std::setprecision(17) << f_seed << "\n";

        std::cout << "\n[COMPARE]\n";
        std::cout << "old_f      = " << old_res.fmin << "\n";
        std::cout << "new_f      = " << new_res.nll_hat << "\n";
        std::cout << "seed_f     = " << f_seed << "\n";
        std::cout << "delta_new_old  = " << (new_res.nll_hat - old_res.fmin) << "\n";
        std::cout << "delta_old_seed = " << (old_res.fmin - f_seed) << "\n";
        std::cout << "delta_new_seed = " << (new_res.nll_hat - f_seed) << "\n";
    }

    return 0;
}