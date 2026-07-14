
#include "StatisticManager.h"
#include "StatisticCachePrinter.h"
#include "StatisticProgress.h"

#include <algorithm>
#include <type_traits>

std::string param_name(const ParamId& pid) {
    std::ostringstream oss;
    oss << pid;
    return oss.str();
}

static void dump_matrix_sanity(const RealMatrix& M,
                               const std::vector<ParamId>& ids,
                               const std::string& name)
{
    std::cout << "[FIT] Matrix " << name
              << " shape=" << M.rows() << "x" << M.cols()
              << " symmetric=" << M.is_symmetric() << "\n";

    double diag_min = std::numeric_limits<double>::infinity();
    double diag_max = -std::numeric_limits<double>::infinity();
    std::size_t nonpos_diag = 0;

    double max_abs_offdiag = 0.0;
    std::size_t n_abs_gt_1 = 0;
    std::pair<std::size_t,std::size_t> argmax_offdiag{0,0};

    for (std::size_t i = 0; i < M.rows(); ++i) {
        const double d = M.at(i,i);
        diag_min = std::min(diag_min, d);
        diag_max = std::max(diag_max, d);
        if (!(d > 0.0)) ++nonpos_diag;

        for (std::size_t j = i + 1; j < M.cols(); ++j) {
            const double a = std::abs(M.at(i,j));
            if (a > max_abs_offdiag) {
                max_abs_offdiag = a;
                argmax_offdiag = {i,j};
            }
            if (a > 1.0 + 1e-10) ++n_abs_gt_1;
        }
    }

    std::cout << "[FIT] " << name
              << " diag_min=" << diag_min
              << " diag_max=" << diag_max
              << " nonpos_diag=" << nonpos_diag
              << " max_abs_offdiag=" << max_abs_offdiag
              << " n_abs_offdiag_gt_1=" << n_abs_gt_1
              << "\n";

    if (M.rows() == M.cols() && M.is_symmetric()) {
        EigenSystem e = M.eig();

        std::size_t imin = 0;
        double evmin = e.D.at(0,0);
        for (std::size_t i = 1; i < M.rows(); ++i) {
            if (e.D.at(i,i) < evmin) {
                evmin = e.D.at(i,i);
                imin = i;
            }
        }

        std::vector<std::pair<double,std::size_t>> comps;
        for (std::size_t i = 0; i < M.rows(); ++i) {
            comps.push_back({std::abs(e.P.at(i, imin)), i});
        }
        std::sort(comps.begin(), comps.end(),
                  [](const auto& a, const auto& b) { return a.first > b.first; });

        std::cout << "[FIT] " << name << " smallest_eigenvalue=" << evmin << "\n";
        std::cout << "[FIT] " << name << " dominant entries of smallest mode:\n";
        for (std::size_t k = 0; k < std::min<std::size_t>(10, comps.size()); ++k) {
            std::size_t i = comps[k].second;
            std::cout << "[FIT]   " << ids[i]
                      << " coeff=" << e.P.at(i, imin) << "\n";
        }
    }
}

namespace {

std::string marginal_type_name(MarginalType mt) {
    switch (mt) {
        case MarginalType::GAUSSIAN:      return "Gaussian";
        case MarginalType::HALF_GAUSSIAN: return "Split Gaussian";
        case MarginalType::FLAT:          return "Flat";
        case MarginalType::LIKELIHOOD:    return "Likelihood";
    }
    return "Unknown";
}

std::string copula_type_name(CopulaType ct) {
    switch (ct) {
        case CopulaType::GAUSSIAN:  return "Gaussian";
        case CopulaType::STUDENT_T: return "Student-t";
    }
    return "Unknown";
}

std::string likelihood_mode_name(StatisticLikelihoodMode mode) {
    switch (mode) {
        case StatisticLikelihoodMode::PROFILED_NUISANCE:   return "profiled nuisance";
        case StatisticLikelihoodMode::CHI2_MC_COVARIANCE:  return "chi-square with MC covariance";
    }
    return "unknown";
}

std::string compact_double(double x) {
    std::ostringstream oss;
    oss << std::setprecision(8) << x;
    return oss.str();
}

std::string marginal_config_summary(const MarginalConfig& cfg) {
    return std::visit([](const auto& c) -> std::string {
        using T = std::decay_t<decltype(c)>;
        std::ostringstream oss;
        oss << std::setprecision(8);

        if constexpr (std::is_same_v<T, FlatMarginalCfg>) {
            oss << "a=" << c.a << ", b=" << c.b;
        } else if constexpr (std::is_same_v<T, GaussianMarginalCfg>) {
            oss << "mu=" << c.mu << ", sigma=" << c.sigma;
        } else if constexpr (std::is_same_v<T, SplitGaussianMarginalCfg>) {
            oss << "mu=" << c.mu << ", sigma_p=" << c.sigma_p
                << ", sigma_m=" << c.sigma_m;
        } else if constexpr (std::is_same_v<T, LikelihoodMarginalCfg>) {
            oss << "values=" << c.values.size()
                << ", weights=" << c.weights.size();
        } else {
            oss << "custom config";
        }
        return oss.str();
    }, cfg);
}

double safe_param_value(const std::shared_ptr<IStatParameterProxy>& pspp, const ParamId& pid) {
    try {
        return pspp->get_param(pid)->get_val().real();
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

double safe_param_sigma(const std::shared_ptr<IStatParameterProxy>& pspp, const ParamId& pid) {
    try {
        return std::abs(pspp->get_param(pid)->get_combined_std().real());
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

void print_nuisance_candidate_table(const std::string& title,
                                    const std::vector<ParamId>& ids,
                                    const std::shared_ptr<IStatParameterProxy>& pspp)
{
    std::cout << "[MC CONFIG] " << title << ": " << ids.size() << "\n";
    for (const auto& pid : ids) {
        std::cout << "[MC CONFIG]   candidate " << pid
                  << " | value=" << compact_double(safe_param_value(pspp, pid))
                  << " | sigma=" << compact_double(safe_param_sigma(pspp, pid))
                  << "\n";
    }
}

void print_nuisance_value_table(const std::string& title,
                                const std::map<ParamId, double>& values,
                                const std::shared_ptr<IStatParameterProxy>& pspp)
{
    std::cout << "[MC CONFIG] " << title << ": " << values.size() << "\n";
    for (const auto& [pid, value] : values) {
        std::cout << "[MC CONFIG]   retained " << pid
                  << " | value=" << compact_double(value)
                  << " | sigma=" << compact_double(safe_param_sigma(pspp, pid))
                  << "\n";
    }
}

void print_nuisance_marginal_line(const ParamId& pid,
                                  double value,
                                  double sigma,
                                  MarginalType mt,
                                  const MarginalConfig& cfg)
{
    std::cout << "[MC CONFIG]   retained " << pid
              << " | value=" << compact_double(value)
              << " | sigma=" << compact_double(sigma)
              << " | marginal=" << marginal_type_name(mt)
              << " (" << marginal_config_summary(cfg) << ")\n";
}

static std::vector<std::map<ParamId, double>> make_sensitivity_contexts(
    const std::map<ParamId, double>& eta0,
    const std::shared_ptr<IStatParameterProxy>& pspp,
    std::size_t n_contexts,
    double context_sigma,
    unsigned seed
) {
    std::vector<std::map<ParamId, double>> contexts;
    contexts.reserve(std::max<std::size_t>(1, n_contexts));

    contexts.push_back(eta0);

    std::mt19937 rng(seed);
    std::normal_distribution<double> normal(0.0, context_sigma);

    for (std::size_t c = 1; c < n_contexts; ++c) {
        auto ctx = eta0;

        for (auto& [pid, val] : ctx) {
            const double sigma =
                std::abs(pspp->get_param(pid)->get_combined_std().real());

            if (!std::isfinite(sigma) || sigma <= 0.0) {
                continue;
            }

            const double z = normal(rng);

            val += z * sigma;
        }

        contexts.push_back(std::move(ctx));
    }

    return contexts;
}

double fit_parameter_offset(const StatisticConfig& config, const ParamId& pid) {
    const auto it = config.fit_parameter_offsets.find(pid);
    return it == config.fit_parameter_offsets.end() ? 0.0 : it->second;
}

std::vector<double> fit_coordinates_from_model_values(
    const StatisticConfig& config,
    const std::vector<ParamId>& p_ids,
    std::vector<double> model_values)
{
    for (std::size_t i = 0; i < p_ids.size() && i < model_values.size(); ++i) {
        model_values[i] += fit_parameter_offset(config, p_ids[i]);
    }
    return model_values;
}

std::vector<double> model_values_from_fit_coordinates(
    const StatisticConfig& config,
    const std::vector<ParamId>& p_ids,
    const std::vector<double>& fit_values)
{
    std::vector<double> model_values = fit_values;
    for (std::size_t i = 0; i < p_ids.size() && i < model_values.size(); ++i) {
        model_values[i] -= fit_parameter_offset(config, p_ids[i]);
    }
    return model_values;
}

fit_app::ParameterDefinition make_fit_param_def(
    const ParamId& pid,
    double value,
    double sigma_hint,
    const StatisticConfig& config)
{
    fit_app::ParameterDefinition out;
    out.name = param_name(pid);
    out.value = value;
    out.step_hint = (std::isfinite(sigma_hint) && sigma_hint > 0.0)
        ? sigma_hint
        : std::max(1e-3, 0.01 * std::max(1.0, std::abs(value)));

    const auto configured_bounds = config.fit_parameter_bounds.find(pid);
    if (configured_bounds != config.fit_parameter_bounds.end()) {
        const auto [low, high] = configured_bounds->second;
        if (!std::isfinite(low) || !std::isfinite(high) || !(low < high)) {
            throw std::invalid_argument("Invalid explicit fit bounds for " + param_name(pid));
        }
        out.limits = configured_bounds->second;
        return out;
    }

    const std::string& nm = out.name;
    if (nm.find("FCONST") != std::string::npos) {
        out.limits = std::make_pair(0.05, 0.35);
    }

    return out;
}

template <class PredMapT>
std::vector<double> ordered_prediction_vector(
    const std::vector<ExperimentObs>& obs_ids,
    const PredMapT& pred_map)
{
    std::vector<double> out;
    out.reserve(obs_ids.size());

    for (const auto& bid : obs_ids) {
        const auto& vec = pred_map.at(bid.obs.s);

        auto it = std::find_if(vec.begin(), vec.end(), [&](const auto& ov) {
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


std::vector<BinnedObservableId> binned_ids_from_experiment_obs(
    const std::vector<ExperimentObs>& obs_ids
) {
    std::vector<BinnedObservableId> out;
    out.reserve(obs_ids.size());
    for (const auto& oid : obs_ids) {
        out.push_back(oid.obs);
    }
    return out;
}


RealMatrix symmetrize_covariance_matrix(RealMatrix cov) {
    if (cov.rows() != cov.cols()) {
        throw std::runtime_error("symmetrize_covariance_matrix: covariance must be square");
    }

    for (std::size_t i = 0; i < cov.rows(); ++i) {
        for (std::size_t j = i + 1; j < cov.cols(); ++j) {
            const double a = cov.at(i, j);
            const double b = cov.at(j, i);

            if (!std::isfinite(a) || !std::isfinite(b)) {
                throw std::runtime_error("symmetrize_covariance_matrix: non-finite covariance entry");
            }

            const double v = 0.5 * (a + b);
            cov.at(i, j) = v;
            cov.at(j, i) = v;
        }

        if (!std::isfinite(cov.at(i, i))) {
            throw std::runtime_error("symmetrize_covariance_matrix: non-finite covariance diagonal");
        }
    }

    return cov;
}

// RealMatrix inverse_covariance_with_ridge(
//     RealMatrix cov,
//     double ridge_rel,
//     double ridge_abs
// ) {
//     cov = symmetrize_covariance_matrix(std::move(cov));

//     double trace = 0.0;
//     for (std::size_t i = 0; i < cov.rows(); ++i) {
//         trace += std::max(0.0, cov.at(i, i));
//     }

//     const double scale = cov.rows() > 0
//         ? trace / static_cast<double>(cov.rows())
//         : 1.0;

//     const double ridge = std::max(
//         ridge_abs,
//         ridge_rel * std::max(scale, 1.0)
//     );

//     for (std::size_t i = 0; i < cov.rows(); ++i) {
//         cov.at(i, i) += ridge;
//     }

//     cov = symmetrize_covariance_matrix(std::move(cov));
//     return cov.inv();
// }

RealMatrix experimental_covariance_matrix(
    const std::vector<ExperimentObs>& obs_ids,
    const std::map<ExperimentObs, std::map<ExperimentObs, double>>& corr_obs,
    const std::map<ExperimentObs, double>& sigma_obs
) {
    const std::size_t n = obs_ids.size();
    RealMatrix cov(n, n);

    for (std::size_t i = 0; i < n; ++i) {
        const auto sig_i_it = sigma_obs.find(obs_ids[i]);
        if (sig_i_it == sigma_obs.end()) {
            std::ostringstream oss;
            oss << "experimental_covariance_matrix: missing sigma for observable "
                << obs_ids[i].str();
            throw std::runtime_error(oss.str());
        }

        const double sigma_i = std::abs(sig_i_it->second);
        if (!std::isfinite(sigma_i)) {
            throw std::runtime_error("experimental_covariance_matrix: non-finite sigma_i");
        }

        for (std::size_t j = 0; j < n; ++j) {
            const auto sig_j_it = sigma_obs.find(obs_ids[j]);
            if (sig_j_it == sigma_obs.end()) {
                std::ostringstream oss;
                oss << "experimental_covariance_matrix: missing sigma for observable "
                    << obs_ids[j].str();
                throw std::runtime_error(oss.str());
            }

            const double sigma_j = std::abs(sig_j_it->second);
            if (!std::isfinite(sigma_j)) {
                throw std::runtime_error("experimental_covariance_matrix: non-finite sigma_j");
            }

            double corr = (i == j) ? 1.0 : 0.0;
            auto row_it = corr_obs.find(obs_ids[i]);
            if (row_it != corr_obs.end()) {
                auto col_it = row_it->second.find(obs_ids[j]);
                if (col_it != row_it->second.end()) {
                    corr = col_it->second;
                }
            }

            if (!std::isfinite(corr)) {
                throw std::runtime_error("experimental_covariance_matrix: non-finite correlation");
            }

            cov.at(i, j) = corr * sigma_i * sigma_j;
        }
    }

    return symmetrize_covariance_matrix(std::move(cov));
}

RealMatrix inverse_covariance_with_ridge(
    RealMatrix cov,
    double ridge_rel,
    double ridge_abs
) {
    cov = symmetrize_covariance_matrix(std::move(cov));

    const std::size_t n = cov.rows();
    if (n != cov.cols()) {
        throw std::runtime_error("inverse_covariance_with_ridge: covariance must be square");
    }

    std::vector<double> sigma(n);

    for (std::size_t i = 0; i < n; ++i) {
        const double vii = cov.at(i, i);

        if (!std::isfinite(vii) || vii <= 0.0) {
            std::ostringstream oss;
            oss << "inverse_covariance_with_ridge: non-positive variance at i="
                << i << ", variance=" << vii;
            throw std::runtime_error(oss.str());
        }

        sigma[i] = std::sqrt(vii);
    }

    // Build dimensionless correlation matrix.
    RealMatrix corr(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            corr.at(i, j) = cov.at(i, j) / (sigma[i] * sigma[j]);
        }
    }

    corr = symmetrize_covariance_matrix(std::move(corr));

    // Ridge in correlation space, dimensionless.
    const double ridge = std::max(ridge_rel, ridge_abs);

    for (std::size_t i = 0; i < n; ++i) {
        corr.at(i, i) += ridge;
    }

    corr = symmetrize_covariance_matrix(std::move(corr));

    RealMatrix corr_inv = corr.inv();

    // Convert back: Cov^{-1} = D^{-1} Corr^{-1} D^{-1}
    RealMatrix cov_inv(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            cov_inv.at(i, j) = corr_inv.at(i, j) / (sigma[i] * sigma[j]);
        }
    }

    return symmetrize_covariance_matrix(std::move(cov_inv));
}

MLFitOptions make_mlfit_options_from_config(const StatisticConfig& config) {
    MLFitOptions fitopt;
    fitopt.run_hesse = config.advanced.MLE_run_hesse;
    fitopt.request_minos = config.advanced.MLE_request_minos;
    fitopt.verbose = config.advanced.MLE_verbose;
    fitopt.strategy = config.advanced.MLE_strategy;
    fitopt.max_fcn = static_cast<unsigned>(config.advanced.MLE_max_iter);
    fitopt.tolerance = config.advanced.MLE_tol;

    fitopt.allow_profile_hessian_fallback = config.advanced.MLE_allow_profile_hessian_fallback;
    fitopt.profile_hessian_step_scale = config.advanced.MLE_profile_hessian_step_scale;
    fitopt.profile_hessian_eig_floor_rel = config.advanced.MLE_profile_hessian_eig_floor_rel;

    fitopt.trace_first_evals = config.advanced.MLE_trace_first_evals;
    fitopt.trace_max_evals = config.advanced.MLE_trace_max_evals;

    return fitopt;
}

MCConfig make_mc_config_from_config(const StatisticConfig& config, bool for_chi2_covariance = false) {
    MCConfig cfg;
    cfg.draws = config.MC_draws;
    cfg.n_threads = config.MC_threads;
    cfg.force_decay_threads_to_one = config.advanced.MC_force_decay_threads_to_one;
    cfg.forced_decay_threads = config.advanced.MC_forced_decay_threads;
    cfg.skew_abs_threshold = config.skew_abs_threshold;
    cfg.covariance_ridge_rel = for_chi2_covariance
        ? config.advanced.chi2_covariance_ridge_rel
        : cfg.covariance_ridge_rel;
    cfg.covariance_ridge_abs = for_chi2_covariance
        ? config.advanced.chi2_covariance_ridge_abs
        : cfg.covariance_ridge_abs;
    cfg.print_progress = config.print_mc_progress;
    cfg.progress_probe_draws = config.mc_progress_probe_draws;
    cfg.progress_update_every = config.mc_progress_update_every;
    cfg.progress_monitor = config.progress_monitor;
    cfg.write_samples_csv = config.write_mc_samples_csv;
    cfg.samples_csv_path = config.mc_samples_csv_path;
    return cfg;
}


}

StatisticManager::StatisticManager(StatisticConfig config,
                                   std::shared_ptr<IModel> obs_int,
                                   std::shared_ptr<IStatCorrelationProxy> pscp,
                                   std::shared_ptr<IStatParameterProxy> pspp,
                                   std::shared_ptr<IStatSourcesProxy> sp,
                                   std::shared_ptr<IStatDependencyPruner> dp,
                                   std::shared_ptr<INuisanceReader> nuisance_reader,
                                   std::shared_ptr<IStatParamOptimizerProxy> spop)
    : obs_int(std::move(obs_int)),
      pscp(std::move(pscp)),
      pspp(std::move(pspp)),
      sp(std::move(sp)),
      dp(std::move(dp)),
      nuisance_reader_(std::move(nuisance_reader)),
      spop(std::move(spop)),
      config(std::move(config))
{
    if (!nuisance_reader_) {
        throw std::invalid_argument("StatisticManager: nuisance_reader is null");
    }

    this->obs_int->compute_observables();
    reload_nuisance_specs();
    invalidate_fit_state();
}

std::vector<std::unique_ptr<IMarginalDistribution>> StatisticManager::build_nuisance_marginal_distributions() {
    unsigned int seed = config.MC_seed;
    std::vector<std::unique_ptr<IMarginalDistribution>> marginals;

    marginals.reserve(cache.eta_specs_real.size());

    const bool print_config = config.print_mc_config || config.print_debug;
    if (print_config) {
        std::cout << "[MC CONFIG] Retained nuisance marginals used by MC: "
                  << cache.eta_specs_real.size() << "\n";
    }

    for (const auto& [pid, value] : cache.eta_specs_real) {
        const MarginalType mt = resolve_nuisance_marginal_type(pid);
        MarginalConfig cfg = make_nuisance_marginal_config(pid, mt);

        if (print_config) {
            print_nuisance_marginal_line(
                pid,
                value,
                safe_param_sigma(pspp, pid),
                mt,
                cfg
            );
        }

        marginals.emplace_back(MarginalFactory::create(mt, cfg, seed));
    }

    return marginals;
}

std::unique_ptr<JointDistribution> StatisticManager::build_nuisance_distribution() {
    unsigned int seed = config.MC_seed;

    if (config.print_mc_config || config.print_debug) {
        std::cout << "[MC CONFIG] Nuisance distribution: "
                  << cache.eta_specs_real.size() << " retained nuisance(s), copula="
                  << copula_type_name(config.advanced.nuisance_copula_type)
                  << "\n";
    }

    std::unique_ptr<ICopula> copula;
    if (config.advanced.nuisance_copula_type == CopulaType::GAUSSIAN) {
        GaussianCopulaConfig copula_cfg;
        copula_cfg.R = RealMatrix(unzip(cache.SigmaEta).vals);
        copula = CopulaFactory::create(config.advanced.nuisance_copula_type, copula_cfg, seed);
    } else if (config.advanced.nuisance_copula_type == CopulaType::STUDENT_T) {
        StudentTCopulaConfig copula_cfg;
        copula_cfg.R = RealMatrix(unzip(cache.SigmaEta).vals);
        copula_cfg.nu = cache.eta_specs_real.size();
        copula = CopulaFactory::create(config.advanced.nuisance_copula_type, copula_cfg, seed);
    }

    return std::make_unique<JointDistribution>(
        std::move(build_nuisance_marginal_distributions()), 
        std::move(copula)
    );
}

std::unique_ptr<JointDistribution> StatisticManager::build_exp_data_distribution() {
    // Keep the experimental-data stream reproducible but decorrelated from the
    // nuisance stream initialized with MC_seed.
    unsigned int seed = config.MC_seed ^ 0x9E3779B9u;
    std::map<ExperimentObs, MarginalType> exp_data_marginals;

    for (auto& [oid, v] : cache.exp_obs) {
        exp_data_marginals.emplace(oid, MarginalType::GAUSSIAN);
    }

    for (auto& [oid, mt] : config.advanced.override_exp_data_marginals) {
        if (!exp_data_marginals.contains(oid)) {
            LOG_WARN("Observable", ObservableMapper::str(oid.obs.s), "is not a taken into account in this fit.");
            continue;
        }
        exp_data_marginals.at(oid) = mt;
    }

    auto unzipped = unzip(exp_data_marginals);
    std::vector<ExperimentObs> obs_ids = unzipped.ids;
    std::vector<std::unique_ptr<IMarginalDistribution>> marginals;

    for (auto& [oid, mt] : exp_data_marginals) {
        MarginalConfig cfg = MarginalConfigFactory().create(oid, mt);
        auto m_ptr = MarginalFactory::create(mt, cfg, seed);
        marginals.emplace_back(std::move(m_ptr));
    }

    std::unique_ptr<ICopula> copula;
    if (config.advanced.exp_data_copula_type == CopulaType::GAUSSIAN) {
        GaussianCopulaConfig copula_cfg;
        copula_cfg.R = RealMatrix(unzip(cache.SigmaObs).vals);
        copula = CopulaFactory::create(config.advanced.exp_data_copula_type, copula_cfg, seed);
    } else if (config.advanced.exp_data_copula_type == CopulaType::STUDENT_T) {
        StudentTCopulaConfig copula_cfg;
        copula_cfg.R = RealMatrix(unzip(cache.SigmaObs).vals);
        copula_cfg.nu = obs_int->n_observables() - 1;
        copula = CopulaFactory::create(config.advanced.exp_data_copula_type, copula_cfg, seed);
    }

    return std::make_unique<JointDistribution>(std::move(marginals), std::move(copula));
}

std::map<BinnedObservableId, GaussianSummary> StatisticManager::compute_uncertainties() {
    auto sums = this->compute_uncertainties_and_sampling();

    if (config.print_debug) {
        std::cout << "[DEBUG] Merged nuisance specification registry: "
                  << this->merged_nuisance_specs_.size() << " entry/entries\n";
        for (const auto& elem : this->merged_nuisance_specs_) {
            std::cout << "[DEBUG]   " << elem.first << "\n";
            std::cout << elem.second << "\n";
        }
    }

    std::map<BinnedObservableId, GaussianSummary> out;
    if (config.print_mc_config || config.print_debug) {
        std::cout << "[MC CONFIG] Observable MC summaries: "
                  << sums.summary.size() << "\n";
    }
    for (const auto& gs : sums.summary) {
        out[gs.id] = gs;
        if (config.print_mc_config || config.print_debug) {
            std::cout << "[MC CONFIG]   " << gs << "\n";
        }
    }
    return out;
}

MCResult StatisticManager::compute_uncertainties_and_sampling() {
    if (this->config.progress_monitor) {
        this->config.progress_monitor->reset("preparing", "Preparing Monte-Carlo uncertainty propagation");
    }
    update_cache();
    auto rvg = build_nuisance_distribution();
    std::vector<ParamId> nuisance_ids = unzip(cache.eta_specs_real).ids;
    RvgNuisanceSampler sampler(nuisance_ids, std::move(rvg));
    MonteCarloEngine mc(this->obs_int, sampler, make_mc_config_from_config(this->config));

    auto sums = mc.summarize(this->cache.p_specs);

    return sums;
}

FitResultWithMaps StatisticManager::compute_MLE(const std::vector<ParamId>& p_specs) {
    if (this->config.progress_monitor) {
        this->config.progress_monitor->reset("preparing", "Preparing chi-square fit");
    }
    update_cache(p_specs);

    if (cache.p_specs.empty()) {
        throw std::invalid_argument("compute_MLE called with an empty fit parameter list.");
    }

    auto unzipped_fit_params = unzip(cache.p_specs);
    auto unzipped_nuisances  = unzip(cache.eta_specs_real);
    auto unzipped_exp_obs    = unzip(cache.exp_obs);

    if (config.print_fit_summary || config.print_debug) {
        auto eta_ids_dbg = unzip(cache.eta_specs_real).ids;
        RealMatrix Reta(unzip(cache.SigmaEta).vals);
        dump_matrix_sanity(Reta, eta_ids_dbg, "SigmaEta");
    }

    const std::vector<ParamId> p_ids = unzipped_fit_params.ids;
    const std::vector<double> p0 = fit_coordinates_from_model_values(
        config, p_ids, unzipped_fit_params.vals
    );

    const std::vector<ParamId> eta_ids = unzipped_nuisances.ids;
    const std::vector<double> eta0 = unzipped_nuisances.vals;
    
    RealMatrix SigmaEta(unzip(cache.SigmaEta).vals);

    const std::vector<ExperimentObs> obs_ids = unzipped_exp_obs.ids;
    const std::vector<double> exp_obs_vals = unzipped_exp_obs.vals;

    if (config.advanced.likelihood_mode == StatisticLikelihoodMode::CHI2_MC_COVARIANCE) {
        if (config.print_fit_summary || config.print_debug) {
            std::cout << "[FIT] Likelihood backend: "
                      << likelihood_mode_name(config.advanced.likelihood_mode) << ".\n";
        }

        const bool show_chi2_progress =
            config.print_chi2_pipeline_progress ||
            config.print_mc_progress ||
            config.print_debug;
        StatisticStageProgressReporter chi2_progress(
            show_chi2_progress,
            "CHI2 workflow",
            7,
            std::cout,
            this->config.progress_monitor,
            "chi2_pipeline"
        );
        chi2_progress.start(
            "Starting chi-square workflow: MC covariance, experimental covariance, likelihood construction, then maximum-likelihood fit. The MC bar estimates only the MC part; later steps are reported separately."
        );

        auto rvg = build_nuisance_distribution();
        std::vector<ParamId> nuisance_ids = unzip(cache.eta_specs_real).ids;
        RvgNuisanceSampler sampler(nuisance_ids, std::move(rvg));

        MonteCarloEngine mc(
            this->obs_int,
            sampler,
            make_mc_config_from_config(this->config, true)
        );

        MCRealization mc_real = mc.sample_predictions(this->cache.p_specs);
        chi2_progress.step(
            1,
            "Monte-Carlo sampling completed; building the MC covariance matrix."
        );

        const std::vector<BinnedObservableId> cov_ids =
            binned_ids_from_experiment_obs(obs_ids);

        MCObservableCovariance cov = covariance_from_obs_samples(
            mc_real.sampled_obss,
            cov_ids,
            this->config.advanced.chi2_covariance_ridge_rel,
            this->config.advanced.chi2_covariance_ridge_abs
        );
        chi2_progress.step(
            2,
            "MC covariance ready; collecting experimental uncertainties."
        );

        std::map<ExperimentObs, double> exp_obs_sigmas;
        for (const auto& model_obs_id : this->obs_int->get_obs_ids()) {
            auto exp_params = pspp->get_obs_param(model_obs_id);
            for (const auto& [exp_obs, param] : exp_params) {
                if (!accepts_experiment_observable(exp_obs)) {
                    continue;
                }
                exp_obs_sigmas[exp_obs] =
                    std::abs(param->get_combined_std().real());
            }
        }

        chi2_progress.step(
            3,
            "Experimental uncertainties ready; assembling the total covariance matrix."
        );
        const RealMatrix covariance_exp = experimental_covariance_matrix(
            obs_ids,
            this->cache.SigmaObs,
            exp_obs_sigmas
        );

        RealMatrix covariance_total =
            symmetrize_covariance_matrix(cov.covariance + covariance_exp);

        // RealMatrix covariance_total =
        //     symmetrize_covariance_matrix(covariance_exp);

        //     for (std::size_t i = 0; i < covariance_total.rows(); ++i) {
        //         for (std::size_t j = 0; j < covariance_total.cols(); ++j) {
        //             if (i != j) {
        //                 covariance_total.at(i, j) = 0.0;
        //             }
        //         }
        //     }

        chi2_progress.step(
            4,
            "Total covariance assembled; inverting the regularized covariance matrix."
        );
        RealMatrix covariance_total_inv = inverse_covariance_with_ridge(
            covariance_total,
            this->config.advanced.chi2_covariance_ridge_rel,
            this->config.advanced.chi2_covariance_ridge_abs
        );
        chi2_progress.step(
            5,
            "Covariance inverse ready; constructing the chi-square likelihood."
        );

        if (config.print_fit_summary || config.print_debug) {
            std::cout << "[FIT] Covariance model: total = MC theory covariance + experimental covariance.\n";
        }

        auto ctx = std::make_shared<LikelihoodContext>();
        ctx->exp_obs_values = exp_obs_vals;
        ctx->nuis_defs.clear();
        ctx->fp_defs.reserve(p_ids.size());

        for (std::size_t i = 0; i < p_ids.size(); ++i) {
            double sigma = std::abs(pspp->get_param(p_ids[i])->get_combined_std().real());
            ctx->fp_defs.emplace_back(make_fit_param_def(p_ids[i], p0[i], sigma, config));
        }

        // auto model_fn = [this, obs_ids, p_ids](
        //     const std::vector<double>& p_vec,
        //     const std::vector<double>& eta_vec) -> std::vector<double>
        // {
        //     if (!eta_vec.empty()) {
        //         throw std::runtime_error("CHI2_MC_COVARIANCE model_fn expects empty eta vector");
        //     }

        //     auto pred_map = this->obs_int->predict_optimized(
        //         zip(p_ids, p_vec),
        //         std::map<ParamId, double>{}
        //     );

        //     return ordered_prediction_vector(obs_ids, pred_map);
        // };

        const std::map<ParamId, double> eta0_map = zip(eta_ids, eta0);

        auto model_fn = [this, obs_ids, p_ids, eta0_map](
            const std::vector<double>& p_vec,
            const std::vector<double>& eta_vec) -> std::vector<double>
        {
            if (!eta_vec.empty()) {
                throw std::runtime_error("CHI2_MC_COVARIANCE model_fn expects empty eta vector");
            }

            // Important : en mode chi2 zéro nuisance, on ne veut pas laisser
            // les nuisances dans l'état du dernier tirage MC.
            // On force explicitement les nuisances au point central eta0.
            const auto model_p = model_values_from_fit_coordinates(this->config, p_ids, p_vec);
            auto pred_map = this->obs_int->predict_optimized(
                zip(p_ids, model_p),
                eta0_map
            );

            return ordered_prediction_vector(obs_ids, pred_map);
        };
        last_ctx_ = ctx;
        last_like_ = std::make_shared<ChiSquaredLikelihood>(
            model_fn,
            ctx,
            p_ids.size(),
            covariance_total_inv
        );

        chi2_progress.step(
            6,
            "Likelihood ready; running the maximum-likelihood fit. This backend-dependent step has no reliable ETA."
        );

        MLFitOptions fitopt = make_mlfit_options_from_config(config);
        last_fitter_ = std::make_shared<MLFitter>(last_like_, fitopt);
        last_fit_raw_ = last_fitter_->maximum_likelihood_fit(p0);
        chi2_progress.finish("Chi-square workflow completed.");

        last_fit_param_ids_ = p_ids;
        last_nuisance_ids_.clear();
        last_fit_param_index_.clear();
        for (std::size_t i = 0; i < p_ids.size(); ++i) {
            last_fit_param_index_[p_ids[i]] = i;
        }

        FitResultWithMaps out;
        out.fit_ok = !last_fit_raw_.p_hat.empty();
        out.ell_hat = last_fit_raw_.ell_hat;
        out.p_hat = zip(p_ids, last_fit_raw_.p_hat);
        out.eta_hat.clear();
        out.p_hat_std = zip(p_ids, last_fit_raw_.p_hat_std);
        out.p_correlations = zip(p_ids, last_fit_raw_.p_hat_correlations);

        cache.mle_result = out;
        return out;
    }

    auto ctx = std::make_shared<LikelihoodContext>();
    ctx->nuisance_dist = build_nuisance_distribution();
    ctx->exp_obs_dist = build_exp_data_distribution();
    ctx->exp_obs_values = exp_obs_vals;

    ctx->fp_defs.reserve(p_ids.size());
    for (std::size_t i = 0; i < p_ids.size(); ++i) {
        double sigma = std::abs(pspp->get_param(p_ids[i])->get_combined_std().real());
        ctx->fp_defs.emplace_back(make_fit_param_def(p_ids[i], p0[i], sigma, config));
    }

    ctx->nuis_defs.reserve(eta_ids.size());
    for (std::size_t i = 0; i < eta_ids.size(); ++i) {
        double sigma = std::abs(pspp->get_param(eta_ids[i])->get_combined_std().real());
        ctx->nuis_defs.emplace_back(
            make_nuisance_parameter_definition(eta_ids[i], eta0[i], sigma)
        );
    }

    auto model_fn = [this, obs_ids, p_ids, eta_ids](
        const std::vector<double>& p_vec,
        const std::vector<double>& eta_vec) -> std::vector<double>
    {
        const auto model_p = model_values_from_fit_coordinates(this->config, p_ids, p_vec);
        auto pred_map = this->obs_int->predict_optimized(
            zip(p_ids, model_p),
            zip(eta_ids, eta_vec)
        );

        return ordered_prediction_vector(obs_ids, pred_map);
    };

    last_ctx_ = ctx;
    last_like_ = std::make_shared<BaseLikelihood>(model_fn, ctx, p_ids.size());

    MLFitOptions fitopt;
    fitopt.run_hesse = config.advanced.MLE_run_hesse;
    fitopt.request_minos = config.advanced.MLE_request_minos;
    fitopt.verbose = config.advanced.MLE_verbose;
    fitopt.strategy = config.advanced.MLE_strategy;
    fitopt.max_fcn = static_cast<unsigned>(config.advanced.MLE_max_iter);
    fitopt.tolerance = config.advanced.MLE_tol;

    fitopt.allow_profile_hessian_fallback = config.advanced.MLE_allow_profile_hessian_fallback;
    fitopt.profile_hessian_step_scale = config.advanced.MLE_profile_hessian_step_scale;
    fitopt.profile_hessian_eig_floor_rel = config.advanced.MLE_profile_hessian_eig_floor_rel;

    fitopt.trace_first_evals = config.advanced.MLE_trace_first_evals;
    fitopt.trace_max_evals = config.advanced.MLE_trace_max_evals;

    last_fitter_ = std::make_shared<MLFitter>(ctx, model_fn, fitopt);
    last_fit_raw_ = last_fitter_->maximum_likelihood_fit(p0);

    last_fit_param_ids_ = p_ids;
    last_nuisance_ids_ = eta_ids;
    last_fit_param_index_.clear();
    for (std::size_t i = 0; i < p_ids.size(); ++i) {
        last_fit_param_index_[p_ids[i]] = i;
    }

    FitResultWithMaps out;
    out.fit_ok = !last_fit_raw_.p_hat.empty();
    out.ell_hat = last_fit_raw_.ell_hat;
    out.p_hat = zip(p_ids, last_fit_raw_.p_hat);
    out.eta_hat = zip(eta_ids, last_fit_raw_.eta_hat);
    out.p_hat_std = zip(p_ids, last_fit_raw_.p_hat_std);
    out.p_correlations = zip(p_ids, last_fit_raw_.p_hat_correlations);

    cache.mle_result = out;
    return out;
}

Contour StatisticManager::confidence_contour(ParamId p1, ParamId p2, double z, std::array<double, 4> bounds, ContourOptions options) {
    if (!cache.mle_result.fit_ok || !last_fitter_) {
        throw std::runtime_error("Please run compute_MLE before requesting a confidence contour.");
    }

    if (!last_fit_param_index_.contains(p1) || !last_fit_param_index_.contains(p2)) {
        throw std::invalid_argument("Contour requested for parameters that are not in the last fitted parameter set.");
    }

    if (p1 == p2) {
        throw std::invalid_argument("Contour requires two distinct parameters.");
    }

    const std::size_t x_id = last_fit_param_index_.at(p1);
    const std::size_t y_id = last_fit_param_index_.at(p2);

    Contour cl = last_fitter_->contour(
        x_id,
        y_id,
        z,
        bounds,
        options
    );

    return cl;
}

void StatisticManager::print_cache()
{
    if (!(config.print_cache_summary || config.print_debug)) {
        return;
    }

    StatisticCachePrinter::print(cache, std::cout);
}

void StatisticManager::update_cache(const std::vector<ParamId>& p_specs) {
    if (selected_experiments_.has_value()) {
        for (const auto& elem : *selected_experiments_) {
            LOG_VERBOSE("USING SELECTED EXPERIMENT: ", elem);
        }
    } else {
        LOG_VERBOSE("USING ALL EXPERIMENTS");
    }

    if (selected_experiment_observables_.has_value()) {
        LOG_VERBOSE("USING EXPLICIT EXPERIMENT-OBSERVABLE SELECTION: ",
                    selected_experiment_observables_->size(), " entries");
    }

    for (const auto& [tp, block] : last_detached_fit_blocks_) {
        dp->reattach_block(tp, block);
    }
    for (const auto& pid : last_detached_fit_params_) {
        if (pid.type.has_value()) {
            dp->reattach_parameter(pid.type.value(), pid.block, pid.code);
        }
    }
    last_detached_fit_blocks_.clear();
    last_detached_fit_params_.clear();

    cache.p_specs = this->get_p_specs(p_specs);

    std::unordered_set<std::string> seen_blocks;

    for (const auto& [pid, _] : cache.p_specs) {
        if (!pid.type.has_value()) {
            continue;
        }

        const auto tp = pid.type.value();
        const std::string block_key =
            std::to_string(static_cast<int>(tp)) + "::" + pid.block.to_string();

        if (!seen_blocks.contains(block_key)) {
            dp->detach_block(tp, pid.block);
            last_detached_fit_blocks_.push_back({tp, pid.block});
            seen_blocks.insert(block_key);
        }

        dp->detach_parameter(tp, pid.block, pid.code);
        last_detached_fit_params_.push_back(pid);
    }

    cache.eta_specs_real = this->get_all_obss_deps();
    for (const auto& [pid, _] : cache.p_specs)
        cache.eta_specs_real.erase(pid);
    
    for (auto it = cache.eta_specs_real.begin(); it != cache.eta_specs_real.end(); ) {
        const ParamId& pid = it->first;

        if (pid.block.to_string().find("__BSM") != std::string::npos) {
            LOG_INFO("Dropping BSM nuisance from cache", pid);
            it = cache.eta_specs_real.erase(it);
        } else {
            ++it;
        }
    }

    cache.SigmaEta = this->get_all_correlations();
    cache.exp_obs = this->get_obs_exp();
    cache.SigmaObs = this->get_all_obs_correlations();

    std::ofstream fs;
    fs.open("covariance.csv");
    for (auto& [pid1, row] : cache.SigmaEta) {
        double sigma_1 = std::abs(pspp->get_param(pid1)->get_combined_std().real());
        for (auto& [pid2, corr] : row) {
            double sigma_2 = std::abs(pspp->get_param(pid2)->get_combined_std().real());
            if (pid2 == (*(--row.end())).first)
                fs << corr * sigma_1 * sigma_2;    
            else
                fs << corr * sigma_1 * sigma_2 << ',';    
        }
        fs << '\n';
    }

    fs.close();
}

void StatisticManager::reload_nuisance_specs() {
    default_nuisance_specs_ = nuisance_reader_->load_default();

    if (current_user_nuisance_file_.has_value()) {
        user_nuisance_specs_ = nuisance_reader_->load_user(*current_user_nuisance_file_);
    } else {
        user_nuisance_specs_ = nuisance_reader_->load_user();
    }

    rebuild_merged_nuisance_specs();
    invalidate_fit_state();
}

void StatisticManager::set_nuisance_user_file(const fs::path& user_yaml_path) {
    current_user_nuisance_file_ = user_yaml_path;
    reload_nuisance_specs();
}

void StatisticManager::clear_nuisance_user_file() {
    current_user_nuisance_file_.reset();
    reload_nuisance_specs();
}

void StatisticManager::rebuild_merged_nuisance_specs() {
    merged_nuisance_specs_ = default_nuisance_specs_;
    for (const auto& [pid, spec] : user_nuisance_specs_) {
        merged_nuisance_specs_[pid] = spec;
    }

}

void StatisticManager::invalidate_fit_state() {
    cache.mle_result = FitResultWithMaps{};

    last_ctx_.reset();
    last_like_.reset();
    last_fitter_.reset();

    last_fit_param_ids_.clear();
    last_nuisance_ids_.clear();
    last_fit_param_index_.clear();
}

const NuisanceSpec* StatisticManager::find_nuisance_spec(const ParamId& pid) const {
    if (const auto it = merged_nuisance_specs_.find(pid);
        it != merged_nuisance_specs_.end()) {
        return &it->second;
    }

    // JSON/YAML nuisance entries identify parameters by block/code only.
    // Runtime dependency ids are usually typed, so retry with the canonical
    // untyped key used by NuisanceReader.
    const ParamId untyped_pid{pid.block, pid.code};
    if (const auto it = merged_nuisance_specs_.find(untyped_pid);
        it != merged_nuisance_specs_.end()) {
        return &it->second;
    }

    return nullptr;
}

MarginalType StatisticManager::resolve_nuisance_marginal_type(const ParamId& pid) const {
    MarginalType mt = MarginalType::GAUSSIAN;

    if (const auto* spec = find_nuisance_spec(pid)) {
        mt = spec->marginal;
    }

    if (config.advanced.override_nuisance_marginals.contains(pid)) {
        mt = config.advanced.override_nuisance_marginals.at(pid);
    }

    return mt;
}

fit_app::ParameterDefinition StatisticManager::make_nuisance_parameter_definition(const ParamId& pid,
                                                                                  double value,
                                                                                  double sigma_hint) const
{
    fit_app::ParameterDefinition out;
    out.name = param_name(pid);
    out.value = value;

    const double s = (std::isfinite(sigma_hint) && sigma_hint > 0.0)
        ? std::abs(sigma_hint)
        : std::max(1e-3, 0.01 * std::abs(value));

    out.step_hint = s;

    if (const auto* spec = find_nuisance_spec(pid)) {
        out.limits = spec->bounds;
        return out;
    }

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

MarginalConfig StatisticManager::make_nuisance_marginal_config(const ParamId& pid,
                                                               MarginalType mt) const
{
    if (const auto* spec = find_nuisance_spec(pid)) {
        return MarginalConfigFactory().create(pid, mt, *spec);
    }

    return MarginalConfigFactory().create(pid, mt);
}

std::map<ParamId, double> StatisticManager::get_all_obss_deps() {

    
    std::unordered_set<ParamId> eta_infos;

    for (const auto& obsId : obs_int->get_obs_ids()) {
        for (auto paramId : obs_int->get_obs_deps(obsId.s)) {
            eta_infos.insert(paramId);
        }
    }

    std::unordered_set<ParamId> eta_infos_leaf = this->sp->get_all_leaf_sources(eta_infos);

    if (config.print_mc_config || config.print_debug) {
        std::vector<ParamId> candidate_ids(eta_infos_leaf.begin(), eta_infos_leaf.end());
        std::sort(candidate_ids.begin(), candidate_ids.end(), [](const ParamId& a, const ParamId& b) {
            return param_name(a) < param_name(b);
        });
        print_nuisance_candidate_table(
            "Potential nuisance candidates before pruning",
            candidate_ids,
            pspp
        );
    }

    std::map<ParamId, double> eta_specs_real_leaf;
    std::map<ParamId, double> delta_rel;

    // 1) Préfiltre existant : relevance basée sur l'incertitude relative
    for (const auto& paramId : eta_infos_leaf) {
        // On ne traite jamais les paramètres de fit comme nuisances
        if (cache.p_specs.contains(paramId)) {
            continue;
        }

        const double u   = std::abs(pspp->get_param(paramId)->get_combined_std().real());
        const double val = pspp->get_param(paramId)->get_val().real();

        if (!std::isfinite(u) || fpeq(u, 0.0)) {
            delta_rel[paramId] = 0.0;
            continue;
        }

        if (fpeq(val, 0.0)) {
            // Cas utile : valeur centrale nulle mais sigma non nulle.
            // On le garde au moins jusqu'au filtre de sensibilité.
            eta_specs_real_leaf[paramId] = val;
            delta_rel[paramId] = std::numeric_limits<double>::infinity();
            continue;
        }

        delta_rel[paramId] = std::abs(u / val);
    }

    double delta_rel_max = 0.0;
    for (const auto& [pid, d] : delta_rel) {
        if (std::isfinite(d)) {
            delta_rel_max = std::max(delta_rel_max, d);
        }
    }
    if (!(delta_rel_max > 0.0)) {
        delta_rel_max = 1.0;
    }

    for (const auto& [pid, d] : delta_rel) {
        const double rel_to_max = std::isfinite(d) ? (d / delta_rel_max) : 1.0;
        if (rel_to_max > config.advanced.nuisance_relevance_cutoff || !std::isfinite(d)) {
            eta_specs_real_leaf[pid] = pspp->get_param(pid)->get_val();
        }
    }

    if (config.print_mc_config || config.print_debug) {
        print_nuisance_value_table(
            "Nuisances after relative-uncertainty preselection",
            eta_specs_real_leaf,
            pspp
        );
    }

    if (config.advanced.nuisance_sensitivity_pruning &&
        config.advanced.nuisance_sensitivity_contexts >= 0 &&
        !eta_specs_real_leaf.empty() &&
        !cache.p_specs.empty())
    {
        const auto exp_obs_map = get_obs_exp();
        const auto unzipped_exp_obs = unzip(exp_obs_map);
        const std::vector<ExperimentObs> obs_ids = unzipped_exp_obs.ids;
        const std::vector<double> exp_obs_vals   = unzipped_exp_obs.vals;

        const std::map<ParamId, double> eta_central_for_restore = eta_specs_real_leaf;
        const std::map<ParamId, double> p_central_for_restore = cache.p_specs;

        struct RestoreModelStateGuard {
            StatisticManager* self = nullptr;
            std::map<ParamId, double> p;
            std::map<ParamId, double> eta;
            bool active = true;

            ~RestoreModelStateGuard() {
                if (!active || self == nullptr) return;

                try {
                    self->obs_int->predict_optimized(p, eta);
                } catch (const std::exception& e) {
                    std::cout << "[FIT] WARNING: failed to restore central model state: "
                            << e.what() << std::endl;
                } catch (...) {
                    std::cout << "[FIT] WARNING: failed to restore central model state: unknown exception"
                            << std::endl;
                }
            }
        };

        RestoreModelStateGuard restore_guard{
            this,
            p_central_for_restore,
            eta_central_for_restore,
            true
        };

        const auto contexts = make_sensitivity_contexts(
            eta_specs_real_leaf,
            pspp,
            config.advanced.nuisance_sensitivity_contexts,
            config.advanced.nuisance_sensitivity_context_sigma,
            config.advanced.nuisance_sensitivity_seed
        );

        std::map<ParamId, double> screened_eta_specs;
        const bool sensitivity_verbose = config.print_fit_summary || config.print_debug;

        if (sensitivity_verbose) {
            std::cout << "[FIT] Model-sensitivity pruning on "
                    << eta_specs_real_leaf.size()
                    << " nuisance candidates using "
                    << contexts.size()
                    << " contexts" << std::endl;
        }

        for (const auto& [pid, nominal] : eta_specs_real_leaf) {
            const double sigma =
                std::abs(pspp->get_param(pid)->get_combined_std().real());

            if (!std::isfinite(sigma) || fpeq(sigma, 0.0)) {
                screened_eta_specs[pid] = nominal;
                if (sensitivity_verbose) {
                    std::cout << "[FIT] sensitivity " << pid
                            << " : sigma not finite or zero -> keep" << std::endl;
                }
                continue;
            }

            const fit_app::ParameterDefinition def =
                make_nuisance_parameter_definition(pid, nominal, sigma);

            const double step = config.advanced.nuisance_sensitivity_probe_sigmas * sigma;

            if (!std::isfinite(step) || fpeq(step, 0.0)) {
                screened_eta_specs[pid] = nominal;
                if (sensitivity_verbose) {
                    std::cout << "[FIT] sensitivity " << pid
                            << " : probe step not finite or zero -> keep" << std::endl;
                }
                continue;
            }

            double best_abs_shift = 0.0;
            double best_rel_shift = 0.0;
            std::size_t best_context = 0;
            bool evaluation_failed = false;

            for (std::size_t c = 0; c < contexts.size(); ++c) {
                auto base_map = contexts[c];

                double center = base_map.at(pid);
                double eta_plus = center + step;
                double eta_minus = center - step;

                if (def.limits.has_value()) {
                    const auto [lo, hi] = def.limits.value();
                    eta_plus  = std::clamp(eta_plus,  lo, hi);
                    eta_minus = std::clamp(eta_minus, lo, hi);
                }

                if (!std::isfinite(eta_plus) ||
                    !std::isfinite(eta_minus) ||
                    std::abs(eta_plus - eta_minus) < 1e-14) {
                    continue;
                }

                try {
                    const auto baseline_pred_map =
                        this->obs_int->predict_optimized(cache.p_specs, base_map);

                    const std::vector<double> baseline_pred =
                        ordered_prediction_vector(obs_ids, baseline_pred_map);

                    auto eta_plus_map = base_map;
                    auto eta_minus_map = base_map;

                    eta_plus_map[pid] = eta_plus;
                    eta_minus_map[pid] = eta_minus;

                    const auto pred_plus_map =
                        this->obs_int->predict_optimized(cache.p_specs, eta_plus_map);

                    const auto pred_minus_map =
                        this->obs_int->predict_optimized(cache.p_specs, eta_minus_map);

                    const std::vector<double> pred_plus =
                        ordered_prediction_vector(obs_ids, pred_plus_map);

                    const std::vector<double> pred_minus =
                        ordered_prediction_vector(obs_ids, pred_minus_map);

                    double max_abs_shift = 0.0;
                    double max_rel_shift = 0.0;

                    for (std::size_t i = 0; i < baseline_pred.size(); ++i) {
                        const double one_sigma_shift =
                            0.5 * std::abs(pred_plus[i] - pred_minus[i]);

                        const double scale = std::max({
                            std::abs(baseline_pred[i]),
                            std::abs(exp_obs_vals[i]),
                            config.advanced.nuisance_sensitivity_scale_floor
                        });

                        max_abs_shift = std::max(max_abs_shift, one_sigma_shift);
                        max_rel_shift = std::max(max_rel_shift, one_sigma_shift / scale);
                    }

                    if (max_rel_shift > best_rel_shift ||
                        max_abs_shift > best_abs_shift) {
                        best_abs_shift = std::max(best_abs_shift, max_abs_shift);
                        best_rel_shift = std::max(best_rel_shift, max_rel_shift);
                        best_context = c;
                    }

                } catch (const std::exception& e) {
                    evaluation_failed = true;

                    if (sensitivity_verbose) {
                        std::cout << "[FIT] sensitivity " << pid
                                << " : context " << c
                                << " failed with exception: "
                                << e.what() << std::endl;
                    }

                    break;
                } catch (...) {
                    evaluation_failed = true;

                    if (sensitivity_verbose) {
                        std::cout << "[FIT] sensitivity " << pid
                                << " : context " << c
                                << " failed with unknown exception"
                                << std::endl;
                    }

                    break;
                }
            }

            if (evaluation_failed && config.advanced.nuisance_sensitivity_keep_on_failure) {
                screened_eta_specs[pid] = nominal;

                if (sensitivity_verbose) {
                    std::cout << "[FIT] sensitivity " << pid
                            << " : evaluation failed -> keep" << std::endl;
                }

                continue;
            }

            const bool keep =
                (best_abs_shift >= config.advanced.nuisance_sensitivity_abs_cutoff) ||
                (best_rel_shift >= config.advanced.nuisance_sensitivity_rel_cutoff);

            if (sensitivity_verbose) {
                std::cout << "[FIT] sensitivity " << pid
                        << " : max_abs_shift=" << best_abs_shift
                        << ", max_rel_shift=" << best_rel_shift
                        << ", best_context=" << best_context
                        << " -> " << (keep ? "keep" : "drop")
                        << std::endl;
            }

            if (keep) {
                screened_eta_specs[pid] = nominal;
            }
        }

        eta_specs_real_leaf = std::move(screened_eta_specs);
    }

    if (config.print_mc_config || config.print_debug) {
        print_nuisance_value_table(
            "Final retained nuisances passed to MC/fit",
            eta_specs_real_leaf,
            pspp
        );
    }

    if (config.print_fit_summary || config.print_debug) {
        std::cout << "[FIT] Significant nuisances retained: "
                  << eta_specs_real_leaf.size() << "\n";
        for (const auto& [pid, val] : eta_specs_real_leaf) {
            std::cout << "[FIT]   " << pid << " = " << compact_double(val) << "\n";
        }
    }

    return eta_specs_real_leaf;
}

std::map<ParamId, double> StatisticManager::get_p_specs(const std::vector<ParamId>& p_specs) {
    std::map<ParamId, double> out;
    for (auto elem : p_specs) {            
        out[elem] = pspp->get_param(elem)->get_val();
    }
    return out;
}
std::map<ParamId, std::map<ParamId, double>> StatisticManager::get_all_correlations() {
    std::map<ParamId, std::map<ParamId, double>> res;
    CovarianceTransformer ct = CovarianceTransformer(pscp, pspp);
    res = ct.transform(this->cache.eta_specs_real);
    return res;
}

std::map<ExperimentObs, std::map<ExperimentObs, double>> StatisticManager::get_all_obs_correlations() {
    std::map<ExperimentObs, std::map<ExperimentObs, double>> res;
    CovarianceTransformer ct = CovarianceTransformer(pscp, pspp);

    res = ct.transform(this->cache.exp_obs);
    return res;
}

std::map<ExperimentObs, double> StatisticManager::get_obs_exp() {
    std::map<ExperimentObs, double> out;

    for (const auto& obsId : obs_int->get_obs_ids()) {
        auto exp_params = pspp->get_obs_param(obsId);

        for (const auto& [exp_obs, param] : exp_params) {
            if (!accepts_experiment_observable(exp_obs)) {
                continue;
            }

            out[exp_obs] = param->get_val();
        }
    }

    return out;
}

void StatisticManager::prepare_likelihood_for_scan(const std::vector<ParamId>& p_specs) {
    update_cache(p_specs);

    if (cache.p_specs.empty()) {
        throw std::invalid_argument("prepare_likelihood_for_scan called with an empty fit parameter list.");
    }

    auto unzipped_fit_params = unzip(cache.p_specs);
    auto unzipped_nuisances  = unzip(cache.eta_specs_real);
    auto unzipped_exp_obs    = unzip(cache.exp_obs);

    const std::vector<ParamId> p_ids = unzipped_fit_params.ids;
    const std::vector<double> p0 = fit_coordinates_from_model_values(
        config, p_ids, unzipped_fit_params.vals
    );

    const std::vector<ParamId> eta_ids = unzipped_nuisances.ids;
    const std::vector<double> eta0     = unzipped_nuisances.vals;

    const std::vector<ExperimentObs> obs_ids = unzipped_exp_obs.ids;
    const std::vector<double> exp_obs_vals   = unzipped_exp_obs.vals;

    auto ctx = std::make_shared<LikelihoodContext>();
    ctx->nuisance_dist = build_nuisance_distribution();
    ctx->exp_obs_dist  = build_exp_data_distribution();
    ctx->exp_obs_values = exp_obs_vals;

    ctx->fp_defs.reserve(p_ids.size());
    for (std::size_t i = 0; i < p_ids.size(); ++i) {
        double sigma = std::abs(pspp->get_param(p_ids[i])->get_combined_std().real());
        ctx->fp_defs.emplace_back(make_fit_param_def(p_ids[i], p0[i], sigma, config));
    }

    ctx->nuis_defs.reserve(eta_ids.size());
    for (std::size_t i = 0; i < eta_ids.size(); ++i) {
        double sigma = std::abs(pspp->get_param(eta_ids[i])->get_combined_std().real());
        ctx->nuis_defs.emplace_back(
            make_nuisance_parameter_definition(eta_ids[i], eta0[i], sigma)
        );
    }

    auto model_fn = [this, obs_ids, p_ids, eta_ids](
        const std::vector<double>& p_vec,
        const std::vector<double>& eta_vec) -> std::vector<double>
    {
        const auto model_p = model_values_from_fit_coordinates(this->config, p_ids, p_vec);
        auto pred_map = this->obs_int->predict_optimized(
            zip(p_ids, model_p),
            zip(eta_ids, eta_vec)
        );

        return ordered_prediction_vector(obs_ids, pred_map);
    };

    last_ctx_ = ctx;
    last_like_ = std::make_shared<BaseLikelihood>(model_fn, ctx, p_ids.size());
    last_fitter_ = std::make_shared<MLFitter>(ctx, model_fn);

    last_fit_param_ids_ = p_ids;
    last_nuisance_ids_ = eta_ids;
    last_fit_param_index_.clear();
    for (std::size_t i = 0; i < p_ids.size(); ++i) {
        last_fit_param_index_[p_ids[i]] = i;
    }

    last_scan_p_ = p0;
    last_scan_eta_ = eta0;
    has_manual_scan_point_ = false;

    if (config.print_scan_summary || config.print_debug) {
        std::cout << "[SCAN] Likelihood prepared without MLE.\n";
        std::cout << "[SCAN] n_fit_params = " << p_ids.size()
                  << ", n_nuisances = " << eta_ids.size() << "\n";
    }
}

void StatisticManager::set_manual_scan_point(const std::map<ParamId, double>& p_hat,
                                             const std::map<ParamId, double>& eta_hat) {
    if (!last_like_) {
        throw std::runtime_error(
            "Please call prepare_likelihood_for_scan(...) or compute_MLE(...) before set_manual_scan_point(...)."
        );
    }

    last_scan_p_.resize(last_fit_param_ids_.size());
    for (std::size_t i = 0; i < last_fit_param_ids_.size(); ++i) {
        auto it = p_hat.find(last_fit_param_ids_[i]);
        if (it == p_hat.end()) {
            throw std::invalid_argument("Missing manual fit-parameter value for one scanned parameter.");
        }
        last_scan_p_[i] = it->second;
    }

    last_scan_eta_.resize(last_nuisance_ids_.size());
    for (std::size_t i = 0; i < last_nuisance_ids_.size(); ++i) {
        auto it = eta_hat.find(last_nuisance_ids_[i]);
        if (it == eta_hat.end()) {
            throw std::invalid_argument("Missing manual nuisance value for one nuisance parameter.");
        }
        last_scan_eta_[i] = it->second;
    }

    has_manual_scan_point_ = true;
    if (config.print_scan_summary || config.print_debug) {
        std::cout << "[SCAN] Manual scan point loaded.\n";
    }
}

LikelihoodScanGrid StatisticManager::scan_likelihood_around_current_point(
    ParamId p1,
    ParamId p2,
    double x_half_width,
    double y_half_width,
    std::size_t nx,
    std::size_t ny
) const {
    if (!last_like_) {
        throw std::runtime_error(
            "Please call prepare_likelihood_for_scan(...) or compute_MLE(...) before requesting a likelihood scan."
        );
    }

    if (!last_fit_param_index_.contains(p1) || !last_fit_param_index_.contains(p2)) {
        throw std::invalid_argument(
            "Likelihood scan requested for parameters that are not in the current prepared parameter set."
        );
    }

    if (p1 == p2) {
        throw std::invalid_argument("Likelihood scan requires two distinct parameters.");
    }

    const std::size_t ix = last_fit_param_index_.at(p1);
    const std::size_t iy = last_fit_param_index_.at(p2);

    std::vector<double> p_ref;
    std::vector<double> eta_ref;

    if (has_manual_scan_point_) {
        p_ref = last_scan_p_;
        eta_ref = last_scan_eta_;
    } else if (!last_fit_raw_.p_hat.empty() && !last_fit_raw_.eta_hat.empty()) {
        p_ref = last_fit_raw_.p_hat;
        eta_ref = last_fit_raw_.eta_hat;
    } else if (!last_scan_p_.empty() && last_scan_eta_.size() == last_nuisance_ids_.size()) {
        // prepare_likelihood_for_scan(...) stores the current central values in
        // last_scan_p_/last_scan_eta_.  This lets lightweight scan examples run
        // without first performing a full MLE. A later compute_MLE(...) or
        // set_manual_scan_point(...) still overrides this default reference point.
        p_ref = last_scan_p_;
        eta_ref = last_scan_eta_;
    } else {
        throw std::runtime_error(
            "No reference point available. Use compute_MLE(...), "
            "set_manual_scan_point(...), or prepare_likelihood_for_scan(...)."
        );
    }

    if (ix >= p_ref.size() || iy >= p_ref.size()) {
        throw std::runtime_error("Internal error: parameter index out of range.");
    }

    std::vector<double> theta0 = p_ref;
    theta0.insert(theta0.end(), eta_ref.begin(), eta_ref.end());

    const double nll0 = last_like_->nll(theta0);

    LikelihoodScanGrid out;
    out.x_param = p1;
    out.y_param = p2;
    out.x_center = p_ref[ix];
    out.y_center = p_ref[iy];
    out.nx = nx;
    out.ny = ny;
    out.points.reserve(nx * ny);

    double nll_min = std::numeric_limits<double>::infinity();

    const double x_min = out.x_center - x_half_width;
    const double x_max = out.x_center + x_half_width;
    const double y_min = out.y_center - y_half_width;
    const double y_max = out.y_center + y_half_width;

    for (std::size_t i = 0; i < nx; ++i) {
        const double x = x_min + (x_max - x_min) * static_cast<double>(i) / static_cast<double>(nx - 1);

        for (std::size_t j = 0; j < ny; ++j) {
            const double y = y_min + (y_max - y_min) * static_cast<double>(j) / static_cast<double>(ny - 1);

            std::vector<double> theta = theta0;
            theta[ix] = x;
            theta[iy] = y;

            const double nll = last_like_->nll(theta);

            LikelihoodScanPoint pt;
            pt.x = x;
            pt.y = y;
            pt.nll = nll;

            nll_min = std::min(nll_min, nll);
            out.points.push_back(pt);
        }
    }

    for (auto& pt : out.points) {
        pt.delta_nll = pt.nll - nll_min;
    }

    if (config.print_scan_summary || config.print_debug) {
        std::cout << "[SCAN] Built likelihood scan around current point\n";
        std::cout << "[SCAN] center x = " << out.x_center << "\n";
        std::cout << "[SCAN] center y = " << out.y_center << "\n";
        std::cout << "[SCAN] nll at reference point = " << nll0 << "\n";
        std::cout << "[SCAN] min nll on grid        = " << nll_min << "\n";
    }

    return out;
}

void StatisticManager::save_likelihood_scan_csv(const std::string& path,
                                                const LikelihoodScanGrid& grid) const {
    std::ofstream out(path);
    out << "# x=" << grid.x_param << "\n";
    out << "# y=" << grid.y_param << "\n";
    out << "# x_center=" << std::setprecision(17) << grid.x_center << "\n";
    out << "# y_center=" << std::setprecision(17) << grid.y_center << "\n";
    out << "# nx=" << grid.nx << "\n";
    out << "# ny=" << grid.ny << "\n";
    out << "x,y,nll,delta_nll\n";

    out << std::setprecision(17);
    for (const auto& pt : grid.points) {
        out << pt.x << ","
            << pt.y << ","
            << pt.nll << ","
            << pt.delta_nll << "\n";
    }
}

void StatisticManager::select_experiment(const std::string& experiment) {
    select_experiments(std::set<std::string>{experiment});
}

void StatisticManager::select_experiments(const std::vector<std::string>& experiments) {
    select_experiments(std::set<std::string>(
        experiments.begin(),
        experiments.end()
    ));
}

void StatisticManager::select_experiments(const std::set<std::string>& experiments) {
    if (experiments.empty()) {
        throw std::invalid_argument(
            "StatisticManager::select_experiments: empty experiment set."
        );
    }

    selected_experiments_ = experiments;
    invalidate_fit_state();
}

void StatisticManager::select_experiments_all() {
    selected_experiments_.reset();
    invalidate_fit_state();
}

bool StatisticManager::has_experiment_selection() const noexcept {
    return selected_experiments_.has_value();
}

std::set<std::string> StatisticManager::selected_experiments() const {
    if (!selected_experiments_) {
        return {};
    }

    return *selected_experiments_;
}

bool StatisticManager::accepts_experiment_observable(const ExperimentObs& exp_obs) const {
    if (selected_experiments_.has_value()
        && !selected_experiments_->contains(exp_obs.experiment)) {
        return false;
    }

    if (selected_experiment_observables_.has_value()
        && !selected_experiment_observables_->contains(exp_obs)) {
        return false;
    }

    return true;
}

void StatisticManager::select_experiment_observables(const std::vector<ExperimentObs>& observables) {
    select_experiment_observables(std::set<ExperimentObs>(
        observables.begin(),
        observables.end()
    ));
}

void StatisticManager::select_experiment_observables(const std::set<ExperimentObs>& observables) {
    if (observables.empty()) {
        throw std::invalid_argument(
            "StatisticManager::select_experiment_observables: empty observable set."
        );
    }

    selected_experiment_observables_ = observables;
    invalidate_fit_state();
}

void StatisticManager::select_experiment_observables_all() {
    selected_experiment_observables_.reset();
    invalidate_fit_state();
}

bool StatisticManager::has_experiment_observable_selection() const noexcept {
    return selected_experiment_observables_.has_value();
}

std::set<ExperimentObs> StatisticManager::selected_experiment_observables() const {
    if (!selected_experiment_observables_) {
        return {};
    }

    return *selected_experiment_observables_;
}
