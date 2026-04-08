// #include "StatisticManager.h"

// namespace {

// std::string param_name(const ParamId& pid) {
//     std::ostringstream oss;
//     oss << pid;
//     return oss.str();
// }

// fit_app::ParameterDefinition make_fit_param_def(const ParamId& pid, double value, double sigma_hint) {
//     fit_app::ParameterDefinition out;
//     out.name = param_name(pid);
//     out.value = value;
//     out.step_hint = (std::isfinite(sigma_hint) && sigma_hint > 0.0)
//         ? sigma_hint
//         : std::max(1e-3, 0.01 * std::abs(value));

//     const std::string& nm = out.name;

//     if (nm.find("FCONST") != std::string::npos) {
//         out.limits = std::make_pair(0.05, 0.35);
//     }

//     return out;
// }

// fit_app::ParameterDefinition make_nuisance_param_def(const ParamId& pid, double value, double sigma_hint) {
//     fit_app::ParameterDefinition out;
//     out.name = param_name(pid);
//     out.value = value;

//     const double s = (std::isfinite(sigma_hint) && sigma_hint > 0.0)
//         ? std::abs(sigma_hint)
//         : std::max(1e-3, 0.01 * std::abs(value));

//     out.step_hint = s;

//     const std::string& nm = out.name;

//     // TODO : Very ugly
//     if (nm.find("SMINPUTS:3") != std::string::npos) {
//         out.limits = std::make_pair(0.05, 0.30);
//     } else if (nm.find("MASS:") != std::string::npos ||
//                nm.find("FLIFE:") != std::string::npos ||
//                nm.find("FCONST:") != std::string::npos ||
//                nm.find("FMASS:") != std::string::npos ||
//                nm.find("SMINPUTS:5") != std::string::npos ||
//                nm.find("SMINPUTS:6") != std::string::npos) {
//         out.limits = std::make_pair(std::max(1e-12, value - 5.0 * s), value + 5.0 * s);
//     }

//     return out;
// }

// template <class PredMapT>
// std::vector<double> ordered_prediction_vector(
//     const std::vector<ExperimentObs>& obs_ids,
//     const PredMapT& pred_map)
// {
//     std::vector<double> out;
//     out.reserve(obs_ids.size());

//     for (const auto& bid : obs_ids) {
//         const auto& vec = pred_map.at(bid.obs.s);

//         auto it = std::find_if(vec.begin(), vec.end(), [&](const auto& ov) {
//             auto bin = ov.bin.value_or(std::pair<double, double>{0.0, 0.0});
//             return bin == bid.obs.p;
//         });

//         if (it == vec.end()) {
//             throw std::runtime_error("Missing predicted observable/bin.");
//         }

//         out.push_back(it->value);
//     }

//     return out;
// }

// }

// StatisticManager::StatisticManager(StatisticConfig config,
//                                    std::shared_ptr<IModel> obs_int,
//                                    std::shared_ptr<IStatCorrelationProxy> pscp,
//                                    std::shared_ptr<IStatParameterProxy> pspp,
//                                    std::shared_ptr<IStatSourcesProxy> sp,
//                                    std::shared_ptr<IStatDependencyPruner> dp,
//                                    std::shared_ptr<INuisanceReader> nuisance_reader)
//     : obs_int(std::move(obs_int)),
//       pscp(std::move(pscp)),
//       pspp(std::move(pspp)),
//       sp(std::move(sp)),
//       dp(std::move(dp)),
//       nuisance_reader_(std::move(nuisance_reader)),
//       config(std::move(config))
// {
//     if (!nuisance_reader_) {
//         throw std::invalid_argument("StatisticManager: nuisance_reader is null");
//     }

//     this->obs_int->compute_observables();
//     reload_nuisance_specs();
//     invalidate_fit_state();
// }

// // std::vector<std::unique_ptr<IMarginalDistribution>> StatisticManager::build_nuisance_marginal_distributions() {
// //     // unsigned int seed = std::random_device{}(); TODO : investigate how random seed leads to unstable results
// //     unsigned int seed = 123456u;
// //     std::map<ParamId, MarginalType> nuisance_marginals;

// //     for (auto& [pid, v] : cache.eta_specs_real) {
// //         nuisance_marginals.emplace(pid, MarginalType::GAUSSIAN);
// //     }

// //     for (auto& [pid, mt] : config.override_nuisance_marginals) {
// //         for (auto elem : nuisance_marginals) {
// //             std::cout << elem.first << std::endl;
// //             std::cout << ParameterTypeMapper::str(elem.first.type.value_or(ParameterType::WILSON)) << std::endl;
// //         }
// //         if (!nuisance_marginals.contains(pid)) {
// //             LOG_WARN("Parameter", pid, "is not a nuisance for the selected observables");
// //             continue;
// //         }
// //         nuisance_marginals.at(pid) = mt;
// //     }

// //     auto unzipped = unzip(nuisance_marginals);
// //     std::vector<ParamId> nuisance_ids = unzipped.ids;
// //     std::vector<std::unique_ptr<IMarginalDistribution>> marginals;

// //     for (auto& [pid, mt] : nuisance_marginals) {
// //         MarginalConfig cfg = MarginalConfigFactory().create(pid, mt);
// //         auto m_ptr = DistributionFactory::create(mt, cfg, seed);
// //         marginals.emplace_back(std::move(m_ptr));
// //     }

// //     return marginals;
// // }

// std::vector<std::unique_ptr<IMarginalDistribution>> StatisticManager::build_nuisance_marginal_distributions() {
//     unsigned int seed = 123456u;
//     std::vector<std::unique_ptr<IMarginalDistribution>> marginals;

//     marginals.reserve(cache.eta_specs_real.size());

//     for (const auto& [pid, _] : cache.eta_specs_real) {
//         const MarginalType mt = resolve_nuisance_marginal_type(pid);
//         MarginalConfig cfg = make_nuisance_marginal_config(pid, mt);
//         marginals.emplace_back(DistributionFactory::create(mt, cfg, seed));
//     }

//     return marginals;
// }

// std::unique_ptr<JointDistribution> StatisticManager::build_nuisance_distribution() {
//     // unsigned int seed = std::random_device{}(); TODO : investigate how random seed leads to unstable results
//     unsigned int seed = 123456u;

//     std::unique_ptr<ICopula> copula;
//     if (config.nuisance_copula_type == CopulaType::GAUSSIAN) {
//         GaussianCopulaConfig copula_cfg;
//         copula_cfg.R = RealMatrix(unzip(cache.SigmaEta).vals);
//         copula = CopulaFactory::create(config.nuisance_copula_type, copula_cfg, seed);
//     } else if (config.nuisance_copula_type == CopulaType::STUDENT_T) {
//         StudentTCopulaConfig copula_cfg;
//         copula_cfg.R = RealMatrix(unzip(cache.SigmaEta).vals);
//         copula_cfg.nu = cache.eta_specs_real.size();
//         copula = CopulaFactory::create(config.nuisance_copula_type, copula_cfg, seed);
//     }

//     return std::make_unique<JointDistribution>(
//         std::move(build_nuisance_marginal_distributions()), 
//         std::move(copula)
//     );
// }

// std::unique_ptr<JointDistribution> StatisticManager::build_exp_data_distribution() {
//     unsigned int seed = std::random_device{}();
//     std::map<ExperimentObs, MarginalType> exp_data_marginals;

//     for (auto& [oid, v] : cache.exp_obs) {
//         exp_data_marginals.emplace(oid, MarginalType::GAUSSIAN);
//     }

//     for (auto& [oid, mt] : config.override_exp_data_marginals) {
//         if (!exp_data_marginals.contains(oid)) {
//             LOG_WARN("Observable", ObservableMapper::str(oid.obs.s), "is not a taken into account in this fit.");
//             continue;
//         }
//         exp_data_marginals.at(oid) = mt;
//     }

//     auto unzipped = unzip(exp_data_marginals);
//     std::vector<ExperimentObs> obs_ids = unzipped.ids;
//     std::vector<std::unique_ptr<IMarginalDistribution>> marginals;

//     for (auto& [oid, mt] : exp_data_marginals) {//TODO : checkkkkkk
//         MarginalConfig cfg = MarginalConfigFactory().create(oid, mt);
//         auto m_ptr = DistributionFactory::create(mt, cfg, seed);
//         marginals.emplace_back(std::move(m_ptr));
//     }

//     std::unique_ptr<ICopula> copula;
//     if (config.exp_data_copula_type == CopulaType::GAUSSIAN) {
//         GaussianCopulaConfig copula_cfg;
//         copula_cfg.R = RealMatrix(unzip(cache.SigmaObs).vals);
//         copula = CopulaFactory::create(config.exp_data_copula_type, copula_cfg, seed);
//     } else if (config.exp_data_copula_type == CopulaType::STUDENT_T) {
//         StudentTCopulaConfig copula_cfg;
//         copula_cfg.R = RealMatrix(unzip(cache.SigmaObs).vals);
//         copula_cfg.nu = obs_int->n_observables() - 1;
//         copula = CopulaFactory::create(config.exp_data_copula_type, copula_cfg, seed);
//     }

//     return std::make_unique<JointDistribution>(std::move(marginals), std::move(copula));
// }

// FitResultWithMaps StatisticManager::compute_MLE(const std::vector<ParamId>& p_specs) {
//     update_cache(p_specs);

//     if (cache.p_specs.empty()) {
//         throw std::invalid_argument("compute_MLE called with an empty fit parameter list.");
//     }

//     auto unzipped_fit_params = unzip(cache.p_specs);
//     auto unzipped_nuisances  = unzip(cache.eta_specs_real);
//     auto unzipped_exp_obs    = unzip(cache.exp_obs);

//     const std::vector<ParamId> p_ids = unzipped_fit_params.ids;
//     const std::vector<double> p0 = unzipped_fit_params.vals;

//     const std::vector<ParamId> eta_ids = unzipped_nuisances.ids;
//     const std::vector<double> eta0 = unzipped_nuisances.vals;

//     const std::vector<ExperimentObs> obs_ids = unzipped_exp_obs.ids;
//     const std::vector<double> exp_obs_vals = unzipped_exp_obs.vals;

//     auto ctx = std::make_shared<LikelihoodContext>();
//     ctx->nuisance_dist = build_nuisance_distribution();
//     ctx->exp_obs_dist = build_exp_data_distribution();
//     ctx->exp_obs_values = exp_obs_vals;

//     // ctx->fp_defs.reserve(p_ids.size());
//     // for (std::size_t i = 0; i < p_ids.size(); ++i) {
//     //     double sigma = std::abs(pspp->get_param(p_ids[i])->get_combined_std().real());
//     //     ctx->fp_defs.emplace_back(make_param_def(p_ids[i], p0[i], sigma));
//     // }

//     // ctx->nuis_defs.reserve(eta_ids.size());
//     // for (std::size_t i = 0; i < eta_ids.size(); ++i) {
//     //     double sigma = std::abs(pspp->get_param(eta_ids[i])->get_combined_std().real());
//     //     ctx->nuis_defs.emplace_back(make_param_def(eta_ids[i], eta0[i], sigma));
//     // }

//     ctx->fp_defs.reserve(p_ids.size());
//     for (std::size_t i = 0; i < p_ids.size(); ++i) {
//         double sigma = std::abs(pspp->get_param(p_ids[i])->get_combined_std().real());
//         ctx->fp_defs.emplace_back(make_fit_param_def(p_ids[i], p0[i], sigma));
//     }

//     // ctx->nuis_defs.reserve(eta_ids.size());
//     // for (std::size_t i = 0; i < eta_ids.size(); ++i) {
//     //     double sigma = std::abs(pspp->get_param(eta_ids[i])->get_combined_std().real());
//     //     ctx->nuis_defs.emplace_back(make_nuisance_param_def(eta_ids[i], eta0[i], sigma));
//     // }
//     ctx->nuis_defs.reserve(eta_ids.size());
//     for (std::size_t i = 0; i < eta_ids.size(); ++i) {
//         double sigma = std::abs(pspp->get_param(eta_ids[i])->get_combined_std().real());
//         ctx->nuis_defs.emplace_back(
//             make_nuisance_parameter_definition(eta_ids[i], eta0[i], sigma)
//         );
//     }

//     auto model_fn = [this, obs_ids, p_ids, eta_ids](
//         const std::vector<double>& p_vec,
//         const std::vector<double>& eta_vec) -> std::vector<double>
//     {
//         auto pred_map = this->obs_int->predict_optimized(
//             zip(p_ids, p_vec),
//             zip(eta_ids, eta_vec)
//         );

//         return ordered_prediction_vector(obs_ids, pred_map);
//     };

//     last_ctx_ = ctx;
//     last_like_ = std::make_shared<BaseLikelihood>(model_fn, ctx, p_ids.size());
//     last_fitter_ = std::make_shared<MLFitter>(ctx, model_fn);

//     last_fit_raw_ = last_fitter_->maximum_likelihood_fit(p0);

//     last_fit_param_ids_ = p_ids;
//     last_nuisance_ids_ = eta_ids;
//     last_fit_param_index_.clear();
//     for (std::size_t i = 0; i < p_ids.size(); ++i) {
//         last_fit_param_index_[p_ids[i]] = i;
//     }

//     FitResultWithMaps out;
//     out.fit_ok = !last_fit_raw_.p_hat.empty();
//     out.ell_hat = last_fit_raw_.ell_hat;
//     out.p_hat = zip(p_ids, last_fit_raw_.p_hat);
//     out.eta_hat = zip(eta_ids, last_fit_raw_.eta_hat);
//     out.p_hat_std = zip(p_ids, last_fit_raw_.p_hat_std);
//     out.p_correlations = zip(p_ids, last_fit_raw_.p_hat_correlations);

//     cache.mle_result = out;
//     return out;
// }

// Contour StatisticManager::confidence_contour(ParamId p1, ParamId p2, double z, std::array<double, 4> bounds, ContourOptions options) {
//     if (!cache.mle_result.fit_ok || !last_fitter_) {
//         throw std::runtime_error("Please run compute_MLE before requesting a confidence contour.");
//     }

//     if (!last_fit_param_index_.contains(p1) || !last_fit_param_index_.contains(p2)) {
//         throw std::invalid_argument("Contour requested for parameters that are not in the last fitted parameter set.");
//     }

//     if (p1 == p2) {
//         throw std::invalid_argument("Contour requires two distinct parameters.");
//     }

//     const std::size_t x_id = last_fit_param_index_.at(p1);
//     const std::size_t y_id = last_fit_param_index_.at(p2);

//     Contour cl = last_fitter_->contour(
//         x_id,
//         y_id,
//         z,
//         bounds,
//         options
//     );

//     return cl;
// }

// void StatisticManager::print_cache()
// {
//     for (auto elem : cache.eta_specs_real) {
//         std::cout << " eta_specs_real : " << elem.first << " = " << elem.second;
//     }
//     std::cout << std::endl;

//     for (auto elem : cache.SigmaEta) {
//         for (auto elem2 : elem.second) {
//             std::cout << " SigmaEta : " << elem.first << " | " << elem2.first << " = " << elem2.second;
//         }
//         std::cout << std::endl;
//     }

//     for (auto elem : cache.exp_obs) {
//         std::cout << " exp_obs : " << elem.first.str() << " = " << elem.second;
//     }
//     std::cout << std::endl;

//     for (auto elem : cache.SigmaObs) {
//         for (auto elem2 : elem.second) {
//             std::cout << " SigmaObs : " << elem.first.str() << " | " << elem2.first.str() << " = " << elem2.second;
//         }
//         std::cout << std::endl;
//     }

//     for (auto elem : cache.p_specs) {
//         std::cout << " p_specs : " << elem.first << " = " << elem.second;
//     }
//     std::cout << std::endl;

//     std::cout << "eta size : " << this->cache.eta_specs_real.size() << std::endl;
//     std::cout << "etasigma size : " << this->cache.eta_specs_real.size() << " | " << this->cache.SigmaEta.begin()->second.size() << std::endl;
//     std::cout << "p_specs size : " << this->cache.p_specs.size() << std::endl;
//     std::cout << "exp_obs : " << this->cache.exp_obs.size() << std::endl;
//     std::cout << "SigmaObs size : " << this->cache.SigmaObs.size() << " | " << this->cache.SigmaObs.begin()->second.size() << std::endl;
// }

// void StatisticManager::reload_nuisance_specs() {
//     default_nuisance_specs_ = nuisance_reader_->load_default();

//     if (current_user_nuisance_file_.has_value()) {
//         user_nuisance_specs_ = nuisance_reader_->load_user(*current_user_nuisance_file_);
//     } else {
//         user_nuisance_specs_ = nuisance_reader_->load_user();
//     }

//     rebuild_merged_nuisance_specs();
//     invalidate_fit_state();
// }

// void StatisticManager::set_nuisance_user_file(const fs::path& user_yaml_path) {
//     current_user_nuisance_file_ = user_yaml_path;
//     reload_nuisance_specs();
// }

// void StatisticManager::clear_nuisance_user_file() {
//     current_user_nuisance_file_.reset();
//     reload_nuisance_specs();
// }

// void StatisticManager::rebuild_merged_nuisance_specs() {
//     merged_nuisance_specs_ = default_nuisance_specs_;
//     for (const auto& [pid, spec] : user_nuisance_specs_) {
//         merged_nuisance_specs_[pid] = spec;
//     }
// }

// void StatisticManager::invalidate_fit_state() {
//     cache.mle_result = FitResultWithMaps{};

//     last_ctx_.reset();
//     last_like_.reset();
//     last_fitter_.reset();

//     last_fit_param_ids_.clear();
//     last_nuisance_ids_.clear();
//     last_fit_param_index_.clear();
// }

// const NuisanceSpec* StatisticManager::find_nuisance_spec(const ParamId& pid) const {
//     auto it = merged_nuisance_specs_.find(pid);
//     if (it == merged_nuisance_specs_.end()) {
//         return nullptr;
//     }
//     return &it->second;
// }

// MarginalType StatisticManager::resolve_nuisance_marginal_type(const ParamId& pid) const {
//     MarginalType mt = MarginalType::GAUSSIAN;

//     if (const auto* spec = find_nuisance_spec(pid)) {
//         mt = spec->marginal;
//     }

//     if (config.override_nuisance_marginals.contains(pid)) {
//         mt = config.override_nuisance_marginals.at(pid);
//     }

//     return mt;
// }

// fit_app::ParameterDefinition StatisticManager::make_nuisance_parameter_definition(const ParamId& pid,
//                                                                                   double value,
//                                                                                   double sigma_hint) const
// {
//     fit_app::ParameterDefinition out;
//     out.name = param_name(pid);
//     out.value = value;

//     const double s = (std::isfinite(sigma_hint) && sigma_hint > 0.0)
//         ? std::abs(sigma_hint)
//         : std::max(1e-3, 0.01 * std::abs(value));

//     out.step_hint = s;

//     if (const auto* spec = find_nuisance_spec(pid)) {
//         out.limits = spec->bounds;
//         return out;
//     }

//     const std::string& nm = out.name;
//     if (nm.find("SMINPUTS:3") != std::string::npos) {
//         out.limits = std::make_pair(0.05, 0.30);
//     } else if (nm.find("MASS:") != std::string::npos ||
//                nm.find("FLIFE:") != std::string::npos ||
//                nm.find("FCONST:") != std::string::npos ||
//                nm.find("FMASS:") != std::string::npos ||
//                nm.find("SMINPUTS:5") != std::string::npos ||
//                nm.find("SMINPUTS:6") != std::string::npos) {
//         out.limits = std::make_pair(std::max(1e-12, value - 5.0 * s), value + 5.0 * s);
//     }

//     return out;
// }

// MarginalConfig StatisticManager::make_nuisance_marginal_config(const ParamId& pid,
//                                                                MarginalType mt) const
// {
//     // if (const auto* spec = find_nuisance_spec(pid)) {
//     //     return MarginalConfigFactory().create(pid, mt, *spec);
//     // }
//     return MarginalConfigFactory().create(pid, mt);
// }

// std::map<ParamId, double> StatisticManager::get_all_obss_deps() {
//     std::unordered_set<ParamId> eta_infos;

//     for (const auto& obsId : obs_int->get_obs_ids()) {
//         for (auto paramId : obs_int->get_obs_deps(obsId.s)) {
//             eta_infos.insert(paramId);
//         }
//     }

//     std::unordered_set<ParamId> eta_infos_leaf = this->sp->get_all_leaf_sources(eta_infos);
//     std::map<ParamId, double> eta_specs_real_leaf;
//     std::map<ParamId, double> delta_rel;

//     for (auto& paramId : eta_infos_leaf) {
//         double u = pspp->get_param(paramId)->get_combined_std().real();
//         double val = pspp->get_param(paramId)->get_val().real();
//         if (fpeq(val, 0.0) && !fpeq(u, 0.0)) {
//             eta_specs_real_leaf[paramId] = val;
//             continue;
//         }
            
//         if (!std::isfinite(u) || fpeq(std::abs(u), 0.0))
//             delta_rel[paramId] = 0;
//         else
//             delta_rel[paramId] = std::abs(u / val);
//     }
        
    
//     using T = std::pair<ParamId, double>;
//     double delta_rel_max = std::max_element(
//         delta_rel.begin(), delta_rel.end(), 
//         [] (const T& p, const T& q) { return p.second < q.second; }
//     )->second;

//     for (auto& paramId : eta_infos_leaf) {
//         LOG_INFO("Compared to max relative uncertainty", paramId, delta_rel[paramId] / delta_rel_max);
//     }

//     for (auto& paramId : eta_infos_leaf) {
//         if (delta_rel[paramId] / delta_rel_max > config.nuisance_relevance_cutoff)
//             eta_specs_real_leaf[paramId] = pspp->get_param(paramId)->get_val();
//     }

//     LOG_INFO("Significant nuisances");
//     for (const auto& [pid, val] : eta_specs_real_leaf) {
//         LOG_INFO(pid, val);
//     }

//     return eta_specs_real_leaf;
// }

// std::map<ParamId, double> StatisticManager::get_p_specs(const std::vector<ParamId>& p_specs) {
//     std::map<ParamId, double> out;
//     for (auto elem : p_specs) {            
//         out[elem] = pspp->get_param(elem)->get_val();
//     }
//     return out;
// }
// std::map<ParamId, std::map<ParamId, double>> StatisticManager::get_all_correlations() {
//     std::map<ParamId, std::map<ParamId, double>> res;
//     CovarianceTransformer ct = CovarianceTransformer(pscp, pspp);
//     res = ct.transform(this->cache.eta_specs_real);
//     return res;
// }

// std::map<ExperimentObs, std::map<ExperimentObs, double>> StatisticManager::get_all_obs_correlations() {
//     std::map<ExperimentObs, std::map<ExperimentObs, double>> res;
//     CovarianceTransformer ct = CovarianceTransformer(pscp, pspp);

//     res = ct.transform(this->cache.exp_obs);
//     return res;
// }

// std::map<ExperimentObs, double> StatisticManager::get_obs_exp() {
//     std::map<ExperimentObs, double> out;

//     for (const auto& obsId : obs_int->get_obs_ids()) {
//         auto _ = pspp->get_obs_param(obsId);
//         for (auto _2 : _) {
//             out[_2.first] = _2.second->get_val(); 
//         }
//         // out[obsId] = pspp->get_obs_param(obsId)->get_val();
//     }
//     return out;
// }

#include "StatisticManager.h"

#include <limits>

namespace {

std::string param_name(const ParamId& pid) {
    std::ostringstream oss;
    oss << pid;
    return oss.str();
}

fit_app::ParameterDefinition make_fit_param_def(const ParamId& pid, double value, double sigma_hint) {
    fit_app::ParameterDefinition out;
    out.name = param_name(pid);
    out.value = value;
    out.step_hint = (std::isfinite(sigma_hint) && sigma_hint > 0.0)
        ? sigma_hint
        : std::max(1e-3, 0.01 * std::abs(value));

    const std::string& nm = out.name;

    if (nm.find("FCONST") != std::string::npos) {
        out.limits = std::make_pair(0.05, 0.35);
    }

    return out;
}

fit_app::ParameterDefinition make_nuisance_param_def(const ParamId& pid, double value, double sigma_hint) {
    fit_app::ParameterDefinition out;
    out.name = param_name(pid);
    out.value = value;

    const double s = (std::isfinite(sigma_hint) && sigma_hint > 0.0)
        ? std::abs(sigma_hint)
        : std::max(1e-3, 0.01 * std::abs(value));

    out.step_hint = s;

    const std::string& nm = out.name;

    // TODO : Very ugly
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


void log_selected_nuisance_diagnostics(const std::vector<ParamId>& eta_ids,
                                      const std::map<ParamId, std::map<ParamId, double>>& sigma_eta) {
    std::cout << "[FIT] Selected nuisance parameters: " << eta_ids.size() << std::endl;
    for (const auto& pid : eta_ids) {
        std::cout << "[FIT]   nuisance " << pid << std::endl;
    }

    if (sigma_eta.empty()) {
        std::cout << "[FIT] SigmaEta is empty" << std::endl;
        return;
    }

    try {
        RealMatrix sigma = RealMatrix(unzip(sigma_eta).vals);
        std::cout << "[FIT] SigmaEta shape=" << sigma.rows() << "x" << sigma.cols()
                  << ", symmetric=" << sigma.is_symmetric() << std::endl;

        RealMatrix sym = 0.5 * (sigma + sigma.transpose());
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

        std::cout << "[FIT] SigmaEta eig_min=" << min_eig
                  << ", eig_max=" << max_eig
                  << ", non_positive=" << non_pos
                  << ", tiny(|eig|<1e-12)=" << tiny
                  << std::endl;
    } catch (const std::exception& e) {
        std::cout << "[FIT] SigmaEta diagnostics failed: " << e.what() << std::endl;
    }
}

}

StatisticManager::StatisticManager(StatisticConfig config,
                                   std::shared_ptr<IModel> obs_int,
                                   std::shared_ptr<IStatCorrelationProxy> pscp,
                                   std::shared_ptr<IStatParameterProxy> pspp,
                                   std::shared_ptr<IStatSourcesProxy> sp,
                                   std::shared_ptr<IStatDependencyPruner> dp,
                                   std::shared_ptr<INuisanceReader> nuisance_reader)
    : obs_int(std::move(obs_int)),
      pscp(std::move(pscp)),
      pspp(std::move(pspp)),
      sp(std::move(sp)),
      dp(std::move(dp)),
      nuisance_reader_(std::move(nuisance_reader)),
      config(std::move(config))
{
    if (!nuisance_reader_) {
        throw std::invalid_argument("StatisticManager: nuisance_reader is null");
    }

    this->obs_int->compute_observables();
    reload_nuisance_specs();
    invalidate_fit_state();
}

// std::vector<std::unique_ptr<IMarginalDistribution>> StatisticManager::build_nuisance_marginal_distributions() {
//     // unsigned int seed = std::random_device{}(); TODO : investigate how random seed leads to unstable results
//     unsigned int seed = 123456u;
//     std::map<ParamId, MarginalType> nuisance_marginals;

//     for (auto& [pid, v] : cache.eta_specs_real) {
//         nuisance_marginals.emplace(pid, MarginalType::GAUSSIAN);
//     }

//     for (auto& [pid, mt] : config.override_nuisance_marginals) {
//         for (auto elem : nuisance_marginals) {
//             std::cout << elem.first << std::endl;
//             std::cout << ParameterTypeMapper::str(elem.first.type.value_or(ParameterType::WILSON)) << std::endl;
//         }
//         if (!nuisance_marginals.contains(pid)) {
//             LOG_WARN("Parameter", pid, "is not a nuisance for the selected observables");
//             continue;
//         }
//         nuisance_marginals.at(pid) = mt;
//     }

//     auto unzipped = unzip(nuisance_marginals);
//     std::vector<ParamId> nuisance_ids = unzipped.ids;
//     std::vector<std::unique_ptr<IMarginalDistribution>> marginals;

//     for (auto& [pid, mt] : nuisance_marginals) {
//         MarginalConfig cfg = MarginalConfigFactory().create(pid, mt);
//         auto m_ptr = DistributionFactory::create(mt, cfg, seed);
//         marginals.emplace_back(std::move(m_ptr));
//     }

//     return marginals;
// }

std::vector<std::unique_ptr<IMarginalDistribution>> StatisticManager::build_nuisance_marginal_distributions() {
    unsigned int seed = 123456u;
    std::vector<std::unique_ptr<IMarginalDistribution>> marginals;

    marginals.reserve(cache.eta_specs_real.size());

    for (const auto& [pid, _] : cache.eta_specs_real) {
        const MarginalType mt = resolve_nuisance_marginal_type(pid);
        MarginalConfig cfg = make_nuisance_marginal_config(pid, mt);
        marginals.emplace_back(DistributionFactory::create(mt, cfg, seed));
    }

    return marginals;
}

std::unique_ptr<JointDistribution> StatisticManager::build_nuisance_distribution() {
    // unsigned int seed = std::random_device{}(); TODO : investigate how random seed leads to unstable results
    unsigned int seed = 123456u;

    std::unique_ptr<ICopula> copula;
    if (config.nuisance_copula_type == CopulaType::GAUSSIAN) {
        GaussianCopulaConfig copula_cfg;
        copula_cfg.R = RealMatrix(unzip(cache.SigmaEta).vals);
        copula = CopulaFactory::create(config.nuisance_copula_type, copula_cfg, seed);
    } else if (config.nuisance_copula_type == CopulaType::STUDENT_T) {
        StudentTCopulaConfig copula_cfg;
        copula_cfg.R = RealMatrix(unzip(cache.SigmaEta).vals);
        copula_cfg.nu = cache.eta_specs_real.size();
        copula = CopulaFactory::create(config.nuisance_copula_type, copula_cfg, seed);
    }

    return std::make_unique<JointDistribution>(
        std::move(build_nuisance_marginal_distributions()), 
        std::move(copula)
    );
}

std::unique_ptr<JointDistribution> StatisticManager::build_exp_data_distribution() {
    unsigned int seed = std::random_device{}();
    std::map<ExperimentObs, MarginalType> exp_data_marginals;

    for (auto& [oid, v] : cache.exp_obs) {
        exp_data_marginals.emplace(oid, MarginalType::GAUSSIAN);
    }

    for (auto& [oid, mt] : config.override_exp_data_marginals) {
        if (!exp_data_marginals.contains(oid)) {
            LOG_WARN("Observable", ObservableMapper::str(oid.obs.s), "is not a taken into account in this fit.");
            continue;
        }
        exp_data_marginals.at(oid) = mt;
    }

    auto unzipped = unzip(exp_data_marginals);
    std::vector<ExperimentObs> obs_ids = unzipped.ids;
    std::vector<std::unique_ptr<IMarginalDistribution>> marginals;

    for (auto& [oid, mt] : exp_data_marginals) {//TODO : checkkkkkk
        MarginalConfig cfg = MarginalConfigFactory().create(oid, mt);
        auto m_ptr = DistributionFactory::create(mt, cfg, seed);
        marginals.emplace_back(std::move(m_ptr));
    }

    std::unique_ptr<ICopula> copula;
    if (config.exp_data_copula_type == CopulaType::GAUSSIAN) {
        GaussianCopulaConfig copula_cfg;
        copula_cfg.R = RealMatrix(unzip(cache.SigmaObs).vals);
        copula = CopulaFactory::create(config.exp_data_copula_type, copula_cfg, seed);
    } else if (config.exp_data_copula_type == CopulaType::STUDENT_T) {
        StudentTCopulaConfig copula_cfg;
        copula_cfg.R = RealMatrix(unzip(cache.SigmaObs).vals);
        copula_cfg.nu = obs_int->n_observables() - 1;
        copula = CopulaFactory::create(config.exp_data_copula_type, copula_cfg, seed);
    }

    return std::make_unique<JointDistribution>(std::move(marginals), std::move(copula));
}

FitResultWithMaps StatisticManager::compute_MLE(const std::vector<ParamId>& p_specs) {
    update_cache(p_specs);

    if (cache.p_specs.empty()) {
        throw std::invalid_argument("compute_MLE called with an empty fit parameter list.");
    }

    auto unzipped_fit_params = unzip(cache.p_specs);
    auto unzipped_nuisances  = unzip(cache.eta_specs_real);
    auto unzipped_exp_obs    = unzip(cache.exp_obs);

    const std::vector<ParamId> p_ids = unzipped_fit_params.ids;
    const std::vector<double> p0 = unzipped_fit_params.vals;

    const std::vector<ParamId> eta_ids = unzipped_nuisances.ids;
    const std::vector<double> eta0 = unzipped_nuisances.vals;

    log_selected_nuisance_diagnostics(eta_ids, cache.SigmaEta);

    const std::vector<ExperimentObs> obs_ids = unzipped_exp_obs.ids;
    const std::vector<double> exp_obs_vals = unzipped_exp_obs.vals;

    auto ctx = std::make_shared<LikelihoodContext>();
    ctx->nuisance_dist = build_nuisance_distribution();
    ctx->exp_obs_dist = build_exp_data_distribution();
    ctx->exp_obs_values = exp_obs_vals;

    // ctx->fp_defs.reserve(p_ids.size());
    // for (std::size_t i = 0; i < p_ids.size(); ++i) {
    //     double sigma = std::abs(pspp->get_param(p_ids[i])->get_combined_std().real());
    //     ctx->fp_defs.emplace_back(make_param_def(p_ids[i], p0[i], sigma));
    // }

    // ctx->nuis_defs.reserve(eta_ids.size());
    // for (std::size_t i = 0; i < eta_ids.size(); ++i) {
    //     double sigma = std::abs(pspp->get_param(eta_ids[i])->get_combined_std().real());
    //     ctx->nuis_defs.emplace_back(make_param_def(eta_ids[i], eta0[i], sigma));
    // }

    ctx->fp_defs.reserve(p_ids.size());
    for (std::size_t i = 0; i < p_ids.size(); ++i) {
        double sigma = std::abs(pspp->get_param(p_ids[i])->get_combined_std().real());
        ctx->fp_defs.emplace_back(make_fit_param_def(p_ids[i], p0[i], sigma));
    }

    // ctx->nuis_defs.reserve(eta_ids.size());
    // for (std::size_t i = 0; i < eta_ids.size(); ++i) {
    //     double sigma = std::abs(pspp->get_param(eta_ids[i])->get_combined_std().real());
    //     ctx->nuis_defs.emplace_back(make_nuisance_param_def(eta_ids[i], eta0[i], sigma));
    // }
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
        auto pred_map = this->obs_int->predict_optimized(
            zip(p_ids, p_vec),
            zip(eta_ids, eta_vec)
        );

        return ordered_prediction_vector(obs_ids, pred_map);
    };

    last_ctx_ = ctx;
    last_like_ = std::make_shared<BaseLikelihood>(model_fn, ctx, p_ids.size());
    last_fitter_ = std::make_shared<MLFitter>(ctx, model_fn);

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
    for (auto elem : cache.eta_specs_real) {
        std::cout << " eta_specs_real : " << elem.first << " = " << elem.second;
    }
    std::cout << std::endl;

    for (auto elem : cache.SigmaEta) {
        for (auto elem2 : elem.second) {
            std::cout << " SigmaEta : " << elem.first << " | " << elem2.first << " = " << elem2.second;
        }
        std::cout << std::endl;
    }

    for (auto elem : cache.exp_obs) {
        std::cout << " exp_obs : " << elem.first.str() << " = " << elem.second;
    }
    std::cout << std::endl;

    for (auto elem : cache.SigmaObs) {
        for (auto elem2 : elem.second) {
            std::cout << " SigmaObs : " << elem.first.str() << " | " << elem2.first.str() << " = " << elem2.second;
        }
        std::cout << std::endl;
    }

    for (auto elem : cache.p_specs) {
        std::cout << " p_specs : " << elem.first << " = " << elem.second;
    }
    std::cout << std::endl;

    std::cout << "eta size : " << this->cache.eta_specs_real.size() << std::endl;
    std::cout << "etasigma size : " << this->cache.eta_specs_real.size() << " | " << this->cache.SigmaEta.begin()->second.size() << std::endl;
    std::cout << "p_specs size : " << this->cache.p_specs.size() << std::endl;
    std::cout << "exp_obs : " << this->cache.exp_obs.size() << std::endl;
    std::cout << "SigmaObs size : " << this->cache.SigmaObs.size() << " | " << this->cache.SigmaObs.begin()->second.size() << std::endl;
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
    auto it = merged_nuisance_specs_.find(pid);
    if (it == merged_nuisance_specs_.end()) {
        return nullptr;
    }
    return &it->second;
}

MarginalType StatisticManager::resolve_nuisance_marginal_type(const ParamId& pid) const {
    MarginalType mt = MarginalType::GAUSSIAN;

    if (const auto* spec = find_nuisance_spec(pid)) {
        mt = spec->marginal;
    }

    if (config.override_nuisance_marginals.contains(pid)) {
        mt = config.override_nuisance_marginals.at(pid);
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
    // if (const auto* spec = find_nuisance_spec(pid)) {
    //     return MarginalConfigFactory().create(pid, mt, *spec);
    // }
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
    std::map<ParamId, double> eta_specs_real_leaf;
    std::map<ParamId, double> delta_rel;

    for (auto& paramId : eta_infos_leaf) {
        double u = pspp->get_param(paramId)->get_combined_std().real();
        double val = pspp->get_param(paramId)->get_val().real();
        if (fpeq(val, 0.0) && !fpeq(u, 0.0)) {
            eta_specs_real_leaf[paramId] = val;
            continue;
        }
            
        if (!std::isfinite(u) || fpeq(std::abs(u), 0.0))
            delta_rel[paramId] = 0;
        else
            delta_rel[paramId] = std::abs(u / val);
    }
        
    
    using T = std::pair<ParamId, double>;
    double delta_rel_max = std::max_element(
        delta_rel.begin(), delta_rel.end(), 
        [] (const T& p, const T& q) { return p.second < q.second; }
    )->second;

    for (auto& paramId : eta_infos_leaf) {
        LOG_INFO("Compared to max relative uncertainty", paramId, delta_rel[paramId] / delta_rel_max);
    }

    for (auto& paramId : eta_infos_leaf) {
        if (delta_rel[paramId] / delta_rel_max > config.nuisance_relevance_cutoff)
            eta_specs_real_leaf[paramId] = pspp->get_param(paramId)->get_val();
    }

    LOG_INFO("Significant nuisances");
    for (const auto& [pid, val] : eta_specs_real_leaf) {
        LOG_INFO(pid, val);
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
        auto _ = pspp->get_obs_param(obsId);
        for (auto _2 : _) {
            out[_2.first] = _2.second->get_val(); 
        }
        // out[obsId] = pspp->get_obs_param(obsId)->get_val();
    }
    return out;
}